/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Reconstruction

Description
    Variation reconstruction approach.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef reconstruction_H
#define reconstruction_H

#include "element.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include <map>
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct Parameter
{
    word   limiter;
    vector vrWeightSqr;
    label  nTroubled, nNegative;
    scalar IS, np, nCells;
    bool   indicator, positive, characteristic;
};

class Reconstruction
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    Reconstruction(const Reconstruction&);
    //- Disallow default bitwise assignment
    void operator=(const Reconstruction&);

    // Member Functions
    void initVrLinearSystem();

    void limitCoefficients();

protected:
    const fluidProperties& fluidProps_;
    //- Reference to mesh
    const fvMesh& mesh_;
    //- Static density
    volScalarField& rho_;
    //- Momentum
    volVectorField& rhoU_;
    //- Energy
    volScalarField& rhoE_;
    Parameter parameter_;
    std::vector<bool> isTrouble_;
    std::vector<Mat9X5> coefs_;
    std::vector<Mat5X9> limitedCoefs_i_;
    std::vector<Mat5X9> limitedCoefs_j_;
    //- Mesh data
    volVectorField rDeltaXYZ_;
    volSymmTensorField basisMean_;
    volScalarField IS_;
    scalarField N_h_;
    scalarField d_ij_;
    //- Linear system
    std::vector<Mat9X9> rA_;
    std::vector<Col9X1> lowerb_;
    std::vector<Col9X1> upperb_;
    std::vector<Mat9X9> B_;
    std::map<label, std::vector<Mat9X9>> coupledB_;
    std::map<label, std::vector<Col9X1>> coupledb_;
    //- Quadrature points
    std::map<label, std::vector<std::unique_ptr<Face>>> patch2ndQuad_;
    std::map<label, std::vector<std::unique_ptr<Face>>> patch4stQuad_;
    std::vector<std::unique_ptr<Face>> faceQuad_;
    std::vector<Col9X1> polyNorm_;
    //- First polynomial coefficients
    volVectorField d1Rho_;
    volTensorField d1RhoU_;
    volVectorField d1RhoE_;
    //- Quadratic polynomial coefficients
    volSymmTensorField d2Rho_;
    volSymmTensorField d2RhoUx_;
    volSymmTensorField d2RhoUy_;
    volSymmTensorField d2RhoUz_;
    volSymmTensorField d2RhoE_;

public:
    Reconstruction
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& rhoU,
        volScalarField& rhoE
    );

    const Parameter& parameter() { return parameter_; }

    void iterationStep();

    // Destructor - default
    virtual ~Reconstruction() {};

};

/*---------------------------------------------------------------------------*\
                           Methods Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * rDeltaXYZ and basisMean * * * * * * * * * * * * //

void evaluateBasisFunctionConst
(
    const fvMesh& mesh,
    const label& cellI,
    vector& rDeltaXYZ,
    symmTensor& basisMean,
    Col9X1& polyNorm
);

// * * * * * * ** * * * Taylor Basis Functions * * * * * * * * * * * //
Col9X1 polynomialDn0
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
);

Col9X1 polynomialDn1
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const vector& normal
);

Col9X1 polynomialDn2
(
    const vector& rDeltaXYZ,
    const vector& normal
);

Mat9X3 polynomialGrad
(
    const vector& delta,
    const vector& rDeltaXYZ
);

// * * * * * * * * * * * * WBAP Functions * * * * * * * * * * * * //

Mat5X3 transformD1Coefs
(
    const Mat5X9& coefs,
    const vector& delta,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
);

Mat5X6 transformD2Coefs
(
    const Mat5X6& d2Coefs,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
);

template<int cols>
Eigen::Matrix<scalar, 5, cols> WBAP_L2_Limiter
(
    const std::vector<Eigen::Matrix<scalar, 5, cols>>& list,
    const scalar& value
)
{
    const label n = list.size();
    Eigen::Matrix<scalar, 5, cols> result = list[0];
    for (label i = 0; i != 5; i++)
    {
        for (label j = 0; j != cols; j++)
        {
            if (mag(list[0](i, j)) < 1e-10) continue;
            scalar temp3 = value;
            scalar temp4 = value;
            for (label k = 1; k != n; k++)
            {
                const scalar theta = list[k](i, j)/list[0](i, j);
                if (theta <= 1e-6)
                {
                    temp3 = 0;
                    break;
                }
                temp3 += 1.0 / Foam::pow3(theta);
                temp4 += 1.0 / Foam::pow4(theta);
            }
            result(i, j) *= temp3 / temp4;
        }
    }
    return result;
}

template<int cols>
Eigen::Matrix<scalar, 5, cols> Eigen_WBAP_L2_Limiter
(
    const std::vector<vector>& normalList,
    const std::vector<Eigen::Matrix<scalar, 5, cols>>& matList,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& value
)
{
    const label n = matList.size();
    std::vector<Eigen::Matrix<scalar, 5, cols>> limited(n-1);
    std::vector<Eigen::Matrix<scalar, 5, cols>> tempList(n);
    Mat5X5 L, R;
    for (label j = 0; j < n-1; j++)
    {
        evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normalList[j]);
        for (label i = 0; i < n; i++) tempList[i] = L*matList[i];
        limited[j] = R*WBAP_L2_Limiter<cols>(tempList, value);
    }
    return WBAP_L2_Limiter<cols>(limited, 1);
}

void New_WBAP_Limiter
(
    Mat5X9& coefs_i,
    Mat5X9& coefs_j,
    const vector& delta,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j,
    const scalar& value
);

Mat5X3 slope_Limiter
(
    const std::vector<vector>& normalList,
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ
);

Mat5X3 slope_Limiter
(
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ
);

// * * * * * * * * * * evaluateVar and Grad Functions * * * * * * * * * * * //

void evaluateVars
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const Mat9X5& coefs,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
);

void evaluateVars
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const vector& d1Rho,
    const tensor& d1RhoU,
    const vector& d1RhoE,
    const symmTensor& d2Rho,
    const symmTensor& d2RhoUx,
    const symmTensor& d2RhoUy,
    const symmTensor& d2RhoUz,
    const symmTensor& d2RhoE,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
);

void evaluateVarsAndGrads
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    vector& rhoGrad,
    tensor& rhoUGrad,
    vector& rhoEGrad,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const Mat9X5& coefs,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
);

void evaluateVarsAndGrads
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    vector& rhoGrad,
    tensor& rhoUGrad,
    vector& rhoEGrad,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const vector& d1Rho,
    const tensor& d1RhoU,
    const vector& d1RhoE,
    const symmTensor& d2Rho,
    const symmTensor& d2RhoUx,
    const symmTensor& d2RhoUy,
    const symmTensor& d2RhoUz,
    const symmTensor& d2RhoE,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //