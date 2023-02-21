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
    euler3rdSolver

Description
    3rd order solver for Euler flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef euler3rdSolver_H
#define euler3rdSolver_H

#include "solver.H"
#include "element.H"
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class euler3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class euler3rdSolver: public solver
{

private:
    // Private Member Functions
    //- Disallow default bitwise copy construct
    euler3rdSolver(const euler3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const euler3rdSolver&);
    
    // Member Functions
    void initVrLinearSystem();

    void limitCoefficients();

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
    ) const;

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
    ) const;

    void evaluateFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho_L,
        const scalar& rho_R,
        const vector& rhoU_L,
        const vector& rhoU_R,
        const scalar& rhoE_L,
        const scalar& rhoE_R,
        const vector& normal,
        const scalar& magSf
    ) const;

protected:
    vector vrWeightSqr_;
    std::vector<Mat9X5> coefs_;
    //- Mesh data
    volVectorField rDeltaXYZ_;
    volSymmTensorField basisMean_;
    volScalarField IS_;
    scalarField N_h_;
    vectorField delta_;
    scalarField d_ij_;
    //- Linear system
    std::vector<Mat9X9> rA_;
    std::vector<Col9X1> lowerb_;
    std::vector<Col9X1> upperb_;
    std::vector<Mat9X9> B_;
    std::map<label, std::vector<Mat9X9>> coupledB_;
    std::map<label, std::vector<Col9X1>> coupledb_;
    //- Quadrature points
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
    
    void reconstructionIter();

public:

    // Constructors - default
    euler3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~euler3rdSolver() {};

    // Member Functions
    void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

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
    symmTensor& basisMean
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

// * * * * * * * * * * * * WBAP Functions * * * * * * * * * * * * //

Mat3X5 transformD1Coefs
(
    const Mat9X5& coefs,
    const vector& delta,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
);

Mat6X5 transformD2Coefs
(
    const Mat6X5& d2Coefs,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
);

void evaluateEigenMatrix
(
    Mat5X5& L,
    Mat5X5& R,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& normal,
    const scalar& gamma
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
                temp3 += 1.0 / Foam::sqr(theta);
                temp4 += 1.0 / Foam::pow3(theta);
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
    const scalar& gamma,
    const scalar& value
)
{
    const label n = matList.size();
    Eigen::Matrix<scalar, 5, cols> limited(matList[0]);
    std::vector<Eigen::Matrix<scalar, 5, cols>> tempList(n);
    Mat5X5 L, R;
    for (const auto& normal : normalList)
    {
        evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal, gamma);
        for (label i = 0; i < n; i++)
        {
            tempList[i] = L*matList[i];
        }
        Eigen::Matrix<scalar, 5, cols> temp = R*WBAP_L2_Limiter<cols>(tempList, value);
        for (label i = 0; i != 5; i++)
        {
            for (label j = 0; j != cols; j++)
            {
                limited(i,j) = MinMod(limited(i,j), temp(i,j));
            }
        }
    }
    return limited;
}

Mat5X3 slope1_Limiter
(
    const std::vector<vector>& normalList,
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ,
    const scalar& gamma
);

Mat5X3 slope1_Limiter
(
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
