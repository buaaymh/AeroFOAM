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

#include "euler2ndSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class euler3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class euler3rdSolver: public solver
{

public:
    using Matrix = Eigen::Matrix<scalar, 6, 6>;
    using Column = Eigen::Matrix<scalar, 6, 1>;

private:
    // Private Member Functions
    //- Disallow default bitwise copy construct
    euler3rdSolver(const euler3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const euler3rdSolver&);

    void evaluateFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho_L,
        const scalar& rho_R,
        const vector& U_L,
        const vector& U_R,
        const scalar& T_L,
        const scalar& T_R,
        const vector& normal
    ) const;

protected:
    //- Density grad
    volVectorField rhoGrad_;
    volVectorField rhoGradLimited_;
    //- Velocity grad
    volTensorField UGrad_;
    volTensorField UGradLimited_;
    //- pressure grad
    volVectorField TGrad_;
    volVectorField TGradLimited_;
    //- Quadratic polynomial coefficients
    volSymmTensorField d2Rho_;
    volSymmTensorField d2Ux_;
    volSymmTensorField d2Uy_;
    volSymmTensorField d2Uz_;
    volSymmTensorField d2T_;
    volSymmTensorField d2RhoLimited_;
    volSymmTensorField d2UxLimited_;
    volSymmTensorField d2UyLimited_;
    volSymmTensorField d2UzLimited_;
    volSymmTensorField d2TLimited_;
    volVectorField rLengthScale_;
    volSymmTensorField basisConst_;
    scalarField N_h_;
    std::vector<bool> p0_;
    std::vector<bool> trouble_;
    //- Linear system
    std::vector<Matrix> rA_;
    std::vector<Matrix> B_;
    //- Quadrature points
    std::vector<std::vector<vector>> quad_;

    void initVrLinearSystem();

    void updateCoefficients();

    void limitCoefficients();

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

// * * * * * * * * * * * *Taylor Basis Functions * * * * * * * * * * * * * * //

symmTensor basisPoly
(
    const vector& delta,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

symmTensor basisPolyDn1
(
    const vector& delta,
    const vector& rLengthScale,
    const vector& normal
);

symmTensor basisPolyDn2
(
    const vector& rLengthScale,
    const vector& normal
);

euler3rdSolver::Column basisPolyEigen
(
    const vector& delta,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

euler3rdSolver::Column basisPolyDn1Eigen
(
    const vector& delta,
    const vector& rLengthScale,
    const vector& normal
);

scalar scalarPoly
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

vector vectorPoly
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

vector scalarPolyGrad
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rLengthScale
);

tensor vectorPolyGrad
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rLengthScale
);

// * * * * * * * * * * * * * * * WBAP Functions * * * * * * * * * * * * * * //

vector vectorWBAP
(
    const std::vector<vector>& list,
    scalar value
);

symmTensor symmTensorWBAP
(
    const std::vector<symmTensor>& list,
    scalar value
);

tensor tensorWBAP
(
    const std::vector<tensor>& list,
    scalar value
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
