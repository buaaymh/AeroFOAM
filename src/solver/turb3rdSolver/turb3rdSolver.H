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
    turb3rdSolver

Description
    3rd order solver for Turbulence flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef turb3rdSolver_H
#define turb3rdSolver_H

#include "euler3rdSolver.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class turb3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class turb3rdSolver: public euler3rdSolver
{
    //- Prandtl numbers
    const scalar Pr_Lam  = 0.72;
    const scalar Pr_Turb = 0.9;
    //- Constant
    // const scalar sigma = 2.0/3.0;
    // const scalar k = 0.41;
    // const scalar Cb1 = 0.1355;
    // const scalar Cb2 = 0.622;
    // const scalar Cv1 = 7.1;
    // const scalar Cn1 = 16.0;
    // const scalar Ct3 = 1.2;
    // const scalar Ct4 = 0.5;
    // const scalar Cw1 = Cb1/sqr(k)+(1.0+Cb2)/sigma;
    // const scalar Cw2 = 0.3;
    // const scalar Cw3 = 2.0;
    
    // const scalar k2 = sqr(k);
    // const scalar Cv1_3 = pow3(Cv1);
    // const scalar Cw3_6 = pow(Cw3, 6.0);
    // const scalar Cb2_sigma = Cb2/sigma;

    // Private Member Functions
    //- Disallow default bitwise copy construct
    turb3rdSolver(const turb3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const turb3rdSolver&);

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
    ) const;

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
        const vector& rhoGrad_L,
        const vector& rhoGrad_R,
        const tensor& rhoUGrad_L,
        const tensor& rhoUGrad_R,
        const vector& rhoEGrad_L,
        const vector& rhoEGrad_R,
        const scalar& deltaTilde,
        const vector& normal,
        const scalar& magSf
    ) const;

    void evaluateWallFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho,
        const vector& rhoU,
        const scalar& rhoE,
        const vector& rhoGrad,
        const tensor& rhoUGrad,
        const vector& rhoEGrad,
        const vector& normal,
        const scalar& magSf
    ) const;

    void evaluateInletOutletFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rhoBound,
        const vector& UBound,
        const scalar& pBound,
        const scalar& rho,
        const vector& rhoU,
        const scalar& rhoE,
        const vector& rhoGrad,
        const tensor& rhoUGrad,
        const vector& rhoEGrad,
        const vector& normal,
        const scalar& magSf
    ) const;
    

    // void evaluateTurbSource
    // (
    //     scalarField& resNuTilde,
    //     scalarField& D
    // ) const;

    // scalar effectiveNuFunc
    // (
    //     const scalar nu,
    //     const scalar nuTilde
    // ) const;

    // void computeProjectedGradient
    // (
    //     const vector& normal,
    //     const vector& edgeVec,
    //     const scalar& var_i,
    //     const scalar& var_j,
    //     const vector& grad_i,
    //     const vector& grad_j,
    //     scalar& proj_vector_ij,
    //     scalar& correctedMeanGrad
    // ) const;

    // void calculateDelta();

public:

    // Constructors - default
    turb3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~turb3rdSolver() {};

    // Member Functions
    // const scalarField& nuTilde() const { return nuTilde_.primitiveField(); }
    // scalarField& nuTilde() { return nuTilde_.primitiveFieldRef(); }

    void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

    // void evaluateTurbRes
    // (
    //     scalarField& resNuTilde,
    //     scalarField& L,
    //     scalarField& U,
    //     scalarField& D
    // );

    // void solveTurbLinearSystem
    // (
    //     scalarField& resNuTilde,
    //     const scalarField& L,
    //     const scalarField& U,
    //     const scalarField& D
    // );

    // void solveTurbPseudoTimeSystem
    // (
    //     const scalar& dt,
    //     const scalar& alpha,
    //     scalarField& resNuTilde,
    //     const scalarField& L,
    //     const scalarField& U,
    //     const scalarField& D
    // );

    void correctFields();

protected:
    
    void evaluateMatrixLDU
    (
        const scalar& alpha,
        scalarField& L,
        scalarField& U,
        scalarField& D
    ) const;

};

// scalar eddyViscosityFunc
// (
//     const scalar nuTilde,
//     const scalar laminarViscosity,
//     const scalar rho
// );

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
