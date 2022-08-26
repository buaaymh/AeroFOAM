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
    const scalar sigma = 2.0/3.0;
    const scalar k = 0.41;
    const scalar Cb1 = 0.1355;
    const scalar Cb2 = 0.622;
    const scalar Cv1 = 7.1;
    const scalar Cn1 = 16.0;
    const scalar Ct3 = 1.2;
    const scalar Ct4 = 0.5;
    const scalar Cw1 = Cb1/sqr(k)+(1.0+Cb2)/sigma;
    const scalar Cw2 = 0.3;
    const scalar Cw3 = 2.0;
    
    const scalar k2 = sqr(k);
    const scalar Cv1_3 = pow3(Cv1);
    const scalar Cw3_6 = pow(Cw3, 6.0);
    const scalar Cb2_sigma = Cb2/sigma;

    scalar ST_inf;
    scalar Ma_Re_inf;
    //- nuTilde
    volScalarField& nuTilde_;
    //- nuTilde grad
    volVectorField nuTildeGrad_;
    //- laminar viscosity
    scalarField laminarViscosity_;
    //- eddy viscosity
    scalarField eddyViscosity_;
    //- nuMax
    scalarField nuMax_;
    //- hybrid delta
    scalarField delta_;
    // Private Member Functions
    //- Disallow default bitwise copy construct
    turb3rdSolver(const turb3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const turb3rdSolver&);

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
        const tensor& UGrad_L,
        const tensor& UGrad_R,
        const vector& TGrad_L,
        const vector& TGrad_R,
        const scalar& nuTilde_L,
        const scalar& nuTilde_R,
        const scalar& dxTilde,
        const vector& normal
    ) const;
    
    void evaluateWallFlux
    (
        vector& rhoUFlux,
        const scalar& p,
        const scalar& T,
        const tensor& UGrad,
        const vector& normal
    ) const;

    void evaluateInletOutletFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho,
        const vector& U,
        const scalar& p,
        const scalar& nuTilde,
        const tensor& UGrad,
        const vector& TGrad,
        const vector& normal
    ) const;

    void evaluateTurbSource
    (
        scalarField& resNuTilde,
        scalarField& D
    ) const;

    scalar eddyViscosityFunc
    (
        const scalar nuTilde,
        const scalar laminarViscosity,
        const scalar rho
    ) const;

    scalar effectiveNuFunc
    (
        const scalar nu,
        const scalar nuTilde
    ) const;

    void computeProjectedGradient
    (
        const vector& normal,
        const vector& edgeVec,
        const scalar& var_i,
        const scalar& var_j,
        const vector& grad_i,
        const vector& grad_j,
        scalar& proj_vector_ij,
        scalar& correctedMeanGrad
    ) const;

    void calculateDelta();

public:

    // Constructors - default
    turb3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p,
        volScalarField& nuTilde
    );

    // Destructor - default
    virtual ~turb3rdSolver() {};

    // Member Functions
    const scalarField& nuTilde() const { return nuTilde_.primitiveField(); }
    scalarField& nuTilde() { return nuTilde_.primitiveFieldRef(); }

    void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

    void evaluateTurbRes
    (
        scalarField& resNuTilde,
        scalarField& L,
        scalarField& U,
        scalarField& D
    );

    void solveTurbLinearSystem
    (
        scalarField& resNuTilde,
        const scalarField& L,
        const scalarField& U,
        const scalarField& D
    );

    void solveTurbPseudoTimeSystem
    (
        const scalar& dt,
        const scalar& alpha,
        scalarField& resNuTilde,
        const scalarField& L,
        const scalarField& U,
        const scalarField& D
    );

    void correctFields();

    void updateLTS();

protected:
    
    void evaluateMatrixLDU
    (
        scalarField& L,
        scalarField& U,
        scalarField& D
    ) const;

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
