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
    scalar Pr_Lam  = 0.72;
    scalar Pr_Turb = 0.9;
    //- Constant
    scalar Cv1 = 7.1;

    scalar ST_inf;
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
    // Private Member Functions
    //- Disallow default bitwise copy construct
    turb3rdSolver(const turb3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const turb3rdSolver&);

    scalar eddyViscosityFunc
    (
        const scalar nuTilde,
        const scalar laminarViscosity,
        const scalar rho
    ) const;

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

    void evaluateTurbRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE,
        scalarField& resNuTilde
    );

    void solveTurbLinearSystem
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE,
        scalarField& resNuTilde
    );

    virtual void solveTurbPseudoTimeSystem
    (
        const scalar& dt,
        const scalar& alpha,
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE,
        scalarField& resNuTilde,
        scalar& resL1
    );

    void correctFields();

    void updateLTS();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
