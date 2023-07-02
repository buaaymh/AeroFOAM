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

\*---------------------------------------------------------------------------*/

#include "turbulence3rdSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::turbulence3rdSolver::turbulence3rdSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p,
    volScalarField& nuTilda
)
:
    turbulenceSolver(fluidProps, rho, U, p, nuTilda),
    Reconstruction(fluidProps, rho, rhoU_, rhoE_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "evaluateFlowRes.H"

#include "functions.H"

void Foam::turbulence3rdSolver::evaluateRhoUForForce
(
    scalar& rho,
    vector& U,
    const vector& coord,
    const label& cellI
)
{
    scalar rhoE; vector rhoU;
    const vector delta = coord - mesh_.C()[cellI];
    evaluateVars
    (
        rho, rhoU, rhoE,
        rho_[cellI], rhoU_[cellI], rhoE_[cellI], coefs_[cellI],
        basisMean_[cellI], rDeltaXYZ_[cellI], delta
    );
    U = rhoU/rho;
}

// ************************************************************************* //