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

#include "euler2ndSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::euler2ndSolver::euler2ndSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& T
)
:
    solver(fluidProps, rho, U, T),
    rhoLimit_
    (
        IOobject
        (
            "rhoLimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    ),
    ULimit_
    (
        IOobject
        (
            "ULimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless, vector::one)
    ),
    TLimit_
    (
        IOobject
        (
            "TLimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    ),
    rhoMin_(scalarField(mesh_.nCells())),
    rhoMax_(scalarField(mesh_.nCells())),
    UMin_(vectorField(mesh_.nCells())),
    UMax_(vectorField(mesh_.nCells())),
    TMin_(scalarField(mesh_.nCells())),
    TMax_(scalarField(mesh_.nCells()))
{
    Info << "Ths solver is 2nd order for Euler flow." << nl
         << "Ths limiter is Venkatakrishnan." << nl
         << "====================================================" << endl << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "limitGrad.H"

#include "evaluateFlowRes.H"

// ************************************************************************* //
