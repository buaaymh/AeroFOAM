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
    volScalarField& p
)
:
    eulerSolver(fluidProps, rho, U, p),
    Limiter(mesh_),
    rhoGrad_
    (
        IOobject
        (
            "rhoGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    UGrad_
    (
        IOobject
        (
            "UGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor(dimless/dimLength, tensor::zero)
    ),
    pGrad_
    (
        IOobject
        (
            "pGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "evaluateFlowRes.H"

vector Foam::euler2ndSolver::sampleVelocity
(
    label cellI,
    const vector& point
) const
{
    return U_[cellI] + (UGrad_[cellI]&(point-mesh_.C()[cellI]));
}

// ************************************************************************* //
