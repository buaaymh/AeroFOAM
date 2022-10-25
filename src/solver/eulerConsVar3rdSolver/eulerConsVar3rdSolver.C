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

#include "eulerConsVar3rdSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::eulerConsVar3rdSolver::eulerConsVar3rdSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p
)
:
    vr3rdSolver(fluidProps, rho, U, p),
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
    rhoGradLimited_
    (
        IOobject
        (
            "rhoGradLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    rhoUGrad_
    (
        IOobject
        (
            "rhoUGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor(dimless/dimLength, tensor::zero)
    ),
    rhoUGradLimited_
    (
        IOobject
        (
            "rhoUGradLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor(dimless/dimLength, tensor::zero)
    ),
    rhoEGrad_
    (
        IOobject
        (
            "rhoEGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    rhoEGradLimited_
    (
        IOobject
        (
            "rhoEGradLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    d2Rho_
    (
        IOobject
        (
            "d2Rho",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoUx_
    (
        IOobject
        (
            "d2RhoUx",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoUy_
    (
        IOobject
        (
            "d2RhoUy",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoUz_
    (
        IOobject
        (
            "d2RhoUz",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoE_
    (
        IOobject
        (
            "d2RhoE",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoLimited_
    (
        IOobject
        (
            "d2RhoLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoUxLimited_
    (
        IOobject
        (
            "d2RhoUxLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoUyLimited_
    (
        IOobject
        (
            "d2RhoUyLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoUzLimited_
    (
        IOobject
        (
            "d2RhoUzLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2RhoELimited_
    (
        IOobject
        (
            "d2RhoELimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    )
{
    Info << "Ths solver is 3rd order for Euler flow." << nl
         << "Ths scheme is VR with AV." << nl
         << "====================================================" << endl << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
#include "reconstructionIter.H"

#include "limitCoefficients.H"

#include "evaluateFlowRes.H"

// ************************************************************************* //
