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

#include "eulerPrimVar3rdSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::eulerPrimVar3rdSolver::eulerPrimVar3rdSolver
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
    UGradLimited_
    (
        IOobject
        (
            "UGradLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor(dimless/dimLength, tensor::zero)
    ),
    TGrad_
    (
        IOobject
        (
            "TGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    TGradLimited_
    (
        IOobject
        (
            "TGradLimited",
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
    d2Ux_
    (
        IOobject
        (
            "d2Ux",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2Uy_
    (
        IOobject
        (
            "d2Uy",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2Uz_
    (
        IOobject
        (
            "d2Uz",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2T_
    (
        IOobject
        (
            "d2T",
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
    d2UxLimited_
    (
        IOobject
        (
            "d2UxLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2UyLimited_
    (
        IOobject
        (
            "d2UyLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2UzLimited_
    (
        IOobject
        (
            "d2UzLimited",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2TLimited_
    (
        IOobject
        (
            "d2TLimited",
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
