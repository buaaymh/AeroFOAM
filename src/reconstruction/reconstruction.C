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

#include "reconstruction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::Reconstruction::Reconstruction
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& rhoU,
    volScalarField& rhoE
)
:
    fluidProps_(fluidProps),
    mesh_(rho.mesh()),
    rho_(rho),
    rhoU_(rhoU),
    rhoE_(rhoE),
    vrWeightSqr_(vector::zero),
    isTrouble_(mesh_.nCells(), false),
    coefs_(mesh_.nCells(), Mat9X5::Zero()),
    limitedCoefs_i_(mesh_.nInternalFaces(), Mat5X9::Zero()),
    limitedCoefs_j_(mesh_.nInternalFaces(), Mat5X9::Zero()),
    rDeltaXYZ_
    (
        IOobject
        (
            "rDeltaXYZ",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless, vector::zero)
    ),
    basisMean_
    (
        IOobject
        (
            "basisMean",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    IS_
    (
        IOobject
        (
            "IS",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    ),
    N_h_(mesh_.nCells(), 0.0),
    d_ij_(mag(mesh_.delta())),
    rA_(mesh_.nCells(), Mat9X9::Zero()),
    lowerb_(mesh_.nInternalFaces(), Col9X1::Zero()),
    upperb_(mesh_.nInternalFaces(), Col9X1::Zero()),
    B_(mesh_.nInternalFaces(), Mat9X9::Zero()),
    polyNorm_(mesh_.nCells(), Col9X1::Zero()),
    d1Rho_
    (
        IOobject
        (
            "d1Rho",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    d1RhoU_
    (
        IOobject
        (
            "d1RhoU",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor(dimless/dimLength, tensor::zero)
    ),
    d1RhoE_
    (
        IOobject
        (
            "d1RhoE",
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
    )
{
    initVrLinearSystem();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "functions.H"

#include "evaluateVarsAndGrads.H"

#include "initVrLinearSystem.H"

#include "iterationStep.H"

#include "limitCoefficients.H"

// ************************************************************************* //