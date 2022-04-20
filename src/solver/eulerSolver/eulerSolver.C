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

#include "fvCFD.H"
#include "eulerSolver.H"
#include "AUSMplusUpFlux.H"
#include "hllcFlux.H"
#include "roeFlux.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::eulerSolver::eulerSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p
)
:
    fluidProps_(fluidProps),
    mesh_(rho.mesh()),
    normal_(mesh_.Sf()/mesh_.magSf()),
    rho_(rho),
    U_(U),
    p_(p),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh_.time().timeName(),
            mesh_
        ),
        rho_*U_
    ),
    rhoE_
    (
        IOobject
        (
            "rhoE",
            mesh_.time().timeName(),
            mesh_
        ),
        p_/(fluidProps_.gamma-1.0) + 0.5*rho_*magSqr(U_)
    ),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluidProps_.gamma*p_/rho_
    ),
    c_
    (
        IOobject
        (
            "c",
            mesh_.time().timeName(),
            mesh_
        ),
        sqrt(T_)
    ),
    Ma_
    (
        IOobject
        (
            "Ma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mag(U_)/Foam::sqrt(T_)
    ),
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
    ),
    rhoLimit_
    (
        IOobject
        (
            "rhoLimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
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
        dimensionedVector(dimless, vector::zero)
    ),
    pLimit_
    (
        IOobject
        (
            "pLimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),
    rhoMin_(scalarField(mesh_.nCells())),
    rhoMax_(scalarField(mesh_.nCells())),
    UMin_(vectorField(mesh_.nCells())),
    UMax_(vectorField(mesh_.nCells())),
    pMin_(scalarField(mesh_.nCells())),
    pMax_(scalarField(mesh_.nCells())),
    volProjections_(vectorField(mesh_.nCells())),
    localDtDv_(scalarField(mesh_.nCells()))
{
    word flux = mesh_.schemesDict().subDict("divSchemes").lookupOrDefault<word>("flux", "roe");
    Info << "The flux scheme is " << flux << endl;
    if (flux == "hllc")
        riemann_ = std::make_unique<hllcFlux>();
    else if (flux == "AUSMplusUp")
        riemann_ = std::make_unique<AUSMplusUpFlux>();
    else if (flux == "roe")
        riemann_ = std::make_unique<roeFlux>();
    else
    {
        Info << "Error in Riemann solver, you should choose hllc, AUSMplusUp or roe" << endl;
        riemann_ = std::make_unique<roeFlux>();
    }
    volProjectionsInit();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "limitGrad.H"

#include "evaluateFlowRes.H"

#include "correctFields.H"

#include "localTimeStep.H"

// ************************************************************************* //
