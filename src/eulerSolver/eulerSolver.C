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
        p_*Gamma/rho_
    ),
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
        p_/(Gamma-1.0) + 0.5*rho_*magSqr(U_)
    ),
    c_(sqrt(T_.primitiveField())),
    volProjections_(vectorField(mesh_.nCells())),
    localDtDv_(scalarField(mesh_.nCells()))
{
    Info << "=================Solver Information==================" << endl;
    word flux = mesh_.schemes().subDict("divSchemes").lookupOrDefault<word>("flux", "roe");
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

void Foam::eulerSolver::volProjectionsInit()
{
    volProjections_ = vectorField(mesh_.nCells(), vector::zero);
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    forAll(owner, faceI)
    {
        const vector faceVector = 0.5*cmptMag(Sf[faceI]);
        volProjections_[owner[faceI]] += faceVector;
        volProjections_[neighbour[faceI]] += faceVector;
    }
    forAll(mesh_.boundary(), patchI)
    {   
        scalar temp = 0.25;
        if (mesh_.boundary()[patchI].coupled()) temp = 0.5;
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const auto& Sf = mesh_.boundary()[patchI].Sf();
        forAll(bfaceCells, faceI)
            volProjections_[bfaceCells[faceI]] += temp*cmptMag(Sf[faceI]);
    }
}

#include "correctFields.H"

#include "functions.H"

#include "evaluateFlowRes.H"

#include "solveFlowLinearSystem.H"

// ************************************************************************* //