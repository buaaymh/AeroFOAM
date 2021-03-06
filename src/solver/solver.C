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
#include "solver.H"
#include "AUSMplusUpFlux.H"
#include "hllcFlux.H"
#include "roeFlux.H"

Foam::solver::solver
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
        p_*fluidProps_.gamma/rho_
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
        p_/(fluidProps_.gamma-1.0) + 0.5*rho_*magSqr(U_)
    ),
    c_(sqrt(T_.primitiveField())),
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
    volProjections_(vectorField(mesh_.nCells())),
    localDtDv_(scalarField(mesh_.nCells()))
{
    Info << "=================Solver Information==================" << endl;
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

void Foam::solver::correctFields()
{
    U_.ref() = rhoU_ / rho_;
    p_.ref() = (rhoE_ - 0.5*rho_*magSqr(U_)) * (fluidProps_.gamma-1.0);
    const bool rhoBool = Foam::positiveCorrect(rho_);
    const bool pBool   = Foam::positiveCorrect(p_);
    if (rhoBool || pBool)
    {
        rhoU_ = rho_ * U_;
        rhoE_ = p_/(fluidProps_.gamma-1.0) + 0.5*rho_*magSqr(U_);
    }
    T_.ref() = p_*fluidProps_.gamma/rho_;
    c_ = sqrt(T_.primitiveField());
    Ma_.primitiveFieldRef() = mag(U_.primitiveFieldRef())/c_;
    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    const volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U_.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p_.boundaryFieldRef();
    volVectorField::Boundary& rhoUBf = rhoU_.boundaryFieldRef();
    volScalarField::Boundary& rhoEBf = rhoE_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T_.boundaryFieldRef();
    volScalarField::Boundary& MaBf = Ma_.boundaryFieldRef();
    rhoUBf = rhoBf * UBf;
    rhoEBf = pBf/(fluidProps_.gamma-1.0) + 0.5*rhoBf*magSqr(UBf);
    TBf = pBf*fluidProps_.gamma/rhoBf;
    MaBf = mag(UBf)/sqrt(TBf);
}

void Foam::solver::updateLTS()
{
    scalar localCFL = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("LocalCFL", 1.0);
    localDtDv_ = localCFL/(volProjections_&(cmptMag(U_.primitiveField())+c_*vector::one));
}

void Foam::solver::volProjectionsInit()
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

#include "lusgsGMRES.H"

#include "solveFlowLinearSystem.H"

// ************************************************************************* //