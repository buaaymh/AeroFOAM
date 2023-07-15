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
        p_*Gamma/rho_
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
    localDtDv_(scalarField(mesh_.nCells())),
    dRho_(scalarField(mesh_.nCells())),
    dRhoU_(vectorField(mesh_.nCells())),
    dRhoE_(scalarField(mesh_.nCells()))
{
    Info << "=================Solver Information==================" << endl;
    word flux = mesh_.schemesDict().subDict("divSchemes").lookupOrDefault<word>("flux", "roe");
    Info << "# Flux scheme            [-] = " << flux << endl;
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
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const auto& Sf = mesh_.boundary()[patchI].Sf();
        forAll(bfaceCells, faceI)
            volProjections_[bfaceCells[faceI]] += 0.5*cmptMag(Sf[faceI]);
    }
}

void Foam::solver::correctFields()
{
    U_.ref() = rhoU_.internalField()/rho_.internalField();
    p_.ref() = (rhoE_.internalField()-0.5*rho_.internalField()*magSqr(U_.internalField()))*(Gamma-1.0);
    const bool rhoBool = Foam::positiveCorrect(rho_);
    const bool pBool   = Foam::positiveCorrect(p_);
    if (rhoBool || pBool)
    {
        rhoU_.ref() = rho_*U_;
        rhoE_.ref() = p_/(Gamma-1.0)+0.5*rho_*magSqr(U_);
    }
    rho_.correctBoundaryConditions();
    rhoU_.correctBoundaryConditions();
    rhoE_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    T_.primitiveFieldRef() = Gamma*p_.primitiveField()/rho_.primitiveField();
    c_ = sqrt(T_.primitiveField());
    Ma_.primitiveFieldRef() = mag(U_.primitiveField())/c_;
    forAll(mesh_.boundary(), patchI)
    {
        const word name = mesh_.boundary()[patchI].name();
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        fvPatchScalarField& rhoBound = rho_.boundaryFieldRef()[patchI];
        fvPatchVectorField& UBound = U_.boundaryFieldRef()[patchI];
        fvPatchScalarField& pBound = p_.boundaryFieldRef()[patchI];
        if (name == "farField")
        {
            if (fluidProps_.Mach_inf < 1.0)
            {
                forAll(bfaceCells, faceI)
                {
                    const label i = bfaceCells[faceI];
                    const scalar c2 = T_[i];
                    const scalar qn = U_[i]&normal[faceI];
                    const scalar rhoc = c_[i]*rho_[i];
                    if (qn < 0.0)
                    {
                        const scalar pB = pBound[faceI];
                        pBound[faceI] = 0.5*(pB+p_[i]-rhoc*(normal[faceI]&(UBound[faceI]-U_[i])));
                        rhoBound[faceI] += (pBound[faceI]-pB)/c2;
                        UBound[faceI] += normal[faceI]*(pBound[faceI]-pB)/rhoc;
                    }
                    else
                    {
                        rhoBound[faceI] = rho_[i]+(pBound[faceI]-p_[i])/c2;
                        UBound[faceI] = U_[i]+normal[faceI]*(p_[i]-pBound[faceI])/rhoc;
                    }
                }
            }
            else
            {
                forAll(bfaceCells, faceI)
                {
                    const label i = bfaceCells[faceI];
                    const scalar qn = U_[i]&normal[faceI];
                    if (qn > 0.0)
                    {
                        rhoBound[faceI] = rho_[i];
                        UBound[faceI]   = U_[i];
                        pBound[faceI]   = p_[i];
                    }
                }
            }
        }
        if (name == "inlet")
        {
            const scalarField Ma = mag(UBound)/sqrt(pBound*Gamma/rhoBound);
            forAll(bfaceCells, faceI)
            {
                if (Ma[faceI] < 1.0)
                {
                    const label i = bfaceCells[faceI];
                    const scalar c2 = T_[i];
                    const scalar rhoc = rho_[i]*c_[i];
                    const scalar pB = pBound[faceI];
                    pBound[faceI] = 0.5*(pB+p_[i]-rhoc*(normal[faceI]&(UBound[faceI]-U_[i])));
                    rhoBound[faceI] += (pBound[faceI]-pB)/c2;
                    UBound[faceI] += normal[faceI]*(pBound[faceI]-pB)/rhoc;
                }
            }
        }
        if (name == "outlet")
        {
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                if (Ma_[i] < 1.0)
                {
                    const scalar c2 = T_[i];
                    const scalar rhoc = rho_[i]*c_[i];
                    rhoBound[faceI] = rho_[i]+(pBound[faceI]-p_[i])/c2;
                    UBound[faceI] = U_[i]+normal[faceI]*(p_[i]-pBound[faceI])/rhoc;
                    if ((UBound[faceI]&normal[faceI]) < 0.0)
                    {
                        const vector signU = vector(sign(U_[i].x()), sign(U_[i].y()), sign(U_[i].z()));
                        UBound[faceI] = cmptMultiply(signU, max(cmptMag(UBound[faceI]), cmptMag(U_[i])));
                    }
                }
            }
        }
    }
    T_.boundaryFieldRef() = p_.boundaryField()*Gamma/rho_.boundaryField();
    Ma_.boundaryFieldRef() = mag(U_.boundaryField())/Foam::sqrt(T_.boundaryField());
}

// ************************************************************************* //