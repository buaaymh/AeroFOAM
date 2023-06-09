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
#include "turbulenceSolver.H"

Foam::turbulenceSolver::turbulenceSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p,
    volScalarField& nuTilda
)
:
    navierStokesSolver(fluidProps, rho, U, p),
    nuTilda_(nuTilda),
    nuLam_
    (
        IOobject
        (
            "nuLam",
            mesh_.time().timeName(),
            mesh_
        ),
        muLam_/rho_ 
    ),
    muTurb_
    (
        IOobject
        (
            "turbulenceViscosity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),
    nuTildaGrad_
    (
        IOobject
        (
            "nuTildaGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    dNuTilda_(scalarField(mesh_.nCells()))
{}

void Foam::turbulenceSolver::correctTurbulenceFields()
{
    nuTilda_.correctBoundaryConditions();
    forAll(mesh_.boundary(), patchI)
    {
        const word name = mesh_.boundary()[patchI].name();
        if (name == "farField")
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const scalarField qn = U_.boundaryField()[patchI]&normal_.boundaryField()[patchI];
            fvPatchScalarField& nuTildaBound = nuTilda_.boundaryFieldRef()[patchI];
            forAll(bfaceCells, faceI)
            {
                if (qn[faceI] > 0.0) nuTildaBound[faceI] = nuTilda_[bfaceCells[faceI]];
            }
        }
    }
    muTurb_ = muTurb(nuTilda_, rho_, muLam_);
}

void Foam::turbulenceSolver::correctFields()
{
    navierStokesSolver::correctFields();
    correctTurbulenceFields();
}

void Foam::turbulenceSolver::smoothFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE,
    scalarField& resNuTilda,
    const scalar& eps
)
{
    volScalarField resRhoStar
    (
        IOobject
        (
            "resRhoStar",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    );
    volVectorField resRhoUStar
    (
        IOobject
        (
            "resRhoUStar",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless, vector::zero)
    );
    volScalarField resRhoEStar
    (
        IOobject
        (
            "resRhoEStar",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    );
    volScalarField resNuTildaStar
    (
        IOobject
        (
            "resNuTildaStar",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    );
    
    resRhoStar.primitiveFieldRef()     = resRho;
    resRhoUStar.primitiveFieldRef()    = resRhoU;
    resRhoEStar.primitiveFieldRef()    = resRhoE;
    resNuTildaStar.primitiveFieldRef() = resNuTilda;

    scalarField bRho(mesh_.nCells(), 0.0);
    vectorField bRhoU(mesh_.nCells(), vector::zero);
    scalarField bRhoE(mesh_.nCells(), 0.0);
    scalarField bNuTilda(mesh_.nCells(), 0.0);
    scalarField coefOfA(mesh_.nCells(), 1.0);

    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(mesh_.owner(), faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        coefOfA[i]  += eps;
        coefOfA[j]  += eps;
        bRho[i]     += resRhoStar[j];
        bRhoU[i]    += resRhoUStar[j];
        bRhoE[i]    += resRhoEStar[j];
        bNuTilda[i] += resNuTildaStar[j];
        bRho[j]     += resRhoStar[i];
        bRhoU[j]    += resRhoUStar[i];
        bRhoE[j]    += resRhoEStar[i];
        bNuTilda[j] += resNuTildaStar[i];
    }
    resRhoStar.correctBoundaryConditions();
    resRhoUStar.correctBoundaryConditions();
    resRhoEStar.correctBoundaryConditions();
    resNuTildaStar.correctBoundaryConditions();
    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const scalarField resRhoStar_neigh  = resRhoStar.boundaryField()[patchI].patchNeighbourField();
            const vectorField resRhoUStar_neigh = resRhoUStar.boundaryField()[patchI].patchNeighbourField();
            const scalarField resRhoEStar_neigh = resRhoEStar.boundaryField()[patchI].patchNeighbourField();
            const scalarField resNuTildaStar_neigh = resNuTildaStar.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                coefOfA[i]  += eps;
                bRho[i]     += resRhoStar_neigh[faceI];
                bRhoU[i]    += resRhoUStar_neigh[faceI];
                bRhoE[i]    += resRhoEStar_neigh[faceI];
                bNuTilda[i] += resNuTildaStar_neigh[faceI];
            }
        }
    }
    resRhoStar.primitiveFieldRef()     = eps*(resRho  + eps*bRho) /coefOfA;
    resRhoUStar.primitiveFieldRef()    = eps*(resRhoU + eps*bRhoU)/coefOfA;
    resRhoEStar.primitiveFieldRef()    = eps*(resRhoE + eps*bRhoE)/coefOfA;
    resNuTildaStar.primitiveFieldRef() = eps*(resNuTilda + eps*bNuTilda)/coefOfA;

    forAll(mesh_.owner(), faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        resRho[i]     += resRhoStar[j];
        resRhoU[i]    += resRhoUStar[j];
        resRhoE[i]    += resRhoEStar[j];
        resNuTilda[i] += resNuTildaStar[j];
        resRho[j]     += resRhoStar[i];
        resRhoU[j]    += resRhoUStar[i];
        resRhoE[j]    += resRhoEStar[i];
        resNuTilda[j] += resNuTildaStar[i];
    }
    resRhoStar.correctBoundaryConditions();
    resRhoUStar.correctBoundaryConditions();
    resRhoEStar.correctBoundaryConditions();
    resNuTildaStar.correctBoundaryConditions();
    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const scalarField resRhoStar_neigh  = resRhoStar.boundaryField()[patchI].patchNeighbourField();
            const vectorField resRhoUStar_neigh = resRhoUStar.boundaryField()[patchI].patchNeighbourField();
            const scalarField resRhoEStar_neigh = resRhoEStar.boundaryField()[patchI].patchNeighbourField();
            const scalarField resNuTildaStar_neigh = resNuTildaStar.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                resRho[i]  +=  resRhoStar_neigh[faceI];
                resRhoU[i] +=  resRhoUStar_neigh[faceI];
                resRhoE[i] +=  resRhoEStar_neigh[faceI];
                resNuTilda[i] += resNuTildaStar_neigh[faceI];
            }
        }
    }
    resRho     /= coefOfA;
    resRhoU    /= coefOfA;
    resRhoE    /= coefOfA;
    resNuTilda /= coefOfA;
}

#include "functions.H"

#include "evaluateFlowRes.H"

#include "evaluateTurbRes.H"

#include "solveTurbLinearSystem.H"