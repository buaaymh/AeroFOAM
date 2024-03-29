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
    Omega_
    (
        IOobject
        (
            "vorticity",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    fr1_
    (
        IOobject
        (
            "rotationFunction",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    ),
    dNuTilda_(scalarField(mesh_.nCells())),
    maxDelta_(scalarField(mesh_.nCells()))
{
    dist_ = wallDist::New(mesh_).y().primitiveField();
    if (fluidProps.simulationType == "DES")
    {
        forAll(mesh_.C(), cellI)
        {
            scalar deltaMaxTmp = 0.0;
            const labelList& cFaces = mesh_.cells()[cellI];
            const point& centrevector = mesh_.C()[cellI];
            forAll(cFaces, faceI)
            {
                label facei = cFaces[faceI];
                const point& facevector = mesh_.Cf()[facei];
                scalar tmp = mag(facevector - centrevector);
                deltaMaxTmp = max(tmp, deltaMaxTmp);
            }
            maxDelta_[cellI] = deltaMaxTmp;
        }
    }
}

void Foam::turbulenceSolver::correctTurbulenceFields()
{
    nuLam_ = muLam_/rho_;
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
    UGrad_ = fvc::grad(U_);
    Omega_ = Foam::sqrt(2.0)*mag(skew(UGrad_));
    if (fluidProps_.simulationType == "SARCM")
    {
        forAll(mesh_.C(), cellI)
        {
            const scalar asymInnerProduct = magSqr(Omega_[cellI]);
            const scalar symInnerProduct = 2.0*magSqr(symm(UGrad_[cellI]));
            const scalar w = atan(0.01*asymInnerProduct)
                           * 2.0/constant::mathematical::pi
                           * (asymInnerProduct-symInnerProduct)+symInnerProduct;
            const scalar rStar  = sqrt(symInnerProduct/max(w, SMALL));
            scalar rTilda = sqrt(w/max(symInnerProduct, SMALL));
            rTilda *= rTilda-1;
            fr1_[cellI] = max
            (
                min
                (
                    (1.0+SA::Cr1)*2.0*rStar/(1.0+rStar)*
                    (1.0-SA::Cr3*atan(SA::Cr2*rTilda))-SA::Cr1,
                    10.0
                ),
                -10.0
            );
        }
    }
    /*--- 
    A New Version of Detached-eddy Simulation, Resistant to Ambiguous Grid Densities.
    Spalart et al. Theoretical and Computational Fluid Dynamics - 2006
    ---*/
    if (fluidProps_.simulationType == "DES")
    {
        dist_ = wallDist::New(mesh_).y().primitiveField();
        forAll(mesh_.C(), cellI)
        {
            scalar r_d = Ma_Re_*nuTilda_[cellI]/(max(mag(UGrad_[cellI]), 1e-10)*SA::k2*sqr(dist_[cellI]));
            scalar f_d = 1.0-Foam::tanh(Foam::pow3(8.0*r_d));
            dist_[cellI] -= f_d*max(0.0, (dist_[cellI]-maxDelta_[cellI]*SA::constDES));
        }
    }
}

#include "functions.H"

#include "evaluateFlowRes.H"

#include "evaluateTurbRes.H"

#include "solveTurbLinearSystem.H"