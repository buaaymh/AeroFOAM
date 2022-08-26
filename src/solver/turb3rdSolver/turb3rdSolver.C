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

#include "turb3rdSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::turb3rdSolver::turb3rdSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p,
    volScalarField& nuTilde
)
:
    euler3rdSolver(fluidProps, rho, U, p),
    ST_inf(110.4/fluidProps.T_inf),
    Ma_Re_inf(fluidProps.Mach_inf/fluidProps.Re_inf),
    nuTilde_(nuTilde),
    nuTildeGrad_
    (
        IOobject
        (
            "nuTildeGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    laminarViscosity_(scalarField(mesh_.nCells())),
    eddyViscosity_(scalarField(mesh_.nCells())),
    nuMax_(scalarField(mesh_.nCells())),
    delta_(scalarField(mesh_.nCells()))
{
    Info << "Ths solver is 3rd order for Turbulence flow." << nl
         << "Ths scheme is VR with AV." << nl
         << "====================================================" << endl << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turb3rdSolver::correctFields()
{
    U_.ref() = rhoU_ / rho_;
    p_.ref() = (rhoE_ - 0.5*rho_*magSqr(U_))*(fluidProps_.gamma-1.0);
    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    nuTilde_.correctBoundaryConditions();
    const bool rhoBool = Foam::positiveCorrect(rho_);
    const bool pBool   = Foam::positiveCorrect(p_);
    if (rhoBool || pBool)
    {
        rhoU_ = rho_ * U_;
        rhoE_ = p_/(fluidProps_.gamma-1.0) + 0.5*rho_*magSqr(U_);
    }
    T_.ref() = p_ * fluidProps_.gamma / rho_;
    c_ = sqrt(T_.primitiveField());
    Ma_.primitiveFieldRef() = mag(U_.primitiveFieldRef())/c_;
    laminarViscosity_ = pow(T_.primitiveField(), 1.5)*(1.0+ST_inf)/(T_.primitiveField()+ST_inf);
    forAll(mesh_.cells(), cellI)
        eddyViscosity_[cellI] = eddyViscosityFunc(nuTilde_[cellI], laminarViscosity_[cellI], rho_[cellI]);
    nuMax_ = max(4.0/3.0, fluidProps_.gamma)*(laminarViscosity_/Pr_Lam + eddyViscosity_/Pr_Turb)/rho_.primitiveField();
    const volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U_.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p_.boundaryFieldRef();
    volVectorField::Boundary& rhoUBf = rhoU_.boundaryFieldRef();
    volScalarField::Boundary& rhoEBf = rhoE_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T_.boundaryFieldRef();
    volScalarField::Boundary& MaBf = Ma_.boundaryFieldRef();
    rhoUBf = rhoBf * UBf;
    rhoEBf = pBf/(fluidProps_.gamma-1.0) + 0.5*rhoBf*magSqr(UBf);
    TBf = pBf * fluidProps_.gamma / rhoBf;
    MaBf = mag(UBf)/sqrt(TBf);
}

void Foam::turb3rdSolver::updateLTS()
{
    scalar localCFL = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("LocalCFL", 1.0);
    scalarField lambdaAc = volProjections_&(cmptMag(U_.primitiveField())+c_*vector::one);
    scalarField lambdaAv = magSqr(volProjections_)*nuMax_/mesh_.V();
    localDtDv_ = localCFL/(lambdaAc+lambdaAv);
}

#include "functions.H"

#include "evaluateFlux.H"

#include "evaluateFlowRes.H"

#include "evaluateTurbRes.H"

#include "evaluateTurbSource.H"

#include "evaluateMatrixLDU.H"

#include "solveTurbLinearSystem.H"

// ************************************************************************* //
