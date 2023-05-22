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

#include "navierStokesSolver.H"

Foam::navierStokesSolver::navierStokesSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p
)
:
    eulerSolver(fluidProps, rho, U, p),
    Ma_Re_(fluidProps.Mach_inf/fluidProps.Re_inf),
    S_T_(110.4/fluidProps.T_inf),
    muLam_
    (
        IOobject
        (
            "laminarViscosity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
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
    nuMax_(scalarField(mesh_.nCells(), 0.0)),
    delta_(mesh_.delta())
{}

void Foam::navierStokesSolver::correctFields()
{
    conservativeToPrimitiveFields();
    muLam_ = pow(T_, 1.5)*(1.0+S_T_)/(T_+S_T_);
    nuMax_ = (max(4.0/3.0, Gamma)/Pr_Lam)*muLam_.primitiveFieldRef()/rho_.primitiveField();
}

#include "evaluateFlowRes.H"

#include "functions.H"

