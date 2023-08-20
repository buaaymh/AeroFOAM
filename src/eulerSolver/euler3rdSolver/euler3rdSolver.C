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

#include "euler3rdSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::euler3rdSolver::euler3rdSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p
)
:
    eulerSolver(fluidProps, rho, U, p),
    Reconstruction(fluidProps, rho, rhoU_, rhoE_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::euler3rdSolver::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& rhoU_L,
    const vector& rhoU_R,
    const scalar& rhoE_L,
    const scalar& rhoE_R,
    const vector& normal,
    const scalar& magSf
) const
{
    scalar rhoFluxTemp;
    vector rhoUFluxTemp;
    scalar rhoEFluxTemp;

    scalar p_L, p_R, T_L, T_R;
    vector U_L, U_R;
    consToPrim(rho_L, rhoU_L, rhoE_L, U_L, p_L, T_L);
    consToPrim(rho_R, rhoU_R, rhoE_R, U_R, p_R, T_R);
    p_L = max(p_L, SMALL);
    p_R = max(p_R, SMALL);
    riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal);
    rhoFlux  += magSf * rhoFluxTemp;
    rhoUFlux += magSf * rhoUFluxTemp;
    rhoEFlux += magSf * rhoEFluxTemp;
}

vector Foam::euler3rdSolver::sampleVelocity
(
    label cellI,
    const vector& point
) const
{
    const vector delta = point-mesh_.C()[cellI];
    scalar rho, rhoE; vector rhoU;
    evaluateVars(rho, rhoU, rhoE, rho_[cellI], rhoU_[cellI], rhoE_[cellI],
                 coefs_[cellI], basisMean_[cellI], rDeltaXYZ_[cellI], delta);
    return rhoU/rho;
}

#include "evaluateFlowRes.H"

// ************************************************************************* //
