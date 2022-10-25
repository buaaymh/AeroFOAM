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
#include "dbnsFlux.H"

void Foam::consToPrim
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& gamma,
    vector& U,
    scalar& p,
    scalar& T
)
{
    U = rhoU/rho;
    p = (rhoE-0.5*rho*magSqr(U))*(gamma-1.0);
    T = p*gamma/rho;
}

void Foam::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho,
    const vector& U,
    const scalar& p,
    const vector& normal,
    const scalar& gamma
)
{
    rhoFlux  = rho*(U&normal);
    rhoUFlux = rhoFlux*U + p*normal;
    rhoEFlux = rhoFlux * (gamma/(gamma-1)*p/max(rho, SMALL) + 0.5*magSqr(U));
}

scalar Foam::entropyCorr
(
    const scalar& z,
    const scalar& d
)
{
    if (z > d)
      return z;
    else
      return 0.5*(z*z+d*d)/d;
}


// ************************************************************************* //
