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

Class
    dbnsFlux

Description
    Riemann solver for Euler fluxes.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef dbnsFlux_H
#define dbnsFlux_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dbnsFlux Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class dbnsFlux
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    dbnsFlux(const dbnsFlux&);
    //- Disallow default bitwise assignment
    void operator=(const dbnsFlux&);

public:

    // Constructors - default
    dbnsFlux() {};

    // Destructor - default
    virtual ~dbnsFlux() {};

    // Member Functions
    virtual void evaluateFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho_L,
        const scalar& rho_R,
        const vector& U_L,
        const vector& U_R,
        const scalar& p_L,
        const scalar& p_R,
        const vector& normal,
        const scalar& gamma
    ) const = 0;
};

/*---------------------------------------------------------------------------*\
                           Methods Declaration
\*---------------------------------------------------------------------------*/

void consToPrim
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    const scalar& gamma
);

void primToCons
(
    scalar& rho,
    vector& U,
    scalar& p,
    const scalar& gamma
);

void evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho,
    const vector& U,
    const scalar& p,
    const vector& normal,
    const scalar& gamma
);

scalar entropyCorr
(
    const scalar& z,
    const scalar& d
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
