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

#include "model.H"

Foam::Model::Model
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    volVectorField& force
)
:
    name_(name),
    mesh_(rho.mesh()),
    rho_(rho),
    U_(U),
    force_(force)
{
    zoneI_ = mesh_.cellZones().findZoneID(name);
    word blade = mesh_.solutionDict().subDict(name).lookup<word>("blade");
    if (blade == "CaradonnaTung")
    {
        blade_ = std::make_unique<CaradonnaTung>();
    }
    else if (blade == "UH60A")
    {
        blade_ = std::make_unique<UH60A>();
    }
    else
    {
        Info << "Error in blade type" << nl
             << "(" << nl
             << " CaradonnaTung" << nl
             << " UH60A" << nl
             << ")" << nl
             << endl;
    }
}

scalar Foam::Model::getAngleOfAttack
(
    scalar u,
    scalar w,
    scalar twist
) const
{
    auto deg = rad2deg(Foam::atan(w / u));  // [-90, 90]
    if (u < 0) {
      if (w > 0) {
        deg += 180/* [-90, 0] -> [90, 180] */;
      } else {
        deg -= 180/* [0, 90] -> [-180, -90] */;
      }
    }
    // deg := angle of inflow
    deg += twist;
    // deg := angle of attack
    if (deg < -180 || 180 < deg) {
      deg += (deg < 0 ? 360 : -360);
    }
    // deg in [-180, 180]
    return deg;
}

scalar Foam::Model::get3DGaussWeight
(
    scalar d2,
    scalar eps
) const
{
    return Foam::exp(-d2/sqr(eps))/(pow3(eps)*sqrt(pow3(constant::mathematical::pi)));
}