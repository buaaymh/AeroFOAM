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

#include "definitions.H"

void Foam::givensRotation
(
    const scalar& h,
    const scalar& beta,
    scalar& c,
    scalar& s
)
{
    if (beta == 0)
    {
        c = 1;
        s = 0;
    }
    else if (mag(beta) > mag(h))
    {
        scalar tau = -h/beta;
        s = 1.0/Foam::sqrt(1.0 + sqr(tau));
        c = s*tau;
    }
    else
    {
        scalar tau = -beta/h;
        c = 1.0/Foam::sqrt(1.0 + sqr(tau));
        s = c*tau;
    }
}

scalar Foam::MinMod
(
    const scalar& a,
    const scalar& b
)
{
    if (a*b <= 0) return 0.0;
    else
    {
        if (mag(a) < mag(b)) return a;
        else return b;
    }
}

scalar Foam::Venkat
(
    const scalar& project,
    const scalar& deltaMin,
    const scalar& deltaMax,
    const scalar& eps2
)
{
    if (project > SMALL)
    {
        scalar y = deltaMax*(deltaMax+project) + eps2;
        return (y + deltaMax*project) / (y + 2*project*project);
    }
    else if (project < -SMALL)
    {
        scalar y = deltaMin*(deltaMin+project) + eps2;
        return (y + deltaMin*project) / (y + 2*project*project);
    }
    return 1.0;
}

scalar Foam::Barth
(
    const scalar& project,
    const scalar& deltaMin,
    const scalar& deltaMax
)
{
    if (project > deltaMax)
    {
        return deltaMax / project;
    }
    else if (project < deltaMin)
    {
        return deltaMin / project;
    }
    return 1.0;
}

vector Foam::vorticity(const tensor& UGrad)
{
    return vector(UGrad.zy()-UGrad.yz(),
                  UGrad.xz()-UGrad.zx(),
                  UGrad.yx()-UGrad.xy());
}

scalar Foam::func
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& dRho,
    const vector& dRhoU,
    const scalar& dRhoE,
    const scalar& x
)
{
    return (Gamma-1.0)*(rhoE+x*dRhoE - 0.5*magSqr(rhoU+x*dRhoU)/(rho+x*dRho));
}

scalar Foam::solveForPressure
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& dRho,
    const vector& dRhoU,
    const scalar& dRhoE,
    const scalar& xMax
)
{
    label iter = 1;
    scalar a = 0.0, b = xMax, c = 0.5*(a+b);
    scalar fa = func(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, a) - 1e-3;
    scalar fc = func(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, c) - 1e-3;
    while (mag(fc) > 1e-5)
    {
        iter++;
        if (fa*fc > 0) a = c;
        else b = c;
        c = 0.5*(a+b);
        fc = func(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, c) - 1e-3;
        if (iter > 10)
        {
            if (fc < 0) c = 0.0;
            break;
        }
    }
    return c;
}