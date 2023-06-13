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

scalar Foam::diff
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
    scalar rhoNew  = rho  + x*dRho;
    vector rhoUNew = rhoU + x*dRhoU;
    return (Gamma-1.0)*(dRhoE - (dRhoU&rhoUNew)/rhoNew + 0.5*dRho*magSqr(rhoUNew)/sqr(rhoNew));
}

scalar Foam::solveForPressure
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& dRho,
    const vector& dRhoU,
    const scalar& dRhoE,
    const scalar& pMin,
    const scalar& xMax
)
{
    label iter = 1;
    scalar x0 = xMax;
    scalar x = x0 - (func(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, x0) - pMin)/
                     diff(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, x0);
    scalar fx = func(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, x) - pMin;
    while((mag(fx) > 1e-5) && (iter < 10))
    {
        x0 = x;
        x = x0 - fx/diff(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, x0);
        fx = func(rho, rhoU, rhoE, dRho, dRhoU, dRhoE, x) - pMin;
        iter++;
    }
    return x;
}

void Foam::evaluateEigenMatrix
(
    Mat5X5& L,
    Mat5X5& R,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& normal
)
{
    scalar rRho = 1.0/rho;
    vector U = rhoU*rRho;
    scalar ek = 0.5*magSqr(U);
    scalar p = (rhoE-rho*ek)*(Gamma-1);
    scalar a = sqrt(p*Gamma*rRho);
    scalar h = (rhoE+p)*rRho;
    scalar Vn = U&normal;
    // build Mat5x5 R
    vector an = a*normal;
    R.col(0) << 1, U.x()-an.x(), U.y()-an.y(), U.z()-an.z(), h-a*Vn;
    R.col(1) << 1, U.x()       , U.y()       , U.z()       , ek;
    R.col(2) << 1, U.x()+an.x(), U.y()+an.y(), U.z()+an.z(), h+a*Vn;
    // build Mat5x5 L
    scalar b1 = (Gamma-1.0)/sqr(a), b2 = b1*ek, ra = 1.0/a;
    vector n_a = normal*ra, b1xU = b1*U;
    L.row(0) << (b2+Vn*ra), -(b1xU.x()+n_a.x()), -(b1xU.y()+n_a.y()), -(b1xU.z()+n_a.z()), b1;
    L.row(1) << -b2+1     ,   b1xU.x()         ,   b1xU.y()         ,   b1xU.z()     ,    -b1;
    L.row(2) << (b2-Vn*ra), -(b1xU.x()-n_a.x()), -(b1xU.y()-n_a.y()), -(b1xU.z()-n_a.z()), b1;
    L.row(0) *= 0.5;
    L.row(2) *= 0.5;
    if (mag(normal.x()) > 1e-8)
    {
        R.col(3) << 0,  normal.y(), -normal.x(), 0, U.x()*normal.y()-U.y()*normal.x();
        R.col(4) << 0, -normal.z(),  0, normal.x(), U.z()*normal.x()-U.x()*normal.z();
        scalar rnx = 1.0/normal.x();
        L.row(3) <<  (U.y()-Vn*normal.y())*rnx,  normal.y(),  (sqr(normal.y())-1.0)*rnx, normal.y()*normal.z()*rnx, 0;
        L.row(4) << -(U.z()-Vn*normal.z())*rnx, -normal.z(), -normal.y()*normal.z()*rnx, (1.0-sqr(normal.z()))*rnx, 0;
    }
    else if (mag(normal.y()) > 1e-8)
    {
        R.col(3) << 0, normal.y(), -normal.x(), 0, U.x()*normal.y()-U.y()*normal.x();
        R.col(4) << 0, 0, normal.z(), -normal.y(), U.y()*normal.z()-U.z()*normal.y();
        scalar rny = 1.0/normal.y();
        L.row(3) << -(U.x()-Vn*normal.x())*rny, (1.0-sqr(normal.x()))*rny, -normal.x(), normal.x()*normal.z()*rny, 0;
        L.row(4) <<  (U.z()-Vn*normal.z())*rny, normal.x()*normal.z()*rny,  normal.z(), (sqr(normal.z())-1.0)*rny, 0;
    }
    else
    {
        R.col(3) << 0, -normal.z(), 0, normal.x(), U.z()*normal.x()-U.x()*normal.z();
        R.col(4) << 0, 0, normal.z(), -normal.y(), U.y()*normal.z()-U.z()*normal.y();
        scalar rnz = 1.0/normal.z();
        L.row(3) <<  (U.x()-Vn*normal.x())*rnz,  (sqr(normal.x())-1.0)*rnz,  normal.x()*normal.y()*rnz,  normal.x(), 0;
        L.row(4) << -(U.y()-Vn*normal.y())*rnz, -normal.x()*normal.y()*rnz,  (1.0-sqr(normal.y()))*rnz, -normal.y(), 0;
    }
}