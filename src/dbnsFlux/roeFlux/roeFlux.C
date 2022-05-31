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
#include "roeFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::roeFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& U_L,
    const vector& U_R,
    const scalar& T_L,
    const scalar& T_R,
    const vector& normal,
    const scalar& gamma
) const
{   
    const scalar ggm1  = gamma/(gamma-1);
    // left & right state
    const scalar p_L = rho_L*T_L/gamma;
    const scalar p_R = rho_R*T_R/gamma;
    const scalar H_L = ggm1*p_L/max(rho_L, SMALL) + 0.5*magSqr(U_L);
    const scalar H_R = ggm1*p_R/max(rho_R, SMALL) + 0.5*magSqr(U_R);
    const scalar qsRho_L = rho_L*(U_L&normal);
    const scalar qsRho_R = rho_R*(U_R&normal);
    const scalar p_A = 0.5 * (p_L + p_R);
    rhoFlux  = 0.5 * (qsRho_L     + qsRho_R    );
    rhoUFlux = 0.5 * (qsRho_L*U_L + qsRho_R*U_R) + p_A*normal;
    rhoEFlux = 0.5 * (qsRho_L*H_L + qsRho_R*H_R);

    // Roe's average
    const scalar rho_A = sqrt(max(rho_L*rho_R, SMALL));
    const scalar dd    = rho_A/max(rho_L, SMALL);
    const scalar dd1   = 1.0/(1.0+dd);
    const vector U_A   = (U_L+dd*U_R)*dd1;
    const scalar H_A   = (H_L+dd*H_R)*dd1;
    const scalar U2_A  = 0.5*magSqr(U_A);
    const scalar c2_A  = (gamma-1)*(H_A-U2_A);
    const scalar c_A   = sqrt(c2_A);
    const scalar Vn_A  = U_A & normal;
    const scalar dU    = (U_R-U_L)&normal;

    // eigenvalues
    scalar h1 = mag(Vn_A - c_A);
    scalar h2 = mag(Vn_A);
    scalar h4 = mag(Vn_A + c_A);
    const scalar delta = 0.1 * c_A;

    const scalar eabs1 = Foam::entropyCorr(h1, delta);
    const scalar eabs2 = Foam::entropyCorr(h2, delta);
    const scalar eabs4 = Foam::entropyCorr(h4, delta);

    // upwind fluxes
    const scalar temp = 1.0/c2_A;
    h1              = rho_A*c_A*dU;
    h2              = eabs1*(p_R-p_L     - h1)*0.5*temp;
    const scalar h3 = eabs2*(rho_R-rho_L - (p_R-p_L)*temp);
    h4              = eabs2*rho_A;
    const scalar h5 = eabs4*(p_R-p_L     + h1)*0.5*temp;

    rhoFlux  -= 0.5*(h2                  + h3                                     + h5);
    rhoUFlux -= 0.5*(h2*(U_A-c_A*normal) + h3*U_A  + h4*(U_R-U_L-dU*normal)       + h5*(U_A+c_A*normal));
    rhoEFlux -= 0.5*(h2*(H_A-c_A*Vn_A)   + h3*U2_A + h4*((U_A&(U_R-U_L))-Vn_A*dU) + h5*(H_A+c_A*Vn_A));
}


// ************************************************************************* //
