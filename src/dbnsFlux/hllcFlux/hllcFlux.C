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
#include "hllcFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::hllcFlux::evaluateFlux
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
    const scalar temp  = 1.0/(gamma-1);
    // left & right state
    scalar rRho = 1.0/rho_L;
    const scalar p_L = T_L*rho_L/gamma;
    const vector rhoU_L = rho_L*U_L;
    const scalar rhoE_L = p_L*temp + 0.5*rho_L*magSqr(U_L);
    const scalar H_L = (rhoE_L + p_L)*rRho;
    const scalar a_L = Foam::sqrt(max(0.0, T_L));
    rRho = 1.0/rho_R;
    const scalar p_R = T_R*rho_R/gamma;
    const vector rhoU_R = rho_R*U_R;
    const scalar rhoE_R = p_R*temp + 0.5*rho_R*magSqr(U_R);
    const scalar H_R = (rhoE_R + p_R)*rRho;
    const scalar a_R = Foam::sqrt(max(0.0, T_R));

    // Compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft  = (U_L & normal);
    const scalar qRight = (U_R & normal);

    // Step 2:
    // Compute Roe weights
    const scalar rho_A = sqrt(max(rho_L*rho_R, SMALL));
    const scalar dd    = rho_A/max(rho_L, SMALL);
    const scalar dd1   = 1.0/(1.0+dd);
    const vector UTilde = (U_L+dd*U_R)*dd1;
    const scalar HTilde = (H_L+dd*H_R)*dd1;
    const scalar contrUTilde  = UTilde & normal;
    const scalar aTilde = Foam::sqrt(max(0.0, (gamma-1)*(HTilde-0.5*sqr(contrUTilde))));

    // Step 3: compute signal speeds for face:
    const scalar SLeft  = min(qLeft - a_L, contrUTilde - aTilde);
    const scalar SRight = max(contrUTilde + aTilde, qRight + a_R);
    const scalar SStar = (rho_R*qRight*(SRight-qRight) -
                          rho_L*qLeft *(SLeft - qLeft) + p_L - p_R )/
                         stabilise((rho_R*(SRight-qRight)-rho_L*(SLeft-qLeft)),VSMALL);
    
    // Compute pressure in star region from the right side
    const scalar pStarRight = rho_R*(qRight-SRight)*(qRight-SStar) + p_R;
    // Should be equal to the left side
    const scalar pStarLeft  = rho_L*(qLeft - SLeft)*(qLeft -SStar) + p_L;
    // Give a warning if this is not the case
    if (mag(pStarRight - pStarLeft) > 1e-6)
    {
        Info << "mag(pStarRight-pStarLeft) > VSMALL " << endl;
    }

    // Use pStarRight for pStar, as in theory, pStarRight == pStarLeft
    const scalar pStar = pStarRight;

    // Step 4: upwinding - compute states:
    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;
    scalar pState = 0.0;

    if (pos(SLeft))
    {
        // compute F_l
        convectionSpeed = qLeft;
        rhoState  = rho_L;
        rhoUState = rhoU_L;
        rhoEState = rhoE_L;
        pState = p_L;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar), VSMALL);

        // Compute left star region
        convectionSpeed = SStar;
        rhoState  = omegaLeft*(SLeft - qLeft)*rho_L;
        rhoUState = omegaLeft*((SLeft - qLeft)*rhoU_L + (pStar - p_L)*normal);
        rhoEState = omegaLeft*((SLeft - qLeft)*rhoE_L - p_L*qLeft + pStar*SStar);
        pState = pStar;
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), VSMALL);

        // compute right star region
        convectionSpeed = SStar;
        rhoState  = omegaRight*(SRight - qRight)*rho_R;
        rhoUState = omegaRight*((SRight - qRight)*rhoU_R + (pStar - p_R)*normal);
        rhoEState = omegaRight*((SRight - qRight)*rhoE_R - p_R*qRight + pStar*SStar);
        pState = pStar;
    }
    else if (neg(SRight))
    {
        // compute F_r
        convectionSpeed = qRight;
        rhoState  = rho_R;
        rhoUState = rhoU_R;
        rhoEState = rhoE_R;
        pState = p_R;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }

    rhoFlux  = convectionSpeed*rhoState;
    rhoUFlux = convectionSpeed*rhoUState+pState*normal;
    rhoEFlux = convectionSpeed*(rhoEState+pState);
}

// ************************************************************************* //
