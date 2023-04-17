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
#include "AUSMplusUpFlux.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::AUSMplusUpFlux::evaluateFlux
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
    const vector& normal
) const
{   
    // left & right state
    const vector rhoU_L = rho_L*U_L;
    const scalar rhoE_L = p_L/(Gamma-1) + 0.5*rho_L*magSqr(U_L);
    const scalar H_L = (rhoE_L + p_L)/rho_L;
    
    const vector rhoU_R = rho_R*U_R;
    const scalar rhoE_R = p_R/(Gamma-1) + 0.5*rho_R*magSqr(U_R);
    const scalar H_R = (rhoE_R + p_R)/rho_R;

    const scalar qLeft  = (U_L & normal);
    const scalar qRight = (U_R & normal);

    scalar temp = 2.0*(Gamma-1)/(Gamma+1)*0.5*(H_L+H_R);
    const scalar aStar = sqrt(temp);
    const scalar aHatLeft  = temp / max(aStar, qLeft);
    const scalar aHatRight = temp / max(aStar,-qRight);
    const scalar aTilde = min(aHatLeft, aHatRight);
    const scalar rhoTilde = 0.5*(rho_L+rho_R);
    
    temp = 1.0/aTilde;
    const scalar sqrMaDash = 0.5*(sqr(qLeft)+sqr(qRight))*sqr(temp);
    
    const scalar MaRelLeft  = qLeft *temp;
    const scalar MaRelRight = qRight*temp;
    
    const scalar magMaRelLeft  = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);
    
    const scalar Ma1PlusLeft   = 0.5*(MaRelLeft +magMaRelLeft );
    const scalar Ma1MinusRight = 0.5*(MaRelRight-magMaRelRight);

    const scalar Ma2PlusLeft   =  0.25*sqr(MaRelLeft +1.0);
    const scalar Ma2PlusRight  =  0.25*sqr(MaRelRight+1.0);
    const scalar Ma2MinusLeft  = -0.25*sqr(MaRelLeft -1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight-1.0);
    
    const scalar Ma4BetaPlusLeft   = ((magMaRelLeft  >= 1.0) ? Ma1PlusLeft   : (Ma2PlusLeft  *(1.0-16.0*beta_*Ma2MinusLeft)));
    const scalar Ma4BetaMinusRight = ((magMaRelRight >= 1.0) ? Ma1MinusRight : (Ma2MinusRight*(1.0+16.0*beta_*Ma2PlusRight)));
        
    const scalar Mp = -Kp_/fa_*max(1.0-sigma_*sqrMaDash,0.0)*(p_R-p_L)/(rhoTilde*sqr(aTilde));

    const scalar P5alphaPlusLeft   = ((magMaRelLeft  >= 1.0) ?
    (Ma1PlusLeft/MaRelLeft)    : (Ma2PlusLeft  *(( 2.0-MaRelLeft) -16.0*alpha_*MaRelLeft *Ma2MinusLeft )));
    const scalar P5alphaMinusRight = ((magMaRelRight >= 1.0) ?
    (Ma1MinusRight/MaRelRight) : (Ma2MinusRight*((-2.0-MaRelRight)+16.0*alpha_*MaRelRight*Ma2PlusRight)));
    
    const scalar pU = -Ku_*P5alphaPlusLeft*P5alphaMinusRight*(rho_L+rho_R)*(fa_*aTilde)*(qRight-qLeft);
    
    const scalar MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight + Mp;
    const scalar pTilde = p_L*P5alphaPlusLeft + p_R*P5alphaMinusRight + pU;
    
    const scalar URelTilde = MaRelTilde*aTilde;
    const scalar magURelTilde = mag(MaRelTilde)*aTilde;
    // There is a typo in Luo et. al, J. Comp. Physics 194 (2004), Chap 4.2 Eq. 4.8
    // refer to the origial Paper from Liou, J. Comp. Physics 129 (1996), Chap4, Eq. 42
    rhoFlux  = 0.5*(URelTilde*(rho_L +rho_R) -magURelTilde*(rho_R -rho_L));
    rhoUFlux = 0.5*(URelTilde*(rhoU_L+rhoU_R)-magURelTilde*(rhoU_R-rhoU_L))+pTilde*normal;
    rhoEFlux = 0.5*(URelTilde*((rhoE_L+p_L)+(rhoE_R+p_R))-magURelTilde*((rhoE_R+p_R)-(rhoE_L+p_L)));
}

// ************************************************************************* //
