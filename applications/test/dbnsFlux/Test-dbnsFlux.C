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

Application
    Test-roeFlux

Description
    This is the OpenFOAM Euler riemann solver test from buaaymh.

\*---------------------------------------------------------------------------*/

#include <memory>

#include "fvCFD.H"
#include "AUSMplusUpFlux.H"
#include "hllcFlux.H"
#include "roeFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    std::unique_ptr<dbnsFlux> roe = std::make_unique<Foam::roeFlux>();
    std::unique_ptr<dbnsFlux> hllc = std::make_unique<Foam::hllcFlux>();
    std::unique_ptr<dbnsFlux> AUSMplusUp = std::make_unique<Foam::AUSMplusUpFlux>();
    vector  normal = vector(1.0, 0.0, 0.0);
    scalar  gamma  = 1.4;

    //Flux
    scalar rhoFlux;
    vector rhoUFlux = vector::zero;
    scalar rhoEFlux;
    
    Info << "Sod problem:" << endl;
    // Left state
    scalar rho_L = 1.0; vector U_L = vector(0.0, 0.0, 0.0); scalar T_L = 1.4;
    // Right state
    scalar rho_R = 0.125; vector U_R = vector(0.0, 0.0, 0.0); scalar T_R = 1.12;
    // Mid state
    scalar rho_M = 0.426319; vector U_M = vector(0.927453, 0.0, 0.0); scalar p_M = 0.303130;
    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_M, U_M, p_M, normal, gamma);
    Info << "Exact flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    roe->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                      rho_L, rho_R, U_L, U_R, T_L, T_R,
                      normal, gamma);
    Info << "roe   flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    hllc->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                       rho_L, rho_R, U_L, U_R, T_L, T_R,
                       normal, gamma);
    Info << "hllc  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    AUSMplusUp->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_L, rho_R, U_L, U_R, T_L, T_R,
                             normal, gamma);
    Info << "AUSM  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    Info << endl;

    Info << "ShockCollision problem:" << endl;
    // Left state
    rho_L = 5.99924; U_L = vector(19.5975, 0.0, 0.0); T_L = 107.555557;
    // Right state
    rho_R = 5.99924; U_R = vector(6.19633, 0.0, 0.0); T_R = 10.7555557;
    // Mid state
    rho_M = 5.99924; U_M = vector(19.5975, 0.0, 0.0); p_M = 460.894;
    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_M, U_M, p_M, normal, gamma);
    Info << "Exact flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    roe->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                      rho_L, rho_R, U_L, U_R, T_L, T_R,
                      normal, gamma);
    Info << "roe   flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    hllc->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                       rho_L, rho_R, U_L, U_R, T_L, T_R,
                       normal, gamma);
    Info << "hllc  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    AUSMplusUp->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_L, rho_R, U_L, U_R, T_L, T_R,
                             normal, gamma);
    Info << "AUSM  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    Info << endl;

    Info << "BlastFromLeft problem:" << endl;
    // Left state
    rho_L = 1.0; U_L = vector(0.0, 0.0, 0.0); T_L = 1400.0;
    // Right state
    rho_R = 1.0; U_R = vector(0.0, 0.0, 0.0); T_R = 0.014;
    // Mid state
    rho_M = 0.575062; U_M = vector(19.5975, 0.0, 0.0); p_M = 460.8938;
    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_M, U_M, p_M, normal, gamma);
    Info << "Exact flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    roe->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                      rho_L, rho_R, U_L, U_R, T_L, T_R,
                      normal, gamma);
    Info << "roe   flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    hllc->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                       rho_L, rho_R, U_L, U_R, T_L, T_R,
                       normal, gamma);
    Info << "hllc  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    AUSMplusUp->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_L, rho_R, U_L, U_R, T_L, T_R,
                             normal, gamma);
    Info << "AUSM  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    Info << endl;

    Info << "BlastFromLeft problem:" << endl;
    // Left state
    rho_L = 1.0; U_L = vector(0.0, 0.0, 0.0); T_L = 0.014;
    // Right state
    rho_R = 1.0; U_R = vector(0.0, 0.0, 0.0); T_R = 140.0;
    // Mid state
    rho_M = 0.575062; U_M = vector(-6.196328, 0.0, 0.0); p_M = 46.09504;
    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_M, U_M, p_M, normal, gamma);
    Info << "Exact flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    roe->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                      rho_L, rho_R, U_L, U_R, T_L, T_R,
                      normal, gamma);
    Info << "roe   flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    hllc->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                       rho_L, rho_R, U_L, U_R, T_L, T_R,
                       normal, gamma);
    Info << "hllc  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    AUSMplusUp->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_L, rho_R, U_L, U_R, T_L, T_R,
                             normal, gamma);
    Info << "AUSM  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    Info << endl;

    Info << "AlmostVaccumed problem:" << endl;
    // Left state
    rho_L = 1.0; U_L = vector(-2.0, 0.0, 0.0); T_L = 0.56;
    // Right state
    rho_R = 1.0; U_R = vector(+2.0, 0.0, 0.0); T_R = 0.56;
    // Mid state
    rho_M = 0.21852; U_M = vector(0.0, 0.0, 0.0); p_M = 0.001894;
    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_M, U_M, p_M, normal, gamma);
    Info << "Exact flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    roe->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                      rho_L, rho_R, U_L, U_R, T_L, T_R,
                      normal, gamma);
    Info << "roe   flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    hllc->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                       rho_L, rho_R, U_L, U_R, T_L, T_R,
                       normal, gamma);
    Info << "hllc  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    AUSMplusUp->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_L, rho_R, U_L, U_R, T_L, T_R,
                             normal, gamma);
    Info << "AUSM  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    Info << endl;

    Info << "Vaccumed problem:" << endl;
    // Left state
    rho_L = 1.0; U_L = vector(-4.0, 0.0, 0.0); T_L = 0.56;
    // Right state
    rho_R = 1.0; U_R = vector(+4.0, 0.0, 0.0); T_R = 0.56;
    // Mid state
    rho_M = 0.0; U_M = vector(0.0, 0.0, 0.0); p_M = 0.0;
    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_M, U_M, p_M, normal, gamma);
    Info << "Exact flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    roe->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                      rho_L, rho_R, U_L, U_R, T_L, T_R,
                      normal, gamma);
    Info << "roe   flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    hllc->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                       rho_L, rho_R, U_L, U_R, T_L, T_R,
                       normal, gamma);
    Info << "hllc  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    AUSMplusUp->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_L, rho_R, U_L, U_R, T_L, T_R,
                             normal, gamma);
    Info << "AUSM  flux = "  << rhoFlux << rhoUFlux << rhoEFlux << endl;
    Info << endl;

    return 0;
}

// ************************************************************************* //