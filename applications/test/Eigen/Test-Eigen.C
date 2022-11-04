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
    Test-Eigen

Description
    This is the OpenFOAM Eigen test from buaaymh.

\*---------------------------------------------------------------------------*/

#include <iostream>
#include "eulerConsVar3rdSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Foam::Mat5X5 L, R;
    const scalar gamma = 1.4;
    const scalar rho  = 1.0;
    const vector U = vector(1.0, 1.0, 1.0);
    const scalar p = 2.0;
    const vector rhoU = rho*U;
    const scalar rhoE = p/(gamma-1.0) + 0.5*rho*magSqr(U);
    const vector normal = vector(0.0, 0.0, 1.0);
    evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal, gamma);
    std::cout << L << nl << R << std::endl;
    Foam::Mat5X5 I = Foam::Mat5X5::Identity();
    Info << (L * R - I).cwiseAbs().maxCoeff() << endl;
    Info << (R * L - I).cwiseAbs().maxCoeff() << endl;

    return 0;
}

// ************************************************************************* //