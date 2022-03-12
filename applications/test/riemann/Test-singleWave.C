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
    Test-singleWave

Description
    This is the OpenFOAM singleWave test from buaaymh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singleWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    vector a(2.0, 0.0, 0.0);
    singleWave riemann1(a);
    Info << "The a of single wave is " << a << endl;
    vector normal(0.6, 0.8, 0.0);
    Info << "The normal of face is " << normal << endl;
    scalar uLeft{-1.0}, uRight{2.0}, flux{0.0};
    Info << "U on left  is " << uLeft << endl;
    Info << "U on right is " << uRight << endl;
    riemann1.evaluateFlux(flux, uLeft, uRight, normal);
    Info << "Flux is " << flux << endl;
    if (mag(flux + 1.2) < 1e-6)
    {
        Info << "Pass!" << endl << endl;
    }
    else
    {
        Info << "Fault!" << endl << endl;
    }

    a = vector(-2.0, 0.0, 0.0);
    singleWave riemann2(a);
    Info << "The a of single wave is " << a << endl;
    normal = vector(0.6, 0.8, 0.0);
    Info << "The normal of face is " << normal << endl;
    uLeft = -4.0, uRight = 2.0;
    Info << "U on left  is " << uLeft << endl;
    Info << "U on right is " << uRight << endl;
    riemann2.evaluateFlux(flux, uLeft, uRight, normal);
    Info << "Flux is " << flux << endl;
    if (mag(flux + 2.4) < 1e-6)
    {
        Info << "Pass!" << endl << endl;
    }
    else
    {
        Info << "Fault!" << endl << endl;
    }

    return 0;
}

// ************************************************************************* //