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
#include "euler3rdSolver.H"

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
    Foam::Mat5X5 I = Foam::Mat5X5::Identity();
    
    vector normal = vector(0.8, 0.6, 0.0);
    evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal);
    Info << "normal = " << normal << endl;
    std::cout << "L = \n" << L << nl << "R = \n" << R << std::endl;
    Info << (L * R - I).cwiseAbs().maxCoeff() << endl;
    Info << (R * L - I).cwiseAbs().maxCoeff() << endl;

    normal = vector(-0.6, 0.0, 0.8);
    evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal);
    Info << "normal = " << normal << endl;
    std::cout << "L = \n" << L << nl << "R = \n" << R << std::endl;
    Info << (L * R - I).cwiseAbs().maxCoeff() << endl;
    Info << (R * L - I).cwiseAbs().maxCoeff() << endl;

    normal = vector(0.0, 0.8, -0.6);
    evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal);
    Info << "normal = " << normal << endl;
    std::cout << "L = \n" << L << nl << "R = \n" << R << std::endl;
    Info << (L * R - I).cwiseAbs().maxCoeff() << endl;
    Info << (R * L - I).cwiseAbs().maxCoeff() << endl;

    Foam::Mat5X5 limiter
    {
        {0.2, 0.0, 0.0, 0.0, 0.0},
        {0.0, 0.4, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.6, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.7, 0.0},
        {0.0, 0.0, 0.0, 0.0, 0.8}
    };
    std::cout << "eigen limiter = \n" << limiter << std::endl;
    std::cout << "conse limiter = \n" << L*limiter*R << std::endl;

    Foam::Col5X1 coef(1.0, 2.0, 3.0, 4.0, 5.0);
    std::cout << "eigen coefficient = \n" << (R*limiter*L)*coef << std::endl;
    std::cout << "conse coefficient = \n" <<    coef            << std::endl;

    Eigen::Matrix<scalar, 3, 3> A;
    Eigen::Matrix<scalar, 3, 1> b;
    A << 2, 1, 3,
         1, 1, 1,
         1, 3, 4;
    b << 13, 6, 19;
    Eigen::Matrix<scalar, 3, 1> x = A.colPivHouseholderQr().solve(b);
    std::cout << "Matrix A = \n" << A << std::endl;
    std::cout << "Vector b = \n" << b << std::endl;
    std::cout << "Solution x = \n" << x << std::endl;
    return 0;
}

// ************************************************************************* //