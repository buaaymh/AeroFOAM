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

#include "fvCFD.H"
#include <iostream>
#include <time.h>
#include <vector>
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Eigen::Matrix<double, 2, 2> mat;
    mat << 2, 0,
           0, 2;
    Eigen::Matrix<float, 2, 2> mat_1 = mat.inverse().cast<float>();
    std::cout << mat_1 << std::endl;

    // clock_t start, end;
    // label n = 2000000;
    // start = std::clock();
    // std::vector<Eigen::Matrix<scalar, 6, 6>> eigenMatsD(n, Eigen::Matrix<scalar, 6, 6>::Ones());
    // std::vector<Eigen::Matrix<float, 6, 6>> eigenMatsF(n, Eigen::Matrix<float, 6, 6>::Ones());
    // std::vector<symmTensor> coefs(n, symmTensor::one);

    // Eigen::Matrix<scalar, 6, 1> temp = Eigen::Matrix<scalar, 6, 1>::Zero();
    // for (int j = 0; j != 500; j++)
    // {
    //     for (int i = 0; i != n; ++i)
    //     {
    //         Eigen::Matrix<scalar, 6, 1> col(coefs[i][0], coefs[i][1], coefs[i][2],
    //                                         coefs[i][3], coefs[i][4], coefs[i][5]);
    //         temp += eigenMatsD[i] * col;
    //         temp += eigenMatsF[i].cast<scalar>() * col;
    //     }
    //     temp *= 0;
    // }
    // end = std::clock();
    // Info <<"CPU Time = "<<scalar(end-start)/CLOCKS_PER_SEC<< "s" <<endl;

    return 0;
}

// ************************************************************************* //