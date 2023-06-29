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

Description
    Definitions of values and funcitons.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef definitions_H
#define definitions_H

#include <array>
#include <map>
#include <memory>
#include <vector>
#include <Eigen/Dense>

#include "fvCFD.H"
#include "IOstreams.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Constants Definitions
\*---------------------------------------------------------------------------*/
static constexpr scalar Gamma   = 1.4;      /* Gas constant */
static constexpr scalar Pr_Lam  = 0.72;     /* Prandtl numbers */
static constexpr scalar Pr_Turb = 0.9;      /* Prandtl numbers */


using Arr4X1 = Eigen::Array<scalar, 4, 1>;
using Arr8X1 = Eigen::Array<scalar, 8, 1>;

using Col3X1 = Eigen::Matrix<scalar, 3, 1>;
using Col4X1 = Eigen::Matrix<scalar, 4, 1>;
using Col5X1 = Eigen::Matrix<scalar, 5, 1>;
using Col6X1 = Eigen::Matrix<scalar, 6, 1>;
using Col9X1 = Eigen::Matrix<scalar, 9, 1>;

using Mat3X2 = Eigen::Matrix<scalar, 3, 2>;
using Mat3X3 = Eigen::Matrix<scalar, 3, 3>;
using Mat3X5 = Eigen::Matrix<scalar, 3, 5>;

using Mat4X2 = Eigen::Matrix<scalar, 4, 2>;
using Mat4X3 = Eigen::Matrix<scalar, 4, 3>;

using Mat5X3 = Eigen::Matrix<scalar, 5, 3>;
using Mat5X5 = Eigen::Matrix<scalar, 5, 5>;
using Mat5X6 = Eigen::Matrix<scalar, 5, 6>;
using Mat5X9 = Eigen::Matrix<scalar, 5, 9>;

using Mat6X3 = Eigen::Matrix<scalar, 6, 3>;
using Mat6X5 = Eigen::Matrix<scalar, 6, 5>;
using Mat6X6 = Eigen::Matrix<scalar, 6, 6>;

using Mat8X3 = Eigen::Matrix<scalar, 8, 3>;

using Mat9X3 = Eigen::Matrix<scalar, 9, 3>;
using Mat9X5 = Eigen::Matrix<scalar, 9, 5>;
using Mat9X9 = Eigen::Matrix<scalar, 9, 9>;

/*---------------------------------------------------------------------------*\
                   Structure fluidProperties Declaration
\*---------------------------------------------------------------------------*/

struct fluidProperties
{
    bool withSourceTerm;
    word simulationType;

    scalar AOA,         /* Angle of attack */
           Mach_inf,    /* Free-stream Mach number */
           Re_inf,      /* Reynolds number */
           T_inf,       /* Free-stream Temperature */
           p_inf;       /* Free-stream pressure */
};

/*---------------------------------------------------------------------------*\
                           Methods Declaration
\*---------------------------------------------------------------------------*/

scalar MinMod
(
    const scalar& a,
    const scalar& b
);

void givensRotation
(
    const scalar& h,
    const scalar& beta,
    scalar& c,
    scalar& s
);

scalar Venkat
(
    const scalar& project,
    const scalar& deltaMin,
    const scalar& deltaMax,
    const scalar& eps2
);

scalar Barth
(
    const scalar& project,
    const scalar& deltaMin,
    const scalar& deltaMax
);

vector vorticity(const tensor& UGrad);

scalar func
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& dRho,
    const vector& dRhoU,
    const scalar& dRhoE,
    const scalar& x
);

scalar diff
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& dRho,
    const vector& dRhoU,
    const scalar& dRhoE,
    const scalar& x
);

scalar solveForPressure
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& dRho,
    const vector& dRhoU,
    const scalar& dRhoE,
    const scalar& pMin,
    const scalar& xMax
);

void evaluateEigenMatrix
(
    Mat5X5& L,
    Mat5X5& R,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& normal
);

scalar deg2rad(scalar deg);

scalar rad2deg(scalar rad);

std::pair<scalar, scalar> cosSin(scalar deg);

void rotateX(vector& y_unit, vector& z_unit, scalar deg);

void rotateY(vector& x_unit, vector& z_unit, scalar deg);

void rotateZ(vector& x_unit, vector& y_unit, scalar deg);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //