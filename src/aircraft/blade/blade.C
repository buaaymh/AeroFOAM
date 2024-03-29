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

#include "blade.H"

Foam::Airfoil2D::Airfoil2D() : Blade() {}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::RectangularWing::RectangularWing() : Blade() {}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::EllipticWing::EllipticWing() : Blade() {}

std::pair<scalar, scalar> Foam::EllipticWing::Cl_Cd(scalar Ma, scalar r, scalar deg) const
{
    return { 2*sqr(constant::mathematical::pi)*deg/180, 0 };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::HorizontalStabilator::HorizontalStabilator() : Blade() {}

const std::array<scalar, 29> HorizontalStabilator::lift_
    = NACA0014::getLiftCoefficients();

const std::array<scalar, 29> HorizontalStabilator::drag_
    = NACA0014::getDragCoefficients();

std::pair<scalar, scalar> Foam::HorizontalStabilator::Cl_Cd(scalar Ma, scalar r, scalar deg) const
{
    scalar Cl, Cd;
    scalar index = mag(deg);
    if (index >= 28) { Cl = lift_[28]; Cd = drag_[28]; }
    else
    {
        label deg_floor = floor(index);
        label deg_ceil = ceil(index);
        if (deg_floor != deg_ceil)
        {
            Cl = (deg_ceil - index) *lift_[deg_floor]
               + (index - deg_floor)*lift_[deg_ceil];
            Cd = (deg_ceil - index) *drag_[deg_floor]
               + (index - deg_floor)*drag_[deg_ceil];
        }
        else { Cl = lift_[deg_floor]; Cd = drag_[deg_floor]; }
    }
    if (deg < 0) Cl = -Cl;
    return { Cl, Cd };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::CaradonnaTung::CaradonnaTung()
:
    Blade()
{}

const std::array<scalar, 23> CaradonnaTung::lift_
    = NACA0012::getLiftCoefficients();

const std::array<scalar, 23> CaradonnaTung::drag_
    = NACA0012::getDragCoefficients();

std::pair<scalar, scalar> Foam::CaradonnaTung::Cl_Cd(scalar Ma, scalar r, scalar deg) const
{
    scalar Cl, Cd;
    scalar index = mag(deg);
    if (index >= 22) { Cl = lift_[22]; Cd = drag_[22]; }
    else
    {
        label deg_floor = floor(index);
        label deg_ceil = ceil(index);
        if (deg_floor != deg_ceil)
        {
            Cl = (deg_ceil - index) *lift_[deg_floor]
               + (index - deg_floor)*lift_[deg_ceil];
            Cd = (deg_ceil - index) *drag_[deg_floor]
               + (index - deg_floor)*drag_[deg_ceil];
        }
        else { Cl = lift_[deg_floor]; Cd = drag_[deg_floor]; }
    }
    if (deg < 0) Cl = -Cl;
    return { Cl, Cd };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::RotorFuselage::RotorFuselage()
:
    Blade()
{}

const std::array<scalar, 29> RotorFuselage::lift_
    = NACA0015::getLiftCoefficients();

const std::array<scalar, 29> RotorFuselage::drag_
    = NACA0015::getDragCoefficients();

std::pair<scalar, scalar> Foam::RotorFuselage::Cl_Cd(scalar Ma, scalar r, scalar deg) const
{
    scalar Cl, Cd;
    scalar index = mag(deg);
    if (index >= 28) { Cl = lift_[28]; Cd = drag_[28]; }
    else
    {
        label deg_floor = floor(index);
        label deg_ceil = ceil(index);
        if (deg_floor != deg_ceil)
        {
            Cl = (deg_ceil - index) *lift_[deg_floor]
               + (index - deg_floor)*lift_[deg_ceil];
            Cd = (deg_ceil - index) *drag_[deg_floor]
               + (index - deg_floor)*drag_[deg_ceil];
        }
        else { Cl = lift_[deg_floor]; Cd = drag_[deg_floor]; }
    }
    if (deg < 0) Cl = -Cl;
    return { Cl, Cd };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::UH60A_Main::UH60A_Main()
:
    Blade()
{}

const std::array<scalar, 181> UH60A_Main::lift_
    = SC1095::getLiftCoefficients();

const std::array<scalar, 181> UH60A_Main::drag_
    = SC1095::getDragCoefficients();

std::pair<scalar, scalar> Foam::UH60A_Main::Cl_Cd(scalar Ma, scalar r, scalar deg) const
{
    scalar Cl, Cd;
    scalar index = 0.5*(deg + 180);
    label deg_floor = floor(index);
    label deg_ceil = ceil(index);
    Cl = (deg_ceil - index) *lift_[deg_floor]
       + (index - deg_floor)*lift_[deg_ceil];
    Cd = (deg_ceil - index) *drag_[deg_floor]
       + (index - deg_floor)*drag_[deg_ceil];
    return { Cl, Cd };
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::UH60A_Tail::UH60A_Tail()
:
    Blade()
{}

const std::array<scalar, 181> UH60A_Tail::lift_
    = SC1095::getLiftCoefficients();

const std::array<scalar, 181> UH60A_Tail::drag_
    = SC1095::getDragCoefficients();

std::pair<scalar, scalar> Foam::UH60A_Tail::Cl_Cd(scalar Ma, scalar r, scalar deg) const
{
    scalar Cl, Cd;
    scalar index = 0.5*(deg + 180);
    label deg_floor = floor(index);
    label deg_ceil = ceil(index);
    Cl = (deg_ceil - index) *lift_[deg_floor]
       + (index - deg_floor)*lift_[deg_ceil];
    Cd = (deg_ceil - index) *drag_[deg_floor]
       + (index - deg_floor)*drag_[deg_ceil];
    return { Cl, Cd };
}