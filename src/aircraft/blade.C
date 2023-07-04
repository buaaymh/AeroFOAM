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

Foam::Blade::Blade
(
    const fvMesh& mesh
)
:
    mesh_(mesh)
{
    nSpans_ = mesh_.solution().subDict("blade").lookup<label>("nSpans");
    minRadius_ = mesh_.solution().subDict("blade").lookup<scalar>("minRadius");
    maxRadius_ = mesh_.solution().subDict("blade").lookup<scalar>("maxRadius");
    chord_ = mesh_.solution().subDict("blade").lookup<scalar>("chord");
    aspectRatio_ = mesh_.solution().subDict("blade").lookup<scalar>("aspectRatio");
    twist_ = mesh_.solution().subDict("blade").lookup<scalar>("twist");
    dSpan_ = (maxRadius_ - minRadius_) / nSpans_;
    eps_cStar_ = 0.02*aspectRatio_*constant::mathematical::pi;
    c0_ = 4*chord_/constant::mathematical::pi;
    eps0_ = eps_cStar_*c0_;
}

scalar Foam::Blade::GaussianRadius
(
    scalar r
) const
{
    r -= 0.5 * maxRadius_;
    const scalar dx = dSpan_/1.5;
    const scalar cStar = c0_*sqrt(1.0-sqr(2*r/maxRadius_));
    const scalar eps = eps_cStar_ * cStar;
    return max(eps, dx);
}

template<class Airfoil>
Foam::Rectangle<Airfoil>::Rectangle
(
    const fvMesh& mesh
)
:
    Blade(mesh)
{}