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

Class
    Face

Description
    Face for 3rd order solver.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef blade_H
#define blade_H

#include "airfoil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Blade Class Declaration
\*---------------------------------------------------------------------------*/

class Blade
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    Blade(const Blade&);
    //- Disallow default bitwise assignment
    void operator=(const Blade&);

public:
    // Protected Member
    const fvMesh& mesh_;
    label  nSpans_;
    scalar minRadius_;
    scalar maxRadius_;
    scalar chord_;
    scalar aspectRatio_;
    scalar twist_;
    scalar dSpan_;
    scalar eps_cStar_;
    scalar c0_;
    scalar eps0_;

    // Constructors
    Blade(const fvMesh& mesh);

    // Destructor - default
    virtual ~Blade() {};
    
    // Member Functions
    label  numOfSpans() const { return nSpans_; }

    scalar dSpan() const { return dSpan_; }

    scalar minRadius() const { return minRadius_; }

    scalar maxRadius() const { return maxRadius_; }

    scalar GaussianRadius(scalar r) const;

    scalar twist() const { return twist_; }

    scalar chord() const { return chord_; }

    scalar eps0() const { return eps0_; }

    virtual scalar Cl(scalar r, scalar deg) const = 0;
    virtual scalar Cd(scalar r, scalar deg) const = 0;
};

template<class Airfoil>
class Rectangle: public Blade
{
public:
    
    Rectangle(const fvMesh& mesh);
    
    // Destructor - default
    ~Rectangle() noexcept = default;

    scalar Cl(scalar r, scalar deg) const
    {
        return lift_[(deg + 180) / 2];
    }

    scalar Cd(scalar r, scalar deg) const
    {
        return drag_[(deg + 180) / 2];
    }

    static const std::array<scalar, 181> lift_;
    static const std::array<scalar, 181> drag_;
};

template<class Airfoil>
const std::array<scalar, 181> Rectangle<Airfoil>::lift_
    = Airfoil::getLiftCoefficients();

template<class Airfoil>
const std::array<scalar, 181> Rectangle<Airfoil>::drag_
    = Airfoil::getDragCoefficients();

template class Rectangle<SC1095>;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //