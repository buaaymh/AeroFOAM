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
    Blade.

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
    // Constructors
    Blade() {};
    // Destructor - default
    virtual ~Blade() {};
    
    // Member Functions
    virtual scalar minRadius() const = 0;

    virtual scalar maxRadius() const = 0;

    virtual scalar chord(scalar r) const = 0;
    
    virtual scalar twist(scalar r) const = 0;

    virtual std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const = 0;
};

class CaradonnaTung: public Blade
{
public:
    // Constructors
    CaradonnaTung();
    // Destructor - default
    ~CaradonnaTung() noexcept = default;

    scalar minRadius() const { return 0.21717; }

    scalar maxRadius() const { return 1.143; }

    scalar chord(scalar r) const { return 0.191; }
    
    scalar twist(scalar r) const { return 0; }

    std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const;

    static const std::array<scalar, 16> lift_;
    static const std::array<scalar, 16> drag_;
};

class UH60A: public Blade
{
public:
    
    UH60A();
    
    // Destructor - default
    ~UH60A() noexcept = default;

    scalar minRadius() const { return 0.1; }

    scalar maxRadius() const { return 1.0; }

    scalar chord(scalar r) const { return 0.1; }
    
    scalar twist(scalar r) const { return 0; }

    std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const;

    static const std::array<scalar, 181> lift_;
    static const std::array<scalar, 181> drag_;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //