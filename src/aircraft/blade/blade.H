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

    virtual scalar aspectRatio() const = 0;
    
    virtual scalar twist(scalar r) const = 0;

    virtual std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const = 0;
};

class Airfoil2D: public Blade
{
public:
    // Constructors
    Airfoil2D();
    // Destructor - default
    ~Airfoil2D() noexcept = default;

    scalar minRadius() const { return -GREAT; }

    scalar maxRadius() const { return  GREAT; }
    
    scalar aspectRatio() const { return GREAT; }

    scalar chord(scalar r) const { return 1; }
    
    scalar twist(scalar r) const { return 0; }

    std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const { return {1, 0}; }
};

class RectangularWing: public Blade
{
public:
    // Constructors
    RectangularWing();
    // Destructor - default
    ~RectangularWing() noexcept = default;

    scalar minRadius() const { return -0.5; }

    scalar maxRadius() const { return  0.5; }
    
    scalar aspectRatio() const { return  10; }

    scalar chord(scalar r) const { return 0.1; }
    
    scalar twist(scalar r) const { return 0; }

    std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const { return {1, 0}; }
};

class EllipticWing: public Blade
{
public:
    // Constructors
    EllipticWing();
    // Destructor - default
    ~EllipticWing() noexcept = default;

    scalar minRadius() const { return -0.5; }

    scalar maxRadius() const { return  0.5; }
    
    scalar aspectRatio() const { return  32.0/constant::mathematical::pi; }

    scalar chord(scalar r) const { return 0.125*sqrt(1.0-sqr(2.0*r)); }
    
    scalar twist(scalar r) const { return 0; }

    std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const;
};

class CaradonnaTung: public Blade
{
public:
    // Constructors
    CaradonnaTung();
    // Destructor - default
    ~CaradonnaTung() noexcept = default;

    scalar minRadius() const { return 0.1143; }

    scalar maxRadius() const { return 1.143; }

    scalar chord(scalar r) const { return 0.191; }
    
    scalar aspectRatio() const { return  6; }
    
    scalar twist(scalar r) const { return 0; }

    std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const;

    static const std::array<scalar, 41> lift_;
    static const std::array<scalar, 41> drag_;
};

class RotorFuselage: public Blade
{
public:
    
    RotorFuselage();
    
    // Destructor - default
    ~RotorFuselage() noexcept = default;

    scalar minRadius() const { return 0.012344; }

    scalar maxRadius() const { return 0.4572; }

    scalar chord(scalar r) const { return 0.086; }
    
    scalar aspectRatio() const { return  5.31628; }
    
    scalar twist(scalar r) const { return 0; }

    std::pair<scalar, scalar> Cl_Cd(scalar Ma, scalar r, scalar deg) const;

    static const std::array<scalar, 181> lift_;
    static const std::array<scalar, 181> drag_;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //