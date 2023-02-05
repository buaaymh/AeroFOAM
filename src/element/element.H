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

#ifndef element_H
#define element_H

#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Gauss Namespace  Declaration
\*---------------------------------------------------------------------------*/

namespace Tria3
{
    static constexpr std::array<scalar, 3> x{0.0, 0.5, 0.5};
    static constexpr std::array<scalar, 3> y{0.5, 0.0, 0.5};
    static constexpr std::array<scalar, 3> w{0.3333333333333333, 0.3333333333333333, 0.3333333333333333};
}

namespace Tria4
{
    static constexpr std::array<scalar, 4> x{0.3333333333333333, 0.6, 0.2, 0.2};
    static constexpr std::array<scalar, 4> y{0.3333333333333333, 0.2, 0.2, 0.6};
    static constexpr std::array<scalar, 4> w{-0.5624999999999998, 0.5208333333333332,
                                              0.5208333333333332, 0.5208333333333332};
}

namespace Quad4
{
    static constexpr std::array<scalar, 4> x{-0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, +0.577350269189626};
    static constexpr std::array<scalar, 4> y{+0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, -0.577350269189626};
    static constexpr std::array<scalar, 4> w{0.25, 0.25, 0.25, 0.25};
}

namespace Quad9
{
    static constexpr std::array<scalar, 9> x{-0.774596669241483,                0.0, +0.774596669241483,
                                             -0.774596669241483,                0.0, +0.774596669241483,
                                             -0.774596669241483,                0.0, +0.774596669241483};
    static constexpr std::array<scalar, 9> y{+0.774596669241483, +0.774596669241483, +0.774596669241483,
                                                            0.0,                0.0,                0.0,
                                             -0.774596669241483, -0.774596669241483, -0.774596669241483};
    static constexpr std::array<scalar, 9> w{0.07716049382716049383, 0.12345679012345679012, 0.07716049382716049383,
                                             0.12345679012345679012, 0.19753086419753086420, 0.12345679012345679012,
                                             0.07716049382716049383, 0.12345679012345679012, 0.07716049382716049383};
}

namespace Hexa8
{
    static constexpr std::array<scalar, 8> x{-0.577350269189626, +0.577350269189626, +0.577350269189626, -0.577350269189626,
                                             -0.577350269189626, +0.577350269189626, +0.577350269189626, -0.577350269189626};
    static constexpr std::array<scalar, 8> y{-0.577350269189626, -0.577350269189626, +0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, -0.577350269189626, +0.577350269189626, +0.577350269189626};
    static constexpr std::array<scalar, 8> z{-0.577350269189626, -0.577350269189626, -0.577350269189626, -0.577350269189626,
                                             +0.577350269189626, +0.577350269189626, +0.577350269189626, +0.577350269189626};
    static constexpr std::array<scalar, 8> w{0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};
}

namespace Prism6
{
    static constexpr std::array<scalar, 6> x{+1.0/6.0, +2.0/3.0, +1.0/6.0,
                                             +1.0/6.0, +2.0/3.0, +1.0/6.0};
    static constexpr std::array<scalar, 6> y{+1.0/6.0, +1.0/6.0, +2.0/3.0,
                                             +1.0/6.0, +1.0/6.0, +2.0/3.0};
    static constexpr std::array<scalar, 6> z{-0.577350269189626, -0.577350269189626, -0.577350269189626,
                                             +0.577350269189626, +0.577350269189626, +0.577350269189626};
    static constexpr std::array<scalar, 6> w{1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};
}

namespace Tetra4
{
    static constexpr std::array<scalar, 4> x{0.1381966011250105, 0.1381966011250105,
                                             0.1381966011250105, 0.5854101966249685};
    static constexpr std::array<scalar, 4> y{0.1381966011250105, 0.1381966011250105,
                                             0.5854101966249685, 0.1381966011250105};
    static constexpr std::array<scalar, 4> z{0.1381966011250105, 0.5854101966249685,
                                             0.1381966011250105, 0.1381966011250105};
    static constexpr std::array<scalar, 4> w{0.25, 0.25, 0.25, 0.25};
}

/*---------------------------------------------------------------------------*\
                           Class Face Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class Face
{
public:

    // Destructor - default
    virtual ~Face() {};

    // Member Functions
    virtual label size() const = 0;

    virtual scalar weight(label i) const = 0;
    
    virtual const vector& at(label i) const = 0;

};

class Triangle3: public Face
{
public:
    // Constructors - default
    Triangle3(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Triangle3() {};

    // Member Functions
    label size() const { return 3; };

    scalar weight(label i) const { return Tria3::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};

class Triangle4: public Face
{
public:
    // Constructors - default
    Triangle4(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Triangle4() {};

    // Member Functions
    label size() const { return 4; };

    scalar weight(label i) const { return Tria4::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};

class Quadrangle4: public Face
{
public:
    // Constructors - default
    Quadrangle4(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Quadrangle4() {};

    // Member Functions
    label size() const { return 4; };

    scalar weight(label i) const { return Quad4::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};

class Quadrangle9: public Face
{
public:
    // Constructors - default
    Quadrangle9(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Quadrangle9() {};

    // Member Functions
    label size() const { return 9; };

    scalar weight(label i) const { return Quad9::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};

void build2ndFace(const fvMesh& mesh, const label& faceI, std::unique_ptr<Face>& face);

void build4stFace(const fvMesh& mesh, const label& faceI, std::unique_ptr<Face>& face);

void gaussHexa8
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<vector>& quadPoints
);

void gaussPrism6
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<vector>& quadPoints
);

void gaussTetra4
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<vector>& quadPoints
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //