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
                           Gauss Class Declaration
\*---------------------------------------------------------------------------*/

template<int nVertex, int nPoints>
class Gauss;

template<>
class Gauss<3, 3>
{
public:
    static constexpr auto buildLocalXCoords() { return std::array<scalar, 3>{0.0, 0.5, 0.5}; }
    static constexpr auto buildLocalYCoords() { return std::array<scalar, 3>{0.5, 0.0, 0.5}; }
    static constexpr auto buildLocalWeights() { return std::array<scalar, 3>{1.0/6.0, 1.0/6.0, 1.0/6.0}; }
};

template<>
class Gauss<3, 6>
{
public:
    static constexpr auto buildLocalXCoords()
    { 
        return std::array<scalar, 6>{0.8168475729804585, 0.0915762135097708, 0.0915762135097708,
                                     0.1081030181680702, 0.4459484909159649, 0.4459484909159649};
    }
    static constexpr auto buildLocalYCoords()
    { 
        return std::array<scalar, 6>{0.0915762135097708, 0.0915762135097708, 0.8168475729804585,
                                     0.4459484909159649, 0.4459484909159649, 0.1081030181680702};
    }
    static constexpr auto buildLocalWeights()
    { 
        return std::array<scalar, 6>{0.054975871827660935, 0.054975871827660935, 0.054975871827660935,
                                     0.111690794839005740, 0.111690794839005740, 0.111690794839005740};
    }
};

template<>
class Gauss<4, 4>
{
public:
    static constexpr auto buildLocalXCoords()
    {
        return std::array<scalar, 4>{-std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0),
                                     +std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0)};
    }
    static constexpr auto buildLocalYCoords()
    {
        return std::array<scalar, 4>{-std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0),
                                     +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0)};
    }
    static constexpr auto buildLocalWeights() { return std::array<scalar, 4>{1.0, 1.0, 1.0, 1.0}; }
};

template<>
class Gauss<4, 9>
{
public:
    static constexpr auto buildLocalXCoords()
    { 
        return std::array<scalar, 9>{-std::sqrt(0.6),             0.0, +std::sqrt(0.6),
                                     -std::sqrt(0.6),             0.0, +std::sqrt(0.6),
                                     -std::sqrt(0.6),             0.0, +std::sqrt(0.6)};
    }
    static constexpr auto buildLocalYCoords()
    { 
        return std::array<scalar, 9>{+std::sqrt(0.6), +std::sqrt(0.6), +std::sqrt(0.6),
                                                0.0 ,            0.0 ,            0.0 ,
                                     -std::sqrt(0.6), -std::sqrt(0.6), -std::sqrt(0.6)};
    }
    static constexpr auto buildLocalWeights()
    { 
        return std::array<scalar, 9>{25.0/81.0, 40.0/81.0, 25.0/81.0,
                                     40.0/81.0, 64.0/81.0, 40.0/81.0,
                                     25.0/81.0, 40.0/81.0, 25.0/81.0};
    }
};

namespace Quad4
{
    static const Arr4X1 x_hexa_i{-1, +1, +1, -1};
    static const Arr4X1 y_hexa_i{-1, -1, +1, +1};
}

namespace Hexa8
{
    static constexpr std::array<scalar, 8> x{-std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0),
                                             -std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0)};
    static constexpr std::array<scalar, 8> y{-std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0),
                                             -std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0)};
    static constexpr std::array<scalar, 8> z{-std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0),
                                             +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0)};
    static constexpr std::array<scalar, 8> w{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    static const Arr8X1 x_hexa_i{-1, +1, +1, -1, -1, +1, +1, -1};
    static const Arr8X1 y_hexa_i{-1, -1, +1, +1, -1, -1, +1, +1};
    static const Arr8X1 z_hexa_i{-1, -1, -1, -1, +1, +1, +1, +1};
}

namespace Prism6
{
    static constexpr std::array<scalar, 6> x{+1.0/6.0, +2.0/3.0, +1.0/6.0,
                                             +1.0/6.0, +2.0/3.0, +1.0/6.0};
    static constexpr std::array<scalar, 6> y{+1.0/6.0, +1.0/6.0, +2.0/3.0,
                                             +1.0/6.0, +1.0/6.0, +2.0/3.0};
    static constexpr std::array<scalar, 6> z{-std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0), -std::sqrt(1.0/3.0),
                                             +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0), +std::sqrt(1.0/3.0)};
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
    static constexpr std::array<scalar, 4> w{1.0/24.0, 1.0/24.0, 1.0/24.0, 1.0/24.0};
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

    virtual const vector& normal(label i) const = 0;

    virtual scalar weight(label i) const = 0;
    
    virtual const vector& at(label i) const = 0;
};

template<int nVertex, int nPoints>
class Quadrature: public Face
{
public:
    // Constructors - default
    Quadrature(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Quadrature() {};

    // Member Functions
    label size() const { return nPoints; };

    const vector& normal(label i) const { return normal_; }

    scalar weight(label i) const { return weights[i]; }
    
    const vector& at(label i) const { return quadPoints[i]; }

protected:
    static const std::array<scalar, nPoints> local_x_;
    static const std::array<scalar, nPoints> local_y_;
    static const std::array<scalar, nPoints> local_w_;
    vector normal_;
    std::vector<scalar> weights;
    std::vector<vector> quadPoints;
};

template<int nVertex, int nPoints>
const std::array<scalar, nPoints> Quadrature<nVertex, nPoints>::local_x_
    = Gauss<nVertex, nPoints>::buildLocalXCoords();

template<int nVertex, int nPoints>
const std::array<scalar, nPoints> Quadrature<nVertex, nPoints>::local_y_
    = Gauss<nVertex, nPoints>::buildLocalYCoords();

template<int nVertex, int nPoints>
const std::array<scalar, nPoints> Quadrature<nVertex, nPoints>::local_w_
    = Gauss<nVertex, nPoints>::buildLocalWeights();

template<int nVertex>
void GetCoordAndWeight
(
    const scalar x_local,
    const scalar y_local,
    const scalar w_local,
    const Eigen::Matrix<scalar, 3, nVertex>& facePoints,
    scalar& weight,
    vector& quadPoint
);

void build2ndFace(const fvMesh& mesh, const label& faceI, std::unique_ptr<Face>& face);

void build4stFace(const fvMesh& mesh, const label& faceI, std::unique_ptr<Face>& face);

void gaussHexa8
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weight,
    std::vector<vector>& quadPoints
);

void gaussPrism6
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weight,
    std::vector<vector>& quadPoints
);

void gaussTetra4
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weight,
    std::vector<vector>& quadPoints
);

void build2ndCell
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weights,
    std::vector<vector>& quadPoints
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //