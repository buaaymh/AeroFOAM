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
    static std::array<scalar, nPoints> x;
    static std::array<scalar, nPoints> y;
    static std::array<scalar, nPoints> w;
    vector normal_;
    std::vector<scalar> weights;
    std::vector<vector> quadPoints;
};

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