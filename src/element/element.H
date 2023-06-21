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

#include "definitions.H"

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

/*---------------------------------------------------------------------------*\
                           namespace Declaration
\*---------------------------------------------------------------------------*/

namespace Local
{
    static const Arr4X1 x_quad_i{-1, +1, +1, -1};
    static const Arr4X1 y_quad_i{-1, -1, +1, +1};
    
    static const Arr8X1 x_hexa_i{-1, +1, +1, -1, -1, +1, +1, -1};
    static const Arr8X1 y_hexa_i{-1, -1, +1, +1, -1, -1, +1, +1};
    static const Arr8X1 z_hexa_i{-1, -1, -1, -1, +1, +1, +1, +1};
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
                           Class Surface Declaration
\*---------------------------------------------------------------------------*/

template<int nVertex>
class Surface;

template<>
class Surface<3>
{
public:

    static void GetCoordAndWeight
    (
        const scalar& x_local,
        const scalar& y_local,
        const scalar& w_local,
        const Eigen::Matrix<scalar, 3, 3>& facePoints,
        scalar& weight,
        vector& quadPoint
    )
    {
        // weight
        Mat3X2 dr = Jacobian(x_local, y_local, facePoints);
        weight = w_local* sqrt((dr.transpose() * dr).determinant());
        // coord
        Col3X1 N{x_local, y_local, 1-x_local-y_local};
        Col3X1 coord = facePoints*N;
        quadPoint = vector(coord(0), coord(1), coord(2));
    }

    static const Mat3X2 Jacobian
    (
        const scalar& x_local,
        const scalar& y_local,
        const Eigen::Matrix<scalar, 3, 3>& facePoints
    )
    {
        return facePoints*Mat3X2{{ 1, 0}, { 0, 1}, {-1,-1}};
    }

    static const auto PolyArray(const scalar& x, const scalar& y)
    {
        return Eigen::Matrix<scalar, 9, 1>{1, x, y, x*x, x*y, y*y, x*x*x, x*x*y, y*y*y};
    }

    static const auto PolyArrayDx(const scalar& x, const scalar& y)
    {
        return Eigen::Matrix<scalar, 9, 1>{0, 1, 0, 2*x, y, 0, 3*x*x, 2*x*y, 0};
    }

    static const auto PolyArrayDy(const scalar& x, const scalar& y)
    {
        return Eigen::Matrix<scalar, 9, 1>{0, 0, 1, 0, x, 2*y, 0, x*x, 3*y*y};
    }

    static const auto buildMatrixAForCurve()
    {
        Eigen::Matrix<scalar, 9, 9> A = Eigen::Matrix<scalar, 9, 9>::Zero();
        A.col(0) -= PolyArray(1, 0); A.col(1) -= PolyArrayDx(1, 0); A.col(2) -= PolyArrayDy(1, 0);
        A.col(3) -= PolyArray(0, 1); A.col(4) -= PolyArrayDx(0, 1); A.col(5) -= PolyArrayDy(0, 1);
        A.col(6) -= PolyArray(0, 0); A.col(7) -= PolyArrayDx(0, 0); A.col(8) -= PolyArrayDy(0, 0);
        A.transposeInPlace();
        return A;
    }

    static const auto buildColumebForCurve
    (
        const Col3X1& faceNormal,
        const Eigen::Matrix<scalar, 3, 3>& faceVertexes,
        const Eigen::Matrix<scalar, 3, 3>& vertexNormal
    )
    {
        Eigen::Matrix<scalar, 9, 1> b = Eigen::Matrix<scalar, 9, 1>::Zero();
        Mat3X2 jacobian = Jacobian(0, 0, faceVertexes);
        for (int i = 0; i != 3; i++)
        {
            Col3X1 F_xyz = vertexNormal.col(i)/(vertexNormal.col(i).dot(faceNormal));
            b.block<2,1>(i*3+1,0) = jacobian.transpose()*F_xyz;
        }
        return b;
    }
};

template<>
class Surface<4>
{
public:

    static void GetCoordAndWeight
    (
        const scalar& x_local,
        const scalar& y_local,
        const scalar& w_local,
        const Eigen::Matrix<scalar, 3, 4>& facePoints,
        scalar& weight,
        vector& quadPoint
    )
    {
        // weight
        Mat3X2 dr = Jacobian(x_local, y_local, facePoints);
        weight = w_local* sqrt((dr.transpose() * dr).determinant());
        // coord
        Arr4X1 N = 0.25*(1+Local::x_quad_i*x_local)*(1+Local::y_quad_i*y_local);
        Col3X1 coord = facePoints*N.matrix();
        quadPoint = vector(coord(0), coord(1), coord(2));
    }

    static const Mat3X2 Jacobian
    (
        const scalar& x_local,
        const scalar& y_local,
        const Eigen::Matrix<scalar, 3, 4>& facePoints
    )
    {
        // weight
        Mat4X2 dn;
        Arr4X1 factor_x = Local::x_quad_i * x_local; factor_x += 1;
        Arr4X1 factor_y = Local::y_quad_i * y_local; factor_y += 1;
        dn.col(0) << Local::x_quad_i * factor_y;
        dn.col(1) << Local::y_quad_i * factor_x;
        dn *= 0.25;
        return facePoints*dn;
    }

    static const auto PolyArray(const scalar& x, const scalar& y)
    {
        return Eigen::Matrix<scalar, 12, 1>{1, x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y, x*x*x*x, y*y*y*y};
    }

    static const auto PolyArrayDx(const scalar& x, const scalar& y)
    {
        return Eigen::Matrix<scalar, 12, 1>{0, 1, 0, 2*x, y, 0, 3*x*x, 2*x*y, y*y, 0, 4*x*x*x, 0};
    }

    static const auto PolyArrayDy(const scalar& x, const scalar& y)
    {
        return Eigen::Matrix<scalar, 12, 1>{0, 0, 1, 0, x, 2*y, 0, x*x, 2*x*y, 3*y*y, 0, 4*y*y*y};
    }

    static const auto buildMatrixAForCurve()
    {
        Eigen::Matrix<scalar, 12, 12> A = Eigen::Matrix<scalar, 12, 12>::Zero();
        A.col(0) -= PolyArray(-1, -1); A.col(1) -= PolyArrayDx(-1, -1); A.col(2) -= PolyArrayDy(-1, -1);
        A.col(3) -= PolyArray( 1, -1); A.col(4) -= PolyArrayDx( 1, -1); A.col(5) -= PolyArrayDy( 1, -1);
        A.col(6) -= PolyArray( 1,  1); A.col(7) -= PolyArrayDx( 1,  1); A.col(8) -= PolyArrayDy( 1,  1);
        A.col(9) -= PolyArray(-1,  1); A.col(10)-= PolyArrayDx(-1,  1); A.col(11)-= PolyArrayDy(-1,  1);
        A.transposeInPlace();
        return A;
    }

    static const auto buildColumebForCurve
    (
        const Col3X1& faceNormal,
        const Eigen::Matrix<scalar, 3, 4>& faceVertexes,
        const Eigen::Matrix<scalar, 3, 4>& vertexNormal
    )
    {
        Eigen::Matrix<scalar, 12, 1> b = Eigen::Matrix<scalar, 12, 1>::Zero();
        for (int i = 0; i != 4; i++)
        {
            Mat3X2 jacobian = Jacobian(Local::x_quad_i(i), Local::y_quad_i(i), faceVertexes);
            Col3X1 F_xyz = vertexNormal.col(i)/(vertexNormal.col(i).dot(faceNormal));
            b.block<2,1>(i*3+1,0) = jacobian.transpose()*F_xyz;
        }
        return b;
    }
};

/*---------------------------------------------------------------------------*\
                           Class Face Declaration
\*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*\
                           Class Quadrature Declaration
\*---------------------------------------------------------------------------*/

template<int nVertex, int nPoints>
class Quadrature: public Face
{
public:
    // Constructors - default
    Quadrature() = default;

    Quadrature(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Quadrature() {};

    // Member Functions
    label size() const { return nPoints; };

    virtual const vector& normal(label i) const { return normal_; }

    scalar weight(label i) const { return weights[i]; }
    
    const vector& at(label i) const { return quadPoints[i]; }

    static const std::array<scalar, nPoints> local_x_;
    static const std::array<scalar, nPoints> local_y_;
    static const std::array<scalar, nPoints> local_w_;
    vector normal_;
    std::array<scalar, nPoints> weights;
    std::array<vector, nPoints> quadPoints;
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

/*---------------------------------------------------------------------------*\
                           Class CurvedQuadrature Declaration
\*---------------------------------------------------------------------------*/

template<int nVertex, int nPoints>
class CurvedQuadrature: public Quadrature<nVertex, nPoints>
{
public:
    using Base = Quadrature<nVertex, nPoints>;
    using Col = Eigen::Matrix<scalar, nVertex*3, 1>;
    using Mat = Eigen::Matrix<scalar, nVertex*3, nVertex*3>;
    // Constructors - default
    CurvedQuadrature(const fvMesh& mesh, const label& patchI, const label& faceI);

    // Destructor - default
    virtual ~CurvedQuadrature() {};
    
    const vector& normal(label i) const { return normalList_[i]; }

    void solveCurvedFace
    (
        const Eigen::Matrix<scalar, 3, nVertex>& faceVertexes,
        const Eigen::Matrix<scalar, 3, nVertex>& vertexNormal
    );

private:
    static const Mat A_;
    std::array<vector, nPoints> normalList_;
};

template<int nVertex, int nPoints>
typename CurvedQuadrature<nVertex, nPoints>::Mat
const CurvedQuadrature<nVertex, nPoints>::A_
    = Surface<nVertex>::buildMatrixAForCurve();

/*---------------------------------------------------------------------------*\
                           Methods Declaration
\*---------------------------------------------------------------------------*/

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

void build2ndFace(const fvMesh& mesh, const label& faceI, std::unique_ptr<Face>& face);

void build4stFace(const fvMesh& mesh, const label& faceI, std::unique_ptr<Face>& face);

void build2ndCurvedBoundary(const fvMesh& mesh, const label& patchI, const label& faceI, std::unique_ptr<Face>& face);

void build4stCurvedBoundary(const fvMesh& mesh, const label& patchI, const label& faceI, std::unique_ptr<Face>& face);

void build2ndCell
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weights,
    std::vector<vector>& quadPoints
);



/*---------------------------------------------------------------------------*\
                           Class Cell Declaration
\*---------------------------------------------------------------------------*/

class Cell
{
public:

    // Destructor - default
    virtual ~Cell() {};
    
    // Constructors - default
    Cell() = default;

    Cell
    (
        const fvMesh& mesh,
        const label& cellI, 
        const scalar& yWall,
        const vector& nWall
    );

    // Member Functions
    label size() const { return weights.size(); };

    const scalar& distToWall(label i) const { return dists[i]; };

    scalar weight(label i) const { return weights[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };
    
private:
    std::vector<scalar> weights;
    std::vector<scalar> dists;
    std::vector<vector> quadPoints;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //