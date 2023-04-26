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

#include "element.H"

Foam::Face::Face
(
    const fvMesh& mesh,
    const label& faceI
)
:
    normal_(mesh.Sf()[faceI]/mesh.magSf()[faceI]),
    weights(6),
    quadPoints(6)
{}

void Foam::Face::GetQuadrangleCoordAndWeight
(
    const scalar x_local,
    const scalar y_local,
    const scalar w_local,
    const Mat4X3& facePoints,
    scalar& weight,
    vector& quadPoint
)
{
    // weight
    Mat4X2 dn;
    Arr4X1 factor_x = Quad4::x_hexa_i * x_local; factor_x += 1;
    Arr4X1 factor_y = Quad4::y_hexa_i * y_local; factor_y += 1;
    dn.col(0) << Quad4::x_hexa_i * factor_y;
    dn.col(1) << Quad4::y_hexa_i * factor_x;
    dn *= 0.25;
    Mat3X2 dr = facePoints.transpose()*dn;
    weight = w_local* sqrt((dr.transpose() * dr).determinant());
    // coord
    Arr4X1 N = 0.25*(1+Quad4::x_hexa_i*x_local)*(1+Quad4::y_hexa_i*y_local);
    Col3X1 coord = facePoints.transpose()*N.matrix();
    quadPoint = vector(coord(0), coord(1), coord(2));
}

void Foam::Face::GetTriangleCoordAndWeight
(
    const scalar x_local,
    const scalar y_local,
    const scalar w_local,
    const Mat3X3& facePoints,
    scalar& weight,
    vector& quadPoint
)
{
    // weight
    Mat3X2 dn
    {
        { 1,  0},
        { 0,  1},
        {-1, -1}
    };
    Mat3X2 dr = facePoints.transpose()*dn;
    weight = w_local* sqrt((dr.transpose() * dr).determinant());
    // coord
    Col3X1 N{x_local, y_local, 1-x_local-y_local};
    Col3X1 coord = facePoints.transpose()*N;
    quadPoint = vector(coord(0), coord(1), coord(2));
}

Foam::Triangle3::Triangle3
(
    const fvMesh& mesh,
    const label& faceI
)
:
    Face(mesh, faceI)
{
    weights.resize(3);
    quadPoints.resize(3);
    const UList<label>& facePointsId = mesh.faces()[faceI];
    Mat3X3 facePoints
    {
        {mesh.points()[facePointsId[0]][0], mesh.points()[facePointsId[0]][1], mesh.points()[facePointsId[0]][2]},
        {mesh.points()[facePointsId[1]][0], mesh.points()[facePointsId[1]][1], mesh.points()[facePointsId[1]][2]},
        {mesh.points()[facePointsId[2]][0], mesh.points()[facePointsId[2]][1], mesh.points()[facePointsId[2]][2]}
    };
    for (label i = 0; i != 3; ++i)
    {
        GetTriangleCoordAndWeight(Tria3::x[i], Tria3::y[i], Tria3::w[i], facePoints,
                                  weights[i], quadPoints[i]);
    }
}

Foam::Triangle6::Triangle6
(
    const fvMesh& mesh,
    const label& faceI
)
:
    Face(mesh, faceI)
{
    weights.resize(6);
    quadPoints.resize(6);
    const UList<label>& facePointsId = mesh.faces()[faceI];
    Mat3X3 facePoints
    {
        {mesh.points()[facePointsId[0]][0], mesh.points()[facePointsId[0]][1], mesh.points()[facePointsId[0]][2]},
        {mesh.points()[facePointsId[1]][0], mesh.points()[facePointsId[1]][1], mesh.points()[facePointsId[1]][2]},
        {mesh.points()[facePointsId[2]][0], mesh.points()[facePointsId[2]][1], mesh.points()[facePointsId[2]][2]}
    };
    for (label i = 0; i != 6; ++i)
    {
        GetTriangleCoordAndWeight(Tria6::x[i], Tria6::y[i], Tria6::w[i], facePoints,
                                  weights[i], quadPoints[i]);
    }
}

Foam::Quadrangle4::Quadrangle4
(
    const fvMesh& mesh,
    const label& faceI
)
:
    Face(mesh, faceI)
{
    weights.resize(4);
    quadPoints.resize(4);
    const UList<label>& facePointsId = mesh.faces()[faceI];
    Mat4X3 facePoints
    {
        {mesh.points()[facePointsId[0]][0], mesh.points()[facePointsId[0]][1], mesh.points()[facePointsId[0]][2]},
        {mesh.points()[facePointsId[1]][0], mesh.points()[facePointsId[1]][1], mesh.points()[facePointsId[1]][2]},
        {mesh.points()[facePointsId[2]][0], mesh.points()[facePointsId[2]][1], mesh.points()[facePointsId[2]][2]},
        {mesh.points()[facePointsId[3]][0], mesh.points()[facePointsId[3]][1], mesh.points()[facePointsId[3]][2]}
    };
    for (label i = 0; i != 4; ++i)
    {
        GetQuadrangleCoordAndWeight(Quad4::x[i], Quad4::y[i], Quad4::w[i], facePoints,
                                    weights[i], quadPoints[i]);
    }
}

Foam::Quadrangle9::Quadrangle9
(
    const fvMesh& mesh,
    const label& faceI
)
:
    Face(mesh, faceI)
{
    weights.resize(9);
    quadPoints.resize(9);
    const UList<label>& facePointsId = mesh.faces()[faceI];
    Mat4X3 facePoints
    {
        {mesh.points()[facePointsId[0]][0], mesh.points()[facePointsId[0]][1], mesh.points()[facePointsId[0]][2]},
        {mesh.points()[facePointsId[1]][0], mesh.points()[facePointsId[1]][1], mesh.points()[facePointsId[1]][2]},
        {mesh.points()[facePointsId[2]][0], mesh.points()[facePointsId[2]][1], mesh.points()[facePointsId[2]][2]},
        {mesh.points()[facePointsId[3]][0], mesh.points()[facePointsId[3]][1], mesh.points()[facePointsId[3]][2]}
    };
    for (label i = 0; i != 9; ++i)
    {
        GetQuadrangleCoordAndWeight(Quad9::x[i], Quad9::y[i], Quad9::w[i], facePoints,
                                    weights[i], quadPoints[i]);
    }
}

void Foam::build2ndFace
(
  const fvMesh& mesh,
  const label& faceI,
  std::unique_ptr<Face>& face
)
{
    const UList<label>& facePointsId = mesh.faces()[faceI];
    const label nNodes = facePointsId.size();
    if (nNodes == 3) face.reset(new Triangle3(mesh, faceI));
    if (nNodes == 4) face.reset(new Quadrangle4(mesh, faceI));
}

void Foam::build4stFace
(
  const fvMesh& mesh,
  const label& faceI,
  std::unique_ptr<Face>& face
)
{
    const UList<label>& facePointsId = mesh.faces()[faceI];
    const label nNodes = facePointsId.size();
    if (nNodes == 3) face.reset(new Triangle6(mesh, faceI));
    if (nNodes == 4) face.reset(new Quadrangle9(mesh, faceI));
}

void Foam::gaussHexa8
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weight,
    std::vector<vector>& quadPoints
)
{
    weight.resize(8);
    quadPoints.resize(8);
    const UList<label>& cellShapesId = mesh.cellShapes()[cellI];
    Mat8X3 cellShapes
    {
        {mesh.points()[cellShapesId[0]][0], mesh.points()[cellShapesId[0]][1], mesh.points()[cellShapesId[0]][2]},
        {mesh.points()[cellShapesId[1]][0], mesh.points()[cellShapesId[1]][1], mesh.points()[cellShapesId[1]][2]},
        {mesh.points()[cellShapesId[2]][0], mesh.points()[cellShapesId[2]][1], mesh.points()[cellShapesId[2]][2]},
        {mesh.points()[cellShapesId[3]][0], mesh.points()[cellShapesId[3]][1], mesh.points()[cellShapesId[3]][2]},
        {mesh.points()[cellShapesId[4]][0], mesh.points()[cellShapesId[4]][1], mesh.points()[cellShapesId[4]][2]},
        {mesh.points()[cellShapesId[5]][0], mesh.points()[cellShapesId[5]][1], mesh.points()[cellShapesId[5]][2]},
        {mesh.points()[cellShapesId[6]][0], mesh.points()[cellShapesId[6]][1], mesh.points()[cellShapesId[6]][2]},
        {mesh.points()[cellShapesId[7]][0], mesh.points()[cellShapesId[7]][1], mesh.points()[cellShapesId[7]][2]}
    };
    for (label i = 0; i != 8; ++i)
    {
        // weight
        Mat8X3 dn;
        Arr8X1 factor_x = Hexa8::x_hexa_i * Hexa8::x[i]; factor_x += 1;
        Arr8X1 factor_y = Hexa8::y_hexa_i * Hexa8::y[i]; factor_y += 1;
        Arr8X1 factor_z = Hexa8::z_hexa_i * Hexa8::z[i]; factor_z += 1;
        dn.col(0) << Hexa8::x_hexa_i * factor_y * factor_z;
        dn.col(1) << Hexa8::y_hexa_i * factor_x * factor_z;
        dn.col(2) << Hexa8::z_hexa_i * factor_x * factor_y;
        dn *= 0.125;
        Mat3X3 dr = cellShapes.transpose()*dn;
        weight[i] = Hexa8::w[i] * mag(dr.determinant());
        // coord
        Arr8X1 N = 0.125*(1+Hexa8::x_hexa_i*Hexa8::x[i])*(1+Hexa8::y_hexa_i*Hexa8::y[i])*(1+Hexa8::z_hexa_i*Hexa8::z[i]);
        Col3X1 coord = cellShapes.transpose()*N.matrix();
        quadPoints[i] = vector(coord(0), coord(1), coord(2));
    }
}

void Foam::gaussPrism6
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weight,
    std::vector<vector>& quadPoints
)
{
    weight.resize(6);
    quadPoints.resize(6);
    const UList<label>& cellShapesId = mesh.cellShapes()[cellI];
    Mat6X3 cellShapes
    {
        {mesh.points()[cellShapesId[0]][0], mesh.points()[cellShapesId[0]][1], mesh.points()[cellShapesId[0]][2]},
        {mesh.points()[cellShapesId[1]][0], mesh.points()[cellShapesId[1]][1], mesh.points()[cellShapesId[1]][2]},
        {mesh.points()[cellShapesId[2]][0], mesh.points()[cellShapesId[2]][1], mesh.points()[cellShapesId[2]][2]},
        {mesh.points()[cellShapesId[3]][0], mesh.points()[cellShapesId[3]][1], mesh.points()[cellShapesId[3]][2]},
        {mesh.points()[cellShapesId[4]][0], mesh.points()[cellShapesId[4]][1], mesh.points()[cellShapesId[4]][2]},
        {mesh.points()[cellShapesId[5]][0], mesh.points()[cellShapesId[5]][1], mesh.points()[cellShapesId[5]][2]}
    };
    for (label i = 0; i != 6; ++i)
    {
        // weight
        Mat6X3 dn
        {
            { 1-Prism6::z[i],               0,               -Prism6::x[i]},
            {              0,  1-Prism6::z[i],               -Prism6::y[i]},
            { Prism6::z[i]-1,  Prism6::z[i]-1, Prism6::x[i]+Prism6::y[i]-1},
            { 1+Prism6::z[i],               0,                Prism6::x[i]},
            {              0,  1+Prism6::z[i],                Prism6::y[i]},
            {-Prism6::z[i]-1, -Prism6::z[i]-1, 1-Prism6::x[i]-Prism6::y[i]},
        };
        dn *= 0.5;
        Mat3X3 dr = cellShapes.transpose()*dn;
        weight[i] = Prism6::w[i] * mag(dr.determinant());
        // coord
        Col6X1 N{0.5*(1-Prism6::z[i])*Prism6::x[i],
                 0.5*(1-Prism6::z[i])*Prism6::y[i],
                 0.5*(1-Prism6::z[i])*(1-Prism6::x[i]-Prism6::y[i]),
                 0.5*(1+Prism6::z[i])*Prism6::x[i],
                 0.5*(1+Prism6::z[i])*Prism6::y[i],
                 0.5*(1+Prism6::z[i])*(1-Prism6::x[i]-Prism6::y[i])};
        Col3X1 coord = cellShapes.transpose()*N;
        quadPoints[i] = vector(coord(0), coord(1), coord(2));
    }
}

void Foam::gaussTetra4
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weight,
    std::vector<vector>& quadPoints
)
{
    weight.resize(4);
    quadPoints.resize(4);
    const UList<label>& cellShapesId = mesh.cellShapes()[cellI];
    Mat4X3 cellShapes
    {
        {mesh.points()[cellShapesId[0]][0], mesh.points()[cellShapesId[0]][1], mesh.points()[cellShapesId[0]][2]},
        {mesh.points()[cellShapesId[1]][0], mesh.points()[cellShapesId[1]][1], mesh.points()[cellShapesId[1]][2]},
        {mesh.points()[cellShapesId[2]][0], mesh.points()[cellShapesId[2]][1], mesh.points()[cellShapesId[2]][2]},
        {mesh.points()[cellShapesId[3]][0], mesh.points()[cellShapesId[3]][1], mesh.points()[cellShapesId[3]][2]}
    };
    for (label i = 0; i != 4; ++i)
    {
        // weight
        Mat4X3 dn
        {
            { 1,  0,  0},
            { 0,  1,  0},
            { 0,  0,  1},
            {-1, -1, -1},
        };
        Mat3X3 dr = cellShapes.transpose()*dn;
        weight[i] = Tetra4::w[i] * mag(dr.determinant());
        // coord
        Col4X1 N{Tetra4::x[i], Tetra4::y[i], Tetra4::z[i], 1-Tetra4::x[i]-Tetra4::y[i]-Tetra4::z[i]};
        Col3X1 coord = cellShapes.transpose()*N;
        quadPoints[i] = vector(coord(0), coord(1), coord(2));
    }
}

void Foam::build2ndCell
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<scalar>& weights,
    std::vector<vector>& quadPoints
)
{
    const label nNodes = mesh.cellShapes()[cellI].size();
    if (nNodes == 8) gaussHexa8(mesh, cellI, weights, quadPoints);
    else if (nNodes == 6) gaussPrism6(mesh, cellI, weights, quadPoints);
    else if (nNodes == 4) gaussTetra4(mesh, cellI, weights, quadPoints);
    else Info << "Wrong element type!" << endl;
}