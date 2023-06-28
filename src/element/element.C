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

template<int nVertex, int nPoints>
Foam::Quadrature<nVertex, nPoints>::Quadrature
(
    const fvMesh& mesh,
    const label& faceI
)
:
    normal_(mesh.Sf()[faceI]/mesh.magSf()[faceI])
{
    const UList<label>& facePointsId = mesh.faces()[faceI];
    Eigen::Matrix<scalar, 3, nVertex> facePoints;
    forAll(facePointsId, vertexI)
    {
        facePoints.col(vertexI) = Col3X1{mesh.points()[facePointsId[vertexI]][0],
                                         mesh.points()[facePointsId[vertexI]][1],
                                         mesh.points()[facePointsId[vertexI]][2]};
    }
    for (int i = 0; i != nPoints; ++i)
    {
        Surface<nVertex>::GetCoordAndWeight(local_x_[i], local_y_[i], local_w_[i],
                                            facePoints, weights[i], quadPoints[i]);
    }
}

template<int nVertex, int nPoints>
Foam::CurvedQuadrature<nVertex, nPoints>::CurvedQuadrature
(
    const fvMesh& mesh,
    const label& patchI,
    const label& faceI
)
:
    Quadrature<nVertex, nPoints>()
{
    Base::normal_ = mesh.Sf()[faceI]/mesh.magSf()[faceI];
    const label start = mesh.boundary()[patchI].start();
    const UList<label>& facePointsId = mesh.faces()[faceI];
    Eigen::Matrix<scalar, 3, nVertex> faceVertexes;
    Eigen::Matrix<scalar, 3, nVertex> vertexNormal;
    forAll(facePointsId, vertexI)
    {
        faceVertexes.col(vertexI) = Col3X1{mesh.points()[facePointsId[vertexI]][0],
                                           mesh.points()[facePointsId[vertexI]][1],
                                           mesh.points()[facePointsId[vertexI]][2]};
        const UList<label>& faceList = mesh.pointFaces()[facePointsId[vertexI]];
        vector tempVector = vector::zero;
        scalar tempScalar = 0;
        forAll(faceList, f)
        {
            if ((faceList[f] >= start) && (faceList[f] < start+mesh.boundary()[patchI].size()))
            {
                if ((Base::normal_&mesh.Sf()[faceList[f]]) <= 0) continue;
                tempVector += mesh.Sf()[faceList[f]];
                tempScalar += mesh.magSf()[faceList[f]];
            }
        }
        tempVector /= tempScalar;
        vertexNormal.col(vertexI) = Col3X1{tempVector[0], tempVector[1], tempVector[2]};
    }
    solveCurvedFace(faceVertexes, vertexNormal);
}


template<int nVertex, int nPoints>
void Foam::CurvedQuadrature<nVertex, nPoints>::solveCurvedFace
(
    const Eigen::Matrix<scalar, 3, nVertex>& faceVertexes,
    const Eigen::Matrix<scalar, 3, nVertex>& vertexNormal
)
{
    // Solving coefficients of the polynomials
    Col3X1 faceNormal{Base::normal_[0], Base::normal_[1], Base::normal_[2]};
    Col b = Surface<nVertex>::buildColumebForCurve(faceNormal, faceVertexes, vertexNormal);
    Col coefList = A_.colPivHouseholderQr().solve(b);
    for (int i = 0; i != nPoints; ++i)
    {
        // Solving quadrature points and weights for linear face
        Surface<nVertex>::GetCoordAndWeight(Base::local_x_[i], Base::local_y_[i], Base::local_w_[i],
                                            faceVertexes, Base::weights[i], Base::quadPoints[i]);
        // Solving new quadrature points
        const scalar h = Surface<nVertex>::PolyArray(Base::local_x_[i], Base::local_y_[i]).dot(coefList);
        Base::quadPoints[i] += Base::normal_*h;
        // Solving new quadrature normals
        Mat3X3 quadA;
        quadA.block<3,2>(0,0) = Surface<nVertex>::Jacobian(Base::local_x_[i], Base::local_y_[i], faceVertexes);
        quadA.col(2) = faceNormal;
        quadA.transposeInPlace();
        Col3X1 quadb{-Surface<nVertex>::PolyArrayDx(Base::local_x_[i], Base::local_y_[i]).dot(coefList),
                     -Surface<nVertex>::PolyArrayDy(Base::local_x_[i], Base::local_y_[i]).dot(coefList), 1};
        Col3X1 quadNormal = quadA.colPivHouseholderQr().solve(quadb);
        quadNormal /= Foam::sqrt(quadNormal.dot(quadNormal));
        normalList_[i] = vector(quadNormal(0), quadNormal(1), quadNormal(2));
        // Solving new quadrature weights
        Base::weights[i] /= (quadNormal.dot(faceNormal));
    }
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
        Arr8X1 factor_x = Local::x_hexa_i * Hexa8::x[i]; factor_x += 1;
        Arr8X1 factor_y = Local::y_hexa_i * Hexa8::y[i]; factor_y += 1;
        Arr8X1 factor_z = Local::z_hexa_i * Hexa8::z[i]; factor_z += 1;
        dn.col(0) << Local::x_hexa_i * factor_y * factor_z;
        dn.col(1) << Local::y_hexa_i * factor_x * factor_z;
        dn.col(2) << Local::z_hexa_i * factor_x * factor_y;
        dn *= 0.125;
        Mat3X3 dr = cellShapes.transpose()*dn;
        weight[i] = Hexa8::w[i] * mag(dr.determinant());
        // coord
        Arr8X1 N = 0.125*(1+Local::x_hexa_i*Hexa8::x[i])*(1+Local::y_hexa_i*Hexa8::y[i])*(1+Local::z_hexa_i*Hexa8::z[i]);
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

void Foam::build2ndFace
(
  const fvMesh& mesh,
  const label& faceI,
  std::unique_ptr<Face>& face
)
{
    const label nNodes = mesh.faces()[faceI].size();
    if (nNodes == 3) face.reset(new Quadrature<3,3>(mesh, faceI));
    if (nNodes == 4) face.reset(new Quadrature<4,4>(mesh, faceI));
}

void Foam::build4stFace
(
  const fvMesh& mesh,
  const label& faceI,
  std::unique_ptr<Face>& face
)
{
    const label nNodes = mesh.faces()[faceI].size();
    if (nNodes == 3) face.reset(new Quadrature<3,6>(mesh, faceI));
    if (nNodes == 4) face.reset(new Quadrature<4,9>(mesh, faceI));
}

void Foam::build2ndCurvedBoundary
(
  const fvMesh& mesh,
  const label& patchI,
  const label& faceI,
  std::unique_ptr<Face>& face
)
{
    const label nNodes = mesh.faces()[faceI].size();
    if (nNodes == 3) face.reset(new CurvedQuadrature<3,3>(mesh, patchI, faceI));
    if (nNodes == 4) face.reset(new CurvedQuadrature<4,4>(mesh, patchI, faceI));
}

void Foam::build4stCurvedBoundary
(
  const fvMesh& mesh,
  const label& patchI,
  const label& faceI,
  std::unique_ptr<Face>& face
)
{
    const label nNodes = mesh.faces()[faceI].size();
    if (nNodes == 3) face.reset(new CurvedQuadrature<3,6>(mesh, patchI, faceI));
    if (nNodes == 4) face.reset(new CurvedQuadrature<4,9>(mesh, patchI, faceI));
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

Foam::Cell::Cell
(
    const fvMesh& mesh,
    const label& cellI,
    const scalar& yWall,
    const vector& nWall
)
{
    const label nNodes = mesh.cellShapes()[cellI].size();
    if (nNodes == 8) gaussHexa8(mesh, cellI, weights, quadPoints);
    else if (nNodes == 6) gaussPrism6(mesh, cellI, weights, quadPoints);
    else if (nNodes == 4) gaussTetra4(mesh, cellI, weights, quadPoints);
    else Info << "Wrong element type!" << endl;
    dists.resize(nNodes);
    for (label i = 0; i != nNodes; ++i)
    {
        dists[i] = yWall - ((quadPoints[i]-mesh.C()[cellI])&nWall);
    }
}