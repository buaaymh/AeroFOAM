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

Foam::Triangle3::Triangle3
(
    const fvMesh& mesh,
    const label& faceI
)
:
    quadPoints(3)
{
    const UList<label>& facePointsId = mesh.faces()[faceI];
    std::array<vector, 3> facePoints{mesh.points()[facePointsId[0]], mesh.points()[facePointsId[1]], mesh.points()[facePointsId[2]]};
    tensor transform
    (
      facePoints[0].x(), facePoints[1].x(), facePoints[2].x(),
      facePoints[0].y(), facePoints[1].y(), facePoints[2].y(),
      facePoints[0].z(), facePoints[1].z(), facePoints[2].z()
    );
    for (label i = 0; i != 3; ++i)
    {
        quadPoints[i] = transform&vector(Tria3::x[i], Tria3::y[i], 1-Tria3::x[i]-Tria3::y[i]);
    }
}

Foam::Triangle4::Triangle4
(
    const fvMesh& mesh,
    const label& faceI
)
:
    quadPoints(4)
{
    const UList<label>& facePointsId = mesh.faces()[faceI];
    std::array<vector, 3> facePoints{mesh.points()[facePointsId[0]], mesh.points()[facePointsId[1]], mesh.points()[facePointsId[2]]};
    tensor transform
    (
      facePoints[0].x(), facePoints[1].x(), facePoints[2].x(),
      facePoints[0].y(), facePoints[1].y(), facePoints[2].y(),
      facePoints[0].z(), facePoints[1].z(), facePoints[2].z()
    );
    for (label i = 0; i != 4; ++i)
    {
        quadPoints[i] = transform&vector(Tria4::x[i], Tria4::y[i], 1-Tria4::x[i]-Tria4::y[i]);
    }
}

Foam::Quadrangle4::Quadrangle4
(
    const fvMesh& mesh,
    const label& faceI
)
:
    quadPoints(4)
{
    const UList<label>& facePointsId = mesh.faces()[faceI];
    std::array<vector, 4> facePoints{mesh.points()[facePointsId[0]], mesh.points()[facePointsId[1]],
                                     mesh.points()[facePointsId[2]], mesh.points()[facePointsId[3]]};
    for (label i = 0; i != 4; ++i)
    {
        std::array<scalar, 4> N{0.25*(1+Quad4::x[i])*(1+Quad4::y[i]), 0.25*(1-Quad4::x[i])*(1+Quad4::y[i]),
                                0.25*(1-Quad4::x[i])*(1-Quad4::y[i]), 0.25*(1+Quad4::x[i])*(1-Quad4::y[i])};
        quadPoints[i] = N[0]*facePoints[0] + N[1]*facePoints[1] + N[2]*facePoints[2] + N[3]*facePoints[3];
    }
}

Foam::Quadrangle9::Quadrangle9
(
    const fvMesh& mesh,
    const label& faceI
)
:
    quadPoints(9)
{
    const UList<label>& facePointsId = mesh.faces()[faceI];
    std::array<vector, 4> facePoints{mesh.points()[facePointsId[0]], mesh.points()[facePointsId[1]],
                                     mesh.points()[facePointsId[2]], mesh.points()[facePointsId[3]]};
    for (label i = 0; i != 9; ++i)
    {
        std::array<scalar, 4> N{0.25*(1+Quad9::x[i])*(1+Quad9::y[i]), 0.25*(1-Quad9::x[i])*(1+Quad9::y[i]),
                                0.25*(1-Quad9::x[i])*(1-Quad9::y[i]), 0.25*(1+Quad9::x[i])*(1-Quad9::y[i])};
        quadPoints[i] = N[0]*facePoints[0] + N[1]*facePoints[1] + N[2]*facePoints[2] + N[3]*facePoints[3];
    }
}