void Foam::gaussQuad4
(
    const fvMesh& mesh,
    const label& faceI,
    std::vector<vector>& quadPoints
)
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

void Foam::gaussQuad9
(
    const fvMesh& mesh,
    const label& faceI,
    std::vector<vector>& quadPoints
)
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

void Foam::gaussHexa8
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<vector>& quadPoints
)
{
    const UList<label>& cellShapesId = mesh.cellShapes()[cellI];
    std::array<vector, 8> cellShapes{mesh.points()[cellShapesId[0]], mesh.points()[cellShapesId[1]],
                                     mesh.points()[cellShapesId[2]], mesh.points()[cellShapesId[3]],
                                     mesh.points()[cellShapesId[4]], mesh.points()[cellShapesId[5]],
                                     mesh.points()[cellShapesId[6]], mesh.points()[cellShapesId[7]]};
    for (label i = 0; i != 8; ++i)
    {
        std::array<scalar, 8> N{0.125*(1-Hexa8::x[i])*(1-Hexa8::y[i])*(1+Hexa8::z[i]),
                                0.125*(1+Hexa8::x[i])*(1-Hexa8::y[i])*(1+Hexa8::z[i]),
                                0.125*(1+Hexa8::x[i])*(1+Hexa8::y[i])*(1+Hexa8::z[i]),
                                0.125*(1-Hexa8::x[i])*(1+Hexa8::y[i])*(1+Hexa8::z[i]),
                                0.125*(1-Hexa8::x[i])*(1-Hexa8::y[i])*(1-Hexa8::z[i]),
                                0.125*(1+Hexa8::x[i])*(1-Hexa8::y[i])*(1-Hexa8::z[i]),
                                0.125*(1+Hexa8::x[i])*(1+Hexa8::y[i])*(1-Hexa8::z[i]),
                                0.125*(1-Hexa8::x[i])*(1+Hexa8::y[i])*(1-Hexa8::z[i])};
        quadPoints[i] = N[0]*cellShapes[0] + N[1]*cellShapes[1] + N[2]*cellShapes[2] + N[3]*cellShapes[3] +
                        N[4]*cellShapes[4] + N[5]*cellShapes[5] + N[6]*cellShapes[6] + N[7]*cellShapes[7];
    }
}

void Foam::gaussPrism6
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<vector>& quadPoints
)
{
    const UList<label>& cellShapesId = mesh.cellShapes()[cellI];
    std::array<vector, 6> cellShapes{mesh.points()[cellShapesId[0]], mesh.points()[cellShapesId[1]],
                                     mesh.points()[cellShapesId[2]], mesh.points()[cellShapesId[3]],
                                     mesh.points()[cellShapesId[4]], mesh.points()[cellShapesId[5]]};
    for (label i = 0; i != 6; ++i)
    {
        std::array<scalar, 6> N{0.5*(1-Prism6::z[i])*Prism6::x[i],
                                0.5*(1-Prism6::z[i])*Prism6::y[i],
                                0.5*(1-Prism6::z[i])*(1-Prism6::x[i]-Prism6::y[i]),
                                0.5*(1+Prism6::z[i])*Prism6::x[i],
                                0.5*(1+Prism6::z[i])*Prism6::y[i],
                                0.5*(1+Prism6::z[i])*(1-Prism6::x[i]-Prism6::y[i])};
        quadPoints[i] = N[0]*cellShapes[0] + N[1]*cellShapes[1] + N[2]*cellShapes[2] +
                        N[3]*cellShapes[3] + N[4]*cellShapes[4] + N[5]*cellShapes[5];
    }
}

void Foam::gaussTetra4
(
    const fvMesh& mesh,
    const label& cellI,
    std::vector<vector>& quadPoints
)
{
    const UList<label>& cellShapesId = mesh.cellShapes()[cellI];
    std::array<vector, 4> cellShapes{mesh.points()[cellShapesId[0]], mesh.points()[cellShapesId[1]],
                                     mesh.points()[cellShapesId[2]], mesh.points()[cellShapesId[3]]};
    for (label i = 0; i != 4; ++i)
    {
        quadPoints[i] = cellShapes[0]*Tetra4::x[i] + cellShapes[1]*Tetra4::y[i] +
                        cellShapes[2]*Tetra4::z[i] + cellShapes[3]*(1.0-Tetra4::x[i]-Tetra4::y[i]-Tetra4::z[i]);
    }
}