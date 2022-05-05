void Foam::vrScheme::gaussQuad4
(
    label faceI,
    std::vector<vector>& quadPoints
)
{
    const UList<label>& facePointsId = mesh_.faces()[faceI];
    std::array<vector, 4> facePoints{mesh_.points()[facePointsId[0]], mesh_.points()[facePointsId[1]],
                                     mesh_.points()[facePointsId[2]], mesh_.points()[facePointsId[3]]};
    for (label i = 0; i != 4; ++i)
    {
        std::array<scalar, 4> N{0.25*(1+Quad4::x[i])*(1+Quad4::y[i]), 0.25*(1-Quad4::x[i])*(1+Quad4::y[i]),
                                0.25*(1-Quad4::x[i])*(1-Quad4::y[i]), 0.25*(1+Quad4::x[i])*(1-Quad4::y[i])};
        quadPoints[i] = N[0]*facePoints[0] + N[1]*facePoints[1] + N[2]*facePoints[2] + N[3]*facePoints[3];
    }
}

void Foam::vrScheme::gaussQuad9
(
    label faceI,
    std::vector<vector>& quadPoints
)
{
    const UList<label>& facePointsId = mesh_.faces()[faceI];
    std::array<vector, 4> facePoints{mesh_.points()[facePointsId[0]], mesh_.points()[facePointsId[1]],
                                     mesh_.points()[facePointsId[2]], mesh_.points()[facePointsId[3]]};
    for (label i = 0; i != 9; ++i)
    {
        std::array<scalar, 4> N{0.25*(1+Quad9::x[i])*(1+Quad9::y[i]), 0.25*(1-Quad9::x[i])*(1+Quad9::y[i]),
                                0.25*(1-Quad9::x[i])*(1-Quad9::y[i]), 0.25*(1+Quad9::x[i])*(1-Quad9::y[i])};
        quadPoints[i] = N[0]*facePoints[0] + N[1]*facePoints[1] + N[2]*facePoints[2] + N[3]*facePoints[3];
    }
}

void Foam::vrScheme::gaussHex8
(
    label cellI,
    std::vector<vector>& quadPoints
)
{
    const UList<label>& cellPointsId = mesh_.cellPoints()[cellI];
    std::array<vector, 8> cellPoints{mesh_.points()[cellPointsId[0]], mesh_.points()[cellPointsId[1]],
                                     mesh_.points()[cellPointsId[2]], mesh_.points()[cellPointsId[3]],
                                     mesh_.points()[cellPointsId[4]], mesh_.points()[cellPointsId[5]],
                                     mesh_.points()[cellPointsId[6]], mesh_.points()[cellPointsId[7]]};
    for (label i = 0; i != 8; ++i)
    {
        std::array<scalar, 8> N{0.125*(1-Hex8::x[i])*(1+Hex8::y[i])*(1+Hex8::z[i]),
                                0.125*(1+Hex8::x[i])*(1+Hex8::y[i])*(1+Hex8::z[i]),
                                0.125*(1-Hex8::x[i])*(1-Hex8::y[i])*(1+Hex8::z[i]),
                                0.125*(1+Hex8::x[i])*(1-Hex8::y[i])*(1+Hex8::z[i]),
                                0.125*(1-Hex8::x[i])*(1+Hex8::y[i])*(1-Hex8::z[i]),
                                0.125*(1+Hex8::x[i])*(1+Hex8::y[i])*(1-Hex8::z[i]),
                                0.125*(1-Hex8::x[i])*(1-Hex8::y[i])*(1-Hex8::z[i]),
                                0.125*(1+Hex8::x[i])*(1-Hex8::y[i])*(1-Hex8::z[i])};
        quadPoints[i] = N[0]*cellPoints[0] + N[1]*cellPoints[1] + N[2]*cellPoints[2] + N[3]*cellPoints[3] +
                        N[4]*cellPoints[4] + N[5]*cellPoints[5] + N[6]*cellPoints[6] + N[7]*cellPoints[7];
    }
}