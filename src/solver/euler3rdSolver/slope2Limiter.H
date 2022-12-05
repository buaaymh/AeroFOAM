std::vector<Mat9X5> coefsLimited(count);
count = 0;
forAll(mesh_.cells(), cellI)
{
    if (!isTrouble[cellI]) continue;
    const labelList& cellFaces = mesh_.cells()[cellI];
    std::vector<vector> normalList;
    std::vector<vector> deltaList;
    std::vector<Col5X1> d0_AdjList;
    std::vector<Mat5X3> d1_AdjList;
    normalList.reserve(mesh_.cellCells()[cellI].size());
    deltaList.reserve(mesh_.cellCells()[cellI].size());
    d0_AdjList.reserve(mesh_.cellCells()[cellI].size());
    d1_AdjList.reserve(mesh_.cellCells()[cellI].size());
    const Col5X1 d0_Own(rho_[cellI], rhoU_[cellI][0], rhoU_[cellI][1], rhoU_[cellI][2], rhoE_[cellI]);
    const Mat5X9 d1_d2_Own = coefs_[cellI].transpose();
    forAll(cellFaces, f)
    {
        const label faceI = cellFaces[f];
        if (mesh_.isInternalFace(faceI))
        {
            label j = owner[faceI];
            if (j == cellI) { j = neighbour[faceI]; }
            normalList.emplace_back(normal_[faceI]);
            deltaList.emplace_back(mesh_.Cf()[faceI]-mesh_.C()[cellI]);
            d0_AdjList.emplace_back(Col5X1(rho_[j], rhoU_[j][0], rhoU_[j][1], rhoU_[j][2], rhoE_[j]));
            Mat5X3 d1_Adj(coefs_[j].block<3,5>(0,0).transpose());
            d1_Adj.col(0) *= rDeltaXYZ_[j][0];
            d1_Adj.col(1) *= rDeltaXYZ_[j][1];
            d1_Adj.col(2) *= rDeltaXYZ_[j][2];
            d1_AdjList.emplace_back(d1_Adj);
        }
    }
    if (characteristic)
    {
        coefsLimited[count++] = slope2_Limiter(normalList, deltaList, d0_AdjList, d1_AdjList, d0_Own, d1_d2_Own, rDeltaXYZ_[cellI], fluidProps_.gamma).transpose();
    }
    else
    {
        coefsLimited[count++] = slope2_Limiter(deltaList, d0_AdjList, d1_AdjList, d0_Own, d1_d2_Own, rDeltaXYZ_[cellI]).transpose();
    }
}
/* Update d2 coefficients of polynomials */
count = 0;
forAll(mesh_.cells(), cellI)
{
    if (isTrouble[cellI])
    {
        coefs_[cellI] = coefsLimited[count++];
    }
}
