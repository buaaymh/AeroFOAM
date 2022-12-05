forAll(mesh_.cells(), cellI)
{
    if (!isTrouble[cellI]) continue;
    const labelList& cellFaces = mesh_.cells()[cellI];
    std::vector<vector> normalList;
    std::vector<vector> deltaList;
    std::vector<Col5X1> d0_AdjList;
    normalList.reserve(mesh_.cellCells()[cellI].size());
    deltaList.reserve(mesh_.cellCells()[cellI].size());
    d0_AdjList.reserve(mesh_.cellCells()[cellI].size());
    const Col5X1 d0_Own(rho_[cellI], rhoU_[cellI][0], rhoU_[cellI][1], rhoU_[cellI][2], rhoE_[cellI]);
    const Mat5X3 d1_Own = coefs_[cellI].block<3,5>(0,0).transpose();
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
        }
    }
    if (characteristic)
    {
        coefs_[cellI].block<3,5>(0,0) = slope1_Limiter(normalList, deltaList, d0_AdjList, d0_Own, d1_Own, rDeltaXYZ_[cellI], fluidProps_.gamma).transpose();
    }
    else
    {
        coefs_[cellI].block<3,5>(0,0) = slope1_Limiter(deltaList, d0_AdjList, d0_Own, d1_Own, rDeltaXYZ_[cellI]).transpose();
    }
    coefs_[cellI].block<6,5>(3,0) = Mat6X5::Zero();
}