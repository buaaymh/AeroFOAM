std::vector<Mat6X5> d2CoefsLimited(parameter_.nTroubled);
parameter_.nTroubled = 0;
forAll(mesh_.cells(), cellI)
{
    if (!isTrouble_[cellI]) continue;
    const labelList& cellFaces = mesh_.cells()[cellI];
    std::vector<vector> normalList;
    std::vector<Mat5X6> d2CoefsList;
    normalList.reserve(mesh_.cellCells()[cellI].size());
    d2CoefsList.reserve(mesh_.cellCells()[cellI].size()+1);
    d2CoefsList.emplace_back(coefs_[cellI].block<6,5>(3,0).transpose());
    forAll(cellFaces, f)
    {
        const label faceI = cellFaces[f];
        if (mesh_.isInternalFace(faceI))
        {
            label j = owner[faceI];
            if (j == cellI) { j = neighbour[faceI]; }
            normalList.emplace_back(normal[faceI]);
            d2CoefsList.emplace_back(transformD2Coefs(coefs_[j].block<6,5>(3,0).transpose(), rDeltaXYZ_[cellI], rDeltaXYZ_[j]));
        }
    }
    if (parameter_.characteristic)
    {
        d2CoefsLimited[parameter_.nTroubled++] = Eigen_WBAP_L2_Limiter<6>(normalList, d2CoefsList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], parameter_.np).transpose();
    }
    else
    {
        d2CoefsLimited[parameter_.nTroubled++] = WBAP_L2_Limiter<6>(d2CoefsList, parameter_.np).transpose();
    }
}
/* Update d2 coefficients of polynomials */
parameter_.nTroubled = 0;
forAll(mesh_.cells(), cellI)
{
    if (isTrouble_[cellI])
    {
        coefs_[cellI].block<6,5>(3,0) = d2CoefsLimited[parameter_.nTroubled++];
    }
}
/* Limiting d1 coefficients of polynomials */
std::vector<Mat3X5> d1CoefsLimited(parameter_.nTroubled);
parameter_.nTroubled = 0;
forAll(mesh_.cells(), cellI)
{
    if (!isTrouble_[cellI]) continue;
    const labelList& cellFaces = mesh_.cells()[cellI];
    std::vector<vector> normalList;
    std::vector<Mat5X3> d1CoefsList;
    normalList.reserve(mesh_.cellCells()[cellI].size());
    d1CoefsList.reserve(mesh_.cellCells()[cellI].size()+1);
    d1CoefsList.emplace_back(coefs_[cellI].block<3,5>(0,0).transpose());
    forAll(cellFaces, f)
    {
        const label faceI = cellFaces[f];
        if (mesh_.isInternalFace(faceI))
        {
            normalList.emplace_back(normal[faceI]);
            label j = owner[faceI];
            vector d_ij = delta_[faceI];
            if (j == cellI)
            {
                j = neighbour[faceI];
                d_ij *= -1;
            }
            d1CoefsList.emplace_back(transformD1Coefs(coefs_[j].transpose(), d_ij, rDeltaXYZ_[cellI], rDeltaXYZ_[j]));
        }
    }
    if (parameter_.characteristic)
    {
        d1CoefsLimited[parameter_.nTroubled++] = Eigen_WBAP_L2_Limiter<3>(normalList, d1CoefsList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], parameter_.np).transpose();
    }
    else
    {
        d1CoefsLimited[parameter_.nTroubled++] = WBAP_L2_Limiter<3>(d1CoefsList, parameter_.np).transpose();
    }
}
/* Update d1 coefficients of polynomials */
parameter_.nTroubled = 0;
forAll(mesh_.cells(), cellI)
{
    if (isTrouble_[cellI])
    {
        coefs_[cellI].block<3,5>(0,0) = d1CoefsLimited[parameter_.nTroubled++];
    }
}