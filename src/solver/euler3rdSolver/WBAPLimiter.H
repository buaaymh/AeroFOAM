const scalar np = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("np", 10);
std::vector<Mat6X5> d2CoefsLimited(count);
count = 0;
forAll(mesh_.cells(), cellI)
{
    if (!isTrouble[cellI]) continue;
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
            normalList.emplace_back(normal_[faceI]);
            d2CoefsList.emplace_back(transformD2Coefs(coefs_[j].block<6,5>(3,0), rDeltaXYZ_[cellI], rDeltaXYZ_[j]).transpose());
        }
    }
    if (characteristic)
    {
        d2CoefsLimited[count++] = Eigen_WBAP_L2_Limiter<6>(normalList, d2CoefsList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], fluidProps_.gamma, np).transpose();
    }
    else
    {
        d2CoefsLimited[count++] = WBAP_L2_Limiter<6>(d2CoefsList, np).transpose();
    }
}
/* Update d2 coefficients of polynomials */
count = 0;
forAll(mesh_.cells(), cellI)
{
    if (isTrouble[cellI])
    {
        coefs_[cellI].block<6,5>(3,0) = d2CoefsLimited[count++];
    }
}
/* Limiting d1 coefficients of polynomials */
std::vector<Mat3X5> d1CoefsLimited(count);
count = 0;
forAll(mesh_.cells(), cellI)
{
    if (!isTrouble[cellI]) continue;
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
            normalList.emplace_back(normal_[faceI]);
            label j = owner[faceI];
            vector d_ij = delta_[faceI];
            if (j == cellI)
            {
                j = neighbour[faceI];
                d_ij *= -1;
            }
            d1CoefsList.emplace_back(transformD1Coefs(coefs_[j], d_ij, rDeltaXYZ_[cellI], rDeltaXYZ_[j]).transpose());
        }
    }
    if (characteristic)
    {
        d1CoefsLimited[count++] = Eigen_WBAP_L2_Limiter<3>(normalList, d1CoefsList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], fluidProps_.gamma, np).transpose();
    }
    else
    {
        d1CoefsLimited[count++] = WBAP_L2_Limiter<3>(d1CoefsList, np).transpose();
    }
}
/* Update d1 coefficients of polynomials */
count = 0;
forAll(mesh_.cells(), cellI)
{
    if (isTrouble[cellI])
    {
        coefs_[cellI].block<3,5>(0,0) = d1CoefsLimited[count++];
    }
}