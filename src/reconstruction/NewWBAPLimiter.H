const scalar np = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<scalar>("np", 2);
const auto& owner = mesh_.owner();
const auto& neighbour = mesh_.neighbour();
Mat5X5 L, R;
forAll(owner, faceI)
{
    const label i = owner[faceI];
    const label j = neighbour[faceI];
    if(!(isTrouble_[i] || isTrouble_[j])) continue;
    evaluateEigenMatrix(L, R, 0.5*(rho_[i]+rho_[j]), 0.5*(rhoU_[i]+rhoU_[j]), 0.5*(rhoE_[i]+rhoE_[j]), normal[faceI]);
    limitedCoefs_i_[faceI] = L*(coefs_[i].transpose());
    limitedCoefs_j_[faceI] = L*(coefs_[j].transpose());
    New_WBAP_Limiter(limitedCoefs_i_[faceI], limitedCoefs_j_[faceI], delta_[faceI], rDeltaXYZ_[i], rDeltaXYZ_[j], np);
    limitedCoefs_i_[faceI] = R*limitedCoefs_i_[faceI];
    limitedCoefs_j_[faceI] = R*limitedCoefs_j_[faceI];
}
forAll(mesh_.cells(), cellI)
{
    if(!isTrouble_[cellI]) continue;
    const labelList& cellFaces = mesh_.cells()[cellI];
    std::vector<Mat5X9> coefsList;
    coefsList.reserve(mesh_.cellCells()[cellI].size());
    forAll(cellFaces, f)
    {
        const label faceI = cellFaces[f];
        if (mesh_.isInternalFace(faceI))
        {
            if(owner[faceI] == cellI) coefsList.emplace_back(limitedCoefs_i_[faceI]);
            else coefsList.emplace_back(limitedCoefs_j_[faceI]);
        }
    }
    coefs_[cellI] = WBAP_L2_Limiter<9>(coefsList, 1).transpose();
}