isTrouble_ = std::vector<bool>(mesh_.nCells(), false);
if (parameter_.indicator)
{
    scalarField rho_error(mesh_.nCells(), 0.0);
    scalarField rho_max(rho_.primitiveField());
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        Col6X1 d2RhoCoefs_i  = coefs_[i].block<6,1>(3,0);
        Col6X1 d2RhoCoefs_j  = coefs_[j].block<6,1>(3,0);
        vector rDeltaXYZ_ij = cmptDivide(rDeltaXYZ_[j], rDeltaXYZ_[i]);
        Col6X1 transfoam_ij = Col6X1(sqr(rDeltaXYZ_ij.x()), rDeltaXYZ_ij.x()*rDeltaXYZ_ij.y(), rDeltaXYZ_ij.x()*rDeltaXYZ_ij.z(),
                                     sqr(rDeltaXYZ_ij.y()), rDeltaXYZ_ij.y()*rDeltaXYZ_ij.z(),
                                     sqr(rDeltaXYZ_ij.z()));
        Col6X1 PolyCenter_i = Col6X1(0.5*basisMean_[i][0], basisMean_[i][1], basisMean_[i][2],
                                     0.5*basisMean_[i][3], basisMean_[i][4],
                                     0.5*basisMean_[i][5]);
        Col6X1 PolyCenter_j = Col6X1(0.5*basisMean_[j][0], basisMean_[j][1], basisMean_[j][2],
                                     0.5*basisMean_[j][3], basisMean_[j][4],
                                     0.5*basisMean_[j][5]);
        rho_error[i] += mag((d2RhoCoefs_i-(d2RhoCoefs_j.array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        rho_error[j] += mag((d2RhoCoefs_j-(d2RhoCoefs_i.array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = rho_error[cellI]/(mesh_.cellCells()[cellI].size()*N_h_[cellI]*rho_max[cellI]);
        if (IS_[cellI] > parameter_.IS)
        {
            isTrouble_[cellI] = true;
            parameter_.nTroubled++;
        }
    }
}
else
{
    scalarField error(mesh_.nCells(), 0.0);
    scalarField rho_max(rho_.primitiveField());
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        scalar error_ij = rho_[i] + coefs_[i].col(0).dot(polynomialDn0(mesh_.Cf()[faceI]-mesh_.C()[i], rDeltaXYZ_[i], basisMean_[i]))
                        - rho_[j] - coefs_[j].col(0).dot(polynomialDn0(mesh_.Cf()[faceI]-mesh_.C()[j], rDeltaXYZ_[j], basisMean_[j]));
        error_ij = mag(error_ij);
        error[i] += error_ij;
        error[j] += error_ij;
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = error[cellI]/(mesh_.cellCells()[cellI].size()*N_h_[cellI]*rho_max[cellI]);
        if (IS_[cellI] > parameter_.IS)
        {
            isTrouble_[cellI] = true;
            parameter_.nTroubled++;
        }
    }
}