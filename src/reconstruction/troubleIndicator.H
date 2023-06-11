const Switch indicator = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<Switch>("indicator", false);
isTrouble_ = std::vector<bool>(mesh_.nCells(), false);
if (indicator)
{
    scalarField rho_error(mesh_.nCells(), 0.0);
    scalarField rho_max(rho_.primitiveField());
    scalarField h_max(mesh_.nCells(), 0.0);
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        h_max[i] = max(h_max[i], d_ij_[faceI]);
        h_max[j] = max(h_max[j], d_ij_[faceI]);
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        if (isP1Cell_[i] && isP1Cell_[j]) continue;
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
        IS_[cellI] = rho_error[cellI]/(mesh_.cellCells()[cellI].size()*pow(h_max[cellI], 1.5)*rho_max[cellI]);
        if (IS_[cellI] > IS)
        {
            isTrouble_[cellI] = true;
            count++;
        }
    }
}
else
{
    scalarField rhoU = magSqr(rhoU_.primitiveField());
    scalarField rho_error(mesh_.nCells(), 0.0);
    scalarField rhoU_error(mesh_.nCells(), 0.0);
    scalarField rhoE_error(mesh_.nCells(), 0.0);
    scalarField rho_max(rho_.primitiveField());
    scalarField rhoU_max(rhoU);
    scalarField rhoE_max(rhoE_.primitiveField());
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const scalar h_k = 1.0/pow(d_ij_[faceI], 1.5);
        rho_max[i]  = max(rho_max[i],  rho_[j]);
        rho_max[j]  = max(rho_max[j],  rho_[i]);
        rhoU_max[i] = max(rhoU_max[i], rhoU[j]);
        rhoU_max[j] = max(rhoU_max[j], rhoU[i]);
        rhoE_max[i] = max(rhoE_max[i], rhoE_[j]);
        rhoE_max[j] = max(rhoE_max[j], rhoE_[i]);
        if (isP1Cell_[i] && isP1Cell_[j]) continue;
        Mat6X5 d2Coefs_i  = coefs_[i].block<6,5>(3,0);
        Mat6X5 d2Coefs_j  = coefs_[j].block<6,5>(3,0);
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
        rho_error[i]  += h_k*mag((d2Coefs_i.col(0)-(d2Coefs_j.col(0).array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        rhoU_error[i] += h_k*sqr((d2Coefs_i.col(1)-(d2Coefs_j.col(1).array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        rhoU_error[i] += h_k*sqr((d2Coefs_i.col(2)-(d2Coefs_j.col(2).array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        rhoU_error[i] += h_k*sqr((d2Coefs_i.col(3)-(d2Coefs_j.col(3).array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        rhoE_error[i] += h_k*mag((d2Coefs_i.col(4)-(d2Coefs_j.col(4).array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        
        rho_error[j]  += h_k*mag((d2Coefs_j.col(0)-(d2Coefs_i.col(0).array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
        rhoU_error[j] += h_k*sqr((d2Coefs_j.col(1)-(d2Coefs_i.col(1).array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
        rhoU_error[j] += h_k*sqr((d2Coefs_j.col(2)-(d2Coefs_i.col(2).array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
        rhoU_error[j] += h_k*sqr((d2Coefs_j.col(3)-(d2Coefs_i.col(3).array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
        rhoE_error[j] += h_k*mag((d2Coefs_j.col(4)-(d2Coefs_i.col(4).array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = max(rho_error[cellI]/(mesh_.cellCells()[cellI].size()*rho_max[cellI]),
                         rhoE_error[cellI]/(mesh_.cellCells()[cellI].size()*rhoE_max[cellI]));
        IS_[cellI] = max(IS_[cellI], rhoU_error[cellI]/(mesh_.cellCells()[cellI].size()*(rhoU_max[cellI]+1e-6)));
        if (IS_[cellI] > IS)
        {
            isTrouble_[cellI] = true;
            count++;
        }
    }
}