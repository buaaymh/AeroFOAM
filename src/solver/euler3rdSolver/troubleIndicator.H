const Switch indicator = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("indicator", false);
isTrouble_ = std::vector<bool>(mesh_.nCells(), false);
if (indicator)
{
    scalarField rho_error(mesh_.nCells(), 0.0);
    scalarField rhoE_error(mesh_.nCells(), 0.0);
    scalarField rho_max(rho_.primitiveField());
    scalarField rhoE_max(rhoE_.primitiveField());
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        rhoE_max[i] = max(rhoE_max[i], rhoE_[j]);
        rhoE_max[j] = max(rhoE_max[j], rhoE_[i]);
        if (isP1Cell_[i] && isP1Cell_[j]) continue;
        Col9X1 u_i = polynomialDn0(mesh_.Cf()[faceI]-mesh_.C()[i], rDeltaXYZ_[i], basisMean_[i]);
        Col9X1 u_j = polynomialDn0(mesh_.Cf()[faceI]-mesh_.C()[j], rDeltaXYZ_[j], basisMean_[j]);
        scalar error_ij = mag(rho_[i] - rho_[j] + coefs_[i].col(0).dot(u_i) - coefs_[j].col(0).dot(u_j));
        rho_error[i] += error_ij;
        rho_error[j] += error_ij;
        error_ij = mag(rhoE_[i] - rhoE_[j] + coefs_[i].col(4).dot(u_i) - coefs_[j].col(4).dot(u_j));
        rhoE_error[i] += error_ij;
        rhoE_error[j] += error_ij;
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = max(rho_error[cellI]/(N_h_[cellI]*rho_max[cellI]),
                         rhoE_error[cellI]/(N_h_[cellI]*rhoE_max[cellI]));
        if (IS_[cellI] > IS)
        {
            isTrouble_[cellI] = true;
            count++;
        }
    }
}
else
{
    scalarField rho_error(mesh_.nCells(), 0.0);
    scalarField rhoE_error(mesh_.nCells(), 0.0);
    scalarField rho_max(rho_.primitiveField());
    scalarField rhoE_max(rhoE_.primitiveField());
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        rhoE_max[i] = max(rhoE_max[i], rhoE_[j]);
        rhoE_max[j] = max(rhoE_max[j], rhoE_[i]);
        if (isP1Cell_[i] && isP1Cell_[j]) continue;
        Col6X1 d2RhoCoefs_i  = coefs_[i].block<6,1>(3,0);
        Col6X1 d2RhoECoefs_i = coefs_[i].block<6,1>(3,4);
        Col6X1 d2RhoCoefs_j  = coefs_[j].block<6,1>(3,0);
        Col6X1 d2RhoECoefs_j = coefs_[j].block<6,1>(3,4);
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
        rho_error[i]  += mag((d2RhoCoefs_i-(d2RhoCoefs_j.array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        rho_error[j]  += mag((d2RhoCoefs_j-(d2RhoCoefs_i.array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
        rhoE_error[i] += mag((d2RhoECoefs_i-(d2RhoECoefs_j.array()*transfoam_ij.array()).matrix()).dot(PolyCenter_i));
        rhoE_error[j] += mag((d2RhoECoefs_j-(d2RhoECoefs_i.array()/transfoam_ij.array()).matrix()).dot(PolyCenter_j));
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = max(rho_error[cellI]/(N_h_[cellI]*rho_max[cellI]),
                         rhoE_error[cellI]/(N_h_[cellI]*rhoE_max[cellI]));
        if (IS_[cellI] > IS)
        {
            isTrouble_[cellI] = true;
            count++;
        }
    }
}