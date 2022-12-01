const Switch indicator = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("indicator", false);
if (indicator)
{
    scalarField error(mesh_.nCells(), 0.0);
    scalarField rhoU_mag(mag(rhoU_.primitiveField()));
    scalarField rhoU_max(rhoU_mag);
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rhoU_max[i] = max(rhoU_max[i], rhoU_mag[j]);
        rhoU_max[j] = max(rhoU_max[j], rhoU_mag[i]);
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        Col3X1 dRhoU_i = coefs_[i].block<9,3>(0,1).transpose()*polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
        Col3X1 dRhoU_j = coefs_[j].block<9,3>(0,1).transpose()*polynomialDn0(delta_j, rDeltaXYZ_[j], basisMean_[j]);
        scalar error_ij = mag((rhoU_[i]+vector(dRhoU_i(0), dRhoU_i(1), dRhoU_i(2)) - 
                               rhoU_[j]-vector(dRhoU_j(0), dRhoU_j(1), dRhoU_j(2)))&normal_[faceI]);
        error[j] += error_ij;
        error[i] += error_ij;
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = error[cellI]/(N_h_[cellI]*rhoU_max[cellI]+1e-6);
        if (IS_[cellI] > IS)
        {
            isTrouble[cellI] = true;
            count++;
        }
    }
}
else
{
    scalarField error(mesh_.nCells(), 0.0);
    scalarField rho_max(rho_.primitiveField());
    scalarField rho_mid(mesh_.nCells());
    forAll(mesh_.cells(), cellI)
    {
        rho_mid[cellI] = rho_[cellI] + coefs_[cellI].col(0).dot(polynomialDn0(vector::zero, rDeltaXYZ_[cellI], basisMean_[cellI]));
    }
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        error[j] += mag(rho_[i] + coefs_[i].col(0).dot(polynomialDn0( delta[faceI], rDeltaXYZ_[i], basisMean_[i])) - rho_mid[j]);
        error[i] += mag(rho_[j] + coefs_[j].col(0).dot(polynomialDn0(-delta[faceI], rDeltaXYZ_[j], basisMean_[j])) - rho_mid[i]);
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = error[cellI]/(N_h_[cellI]*rho_max[cellI]);
        if (IS_[cellI] > IS)
        {
            isTrouble[cellI] = true;
            count++;
        }
    }
}