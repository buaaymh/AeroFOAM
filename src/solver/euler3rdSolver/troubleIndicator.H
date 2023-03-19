const Switch indicator = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("indicator", false);
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
        Col9X1 u_i = polynomialDn0(mesh_.Cf()[faceI]-mesh_.C()[i], rDeltaXYZ_[i], basisMean_[i]);
        Col9X1 u_j = polynomialDn0(mesh_.Cf()[faceI]-mesh_.C()[j], rDeltaXYZ_[j], basisMean_[j]);
        scalar error_ij = mag(rho_[i] - rho_[j] + coefs_[i].col(0).dot(u_i) - coefs_[j].col(0).dot(u_j));
        rho_error[i] += error_ij;
        rho_error[j] += error_ij;
        error_ij = mag(rhoE_[i] - rhoE_[j] + coefs_[i].col(4).dot(u_i) - coefs_[j].col(4).dot(u_j));
        rhoE_error[i] += error_ij;
        rhoE_error[j] += error_ij;
    }
    forAll(mesh_.boundary(), patchI)
    {
        const word type = mesh_.boundary()[patchI].type();
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const fvPatchScalarField& rhoEBound = rhoE_.boundaryField()[patchI];
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const vectorField& Cf = mesh_.boundary()[patchI].Cf();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                Col9X1 u_i = polynomialDn0(Cf[faceI]-mesh_.C()[i], rDeltaXYZ_[i], basisMean_[i]);
                rhoE_error[i] += mag(rhoE_[i] - rhoEBound[faceI] + coefs_[i].col(4).dot(u_i));
            }
        }
    }
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = max(rho_error[cellI]/(N_h_[cellI]*rho_max[cellI]),
                         rhoE_error[cellI]/(N_h_[cellI]*rhoE_max[cellI]));
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
        error[j] += mag(rho_[i] + coefs_[i].col(0).dot(polynomialDn0( delta_[faceI], rDeltaXYZ_[i], basisMean_[i])) - rho_mid[j]);
        error[i] += mag(rho_[j] + coefs_[j].col(0).dot(polynomialDn0(-delta_[faceI], rDeltaXYZ_[j], basisMean_[j])) - rho_mid[i]);
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