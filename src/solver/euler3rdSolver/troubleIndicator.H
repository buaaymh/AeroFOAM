const Switch indicator = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("indicator", false);
if (indicator)
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
    forAll(mesh_.boundary(), patchI)
    {
        const word type = mesh_.boundary()[patchI].type();
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const vectorField& normal = normal_.boundaryField()[patchI];
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const vectorField& Cf = mesh_.boundary()[patchI].Cf();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                Col3X1 dRhoU_i = coefs_[i].block<9,3>(0,1).transpose()*polynomialDn0(Cf[faceI]-mesh_.C()[i], rDeltaXYZ_[i], basisMean_[i]);
                error[i] += mag((rhoU_[i]+vector(dRhoU_i(0),dRhoU_i(1),dRhoU_i(2)))&normal[faceI]);
            }
        }
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