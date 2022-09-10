void Foam::vrScheme::updateCoefficients
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    volVectorField& rhoGrad,
    volTensorField& UGrad,
    volVectorField& TGrad,
    volSymmTensorField& d2Rho,
    volSymmTensorField& d2Ux,
    volSymmTensorField& d2Uy,
    volSymmTensorField& d2Uz,
    volSymmTensorField& d2T,
    const scalar gamma
)
{
    if (adaptive_) { detectP0(rhoGrad); }
    bRho_ = std::vector<vrScheme::Column>(mesh_.nCells(), vrScheme::Column::Zero());
    bUx_ = std::vector<vrScheme::Column>(mesh_.nCells(), vrScheme::Column::Zero());
    bUy_ = std::vector<vrScheme::Column>(mesh_.nCells(), vrScheme::Column::Zero());
    bUz_ = std::vector<vrScheme::Column>(mesh_.nCells(), vrScheme::Column::Zero());
    bT_ = std::vector<vrScheme::Column>(mesh_.nCells(), vrScheme::Column::Zero());
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        if (p0_[i] && p0_[j]) continue;
        const vector normal = mesh_.Sf()[faceI]/mesh_.magSf()[faceI];
        vrScheme::Column bRho_i = vrScheme::Column::Zero();
        vrScheme::Column bUx_i = vrScheme::Column::Zero();
        vrScheme::Column bUy_i = vrScheme::Column::Zero();
        vrScheme::Column bUz_i = vrScheme::Column::Zero();
        vrScheme::Column bT_i = vrScheme::Column::Zero();
        vrScheme::Column bRho_j = vrScheme::Column::Zero();
        vrScheme::Column bUx_j = vrScheme::Column::Zero();
        vrScheme::Column bUy_j = vrScheme::Column::Zero();
        vrScheme::Column bUz_j = vrScheme::Column::Zero();
        vrScheme::Column bT_j = vrScheme::Column::Zero();
        const scalar dRho0 = rho[j] - rho[i];
        const vector dU0   = U[j] - U[i];
        const scalar dT0   = T[j] - T[i];
        scalar temp = 0.25 * mesh_.magSf()[faceI] * vrWeight_[0] / delta_[faceI];
        for (label gaussI = 0; gaussI != 4; ++gaussI)
        {
            const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
            const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
            const scalar wTemp = Quad4::w[gaussI] * temp;
            const scalar dRho = (dRho0 + (rhoGrad[j]&delta_j) - (rhoGrad[i]&delta_i)) * wTemp;
            const vector dU   = (dU0 + (UGrad[j]&delta_j) - (UGrad[i]&delta_i)) * wTemp;
            const scalar dT   = (dT0 + (TGrad[j]&delta_j) - (TGrad[i]&delta_i)) * wTemp;
            //cell i
            vrScheme::Column Dn0 = Foam::basisPolyEigen(delta_i, rLengthScale_[i], basisConst_[i]);
            bRho_i += dRho   * Dn0;
            bUx_i  += dU.x() * Dn0;
            bUy_i  += dU.y() * Dn0;
            bUx_i  += dU.x() * Dn0;
            bT_i   += dT     * Dn0;
            //cell j
            Dn0 = Foam::basisPolyEigen(delta_j, rLengthScale_[j], basisConst_[j]);
            bRho_j -= dRho   * Dn0;
            bUx_j  -= dU.x() * Dn0;
            bUy_j  -= dU.y() * Dn0;
            bUx_j  -= dU.x() * Dn0;
            bT_j   -= dT     * Dn0;
        }
        bRho_[i] += bRho_i;
        bUx_[i] += bUx_i;
        bUy_[i] += bUy_i;
        bUz_[i] += bUz_i;
        bT_[i] += bT_i;
        bRho_[j] += bRho_j;
        bUx_[j] += bUx_j;
        bUy_[j] += bUy_j;
        bUz_[j] += bUz_j;
        bT_[j] += bT_j;
        
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        temp = delta_[faceI] * mesh_.magSf()[faceI] * vrWeight_[1];
        const scalar dRho1 = (rhoGrad[j]-rhoGrad[i])&normal * temp;
        const vector dU1   = (UGrad[j]-UGrad[i])&normal * temp;
        const scalar dT1   = (TGrad[j]-TGrad[i])&normal * temp;
        //cell i
        vrScheme::Column Dn1 = basisPolyDn1Eigen(delta_i, rLengthScale_[i], normal);
        bRho_[i] += dRho1 * Dn1;
        bUx_[i] += dU1.x() * Dn1;
        bUy_[i] += dU1.y() * Dn1;
        bUz_[i] += dU1.z() * Dn1;
        bT_[i] += dT1 * Dn1;
        //cell j
        Dn1 = basisPolyDn1Eigen(delta_j, rLengthScale_[j], normal);
        bRho_[j] -= dRho1 * Dn1;
        bUx_[j] -= dU1.x() * Dn1;
        bUy_[j] -= dU1.y() * Dn1;
        bUz_[j] -= dU1.z() * Dn1;
        bT_[j] -= dT1 * Dn1;
    }
    forAll(mesh_.cells(), cellI)
    {
        if (p0_[cellI])
        {
            d2Rho[cellI] = symmTensor::zero;
            d2Ux[cellI] = symmTensor::zero;
            d2Uy[cellI] = symmTensor::zero;
            d2Uz[cellI] = symmTensor::zero;
            d2T[cellI] = symmTensor::zero;
            continue;
        }
        const labelList& cellFaces = mesh_.cells()[cellI];
        vrScheme::Column d2RhoI = bRho_[cellI];
        vrScheme::Column d2UxI = bUx_[cellI];
        vrScheme::Column d2UyI = bUy_[cellI];
        vrScheme::Column d2UzI = bUz_[cellI];
        vrScheme::Column d2TI = bT_[cellI];
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j != cellI) /* j < cellI */
                {
                    d2RhoI += B_[faceI] * vrScheme::Column(d2Rho[j][0], d2Rho[j][1], d2Rho[j][2],
                                                            d2Rho[j][3], d2Rho[j][4], d2Rho[j][5]);
                    d2UxI  += B_[faceI] * vrScheme::Column(d2Ux[j][0], d2Ux[j][1], d2Ux[j][2],
                                                            d2Ux[j][3], d2Ux[j][4], d2Ux[j][5]);
                    d2UyI  += B_[faceI] * vrScheme::Column(d2Uy[j][0], d2Uy[j][1], d2Uy[j][2],
                                                            d2Uy[j][3], d2Uy[j][4], d2Uy[j][5]);
                    d2UzI  += B_[faceI] * vrScheme::Column(d2Uz[j][0], d2Uz[j][1], d2Uz[j][2],
                                                            d2Uz[j][3], d2Uz[j][4], d2Uz[j][5]);
                    d2TI += B_[faceI] * vrScheme::Column(d2T[j][0], d2T[j][1], d2T[j][2],
                                                            d2T[j][3], d2T[j][4], d2T[j][5]);
                }
                else /* j > cellI */
                {
                    j = neighbour[faceI];
                    d2RhoI += B_[faceI].transpose() * vrScheme::Column(d2Rho[j][0], d2Rho[j][1], d2Rho[j][2],
                                                                    d2Rho[j][3], d2Rho[j][4], d2Rho[j][5]);
                    d2UxI  += B_[faceI].transpose() * vrScheme::Column(d2Ux[j][0], d2Ux[j][1], d2Ux[j][2],
                                                                    d2Ux[j][3], d2Ux[j][4], d2Ux[j][5]);
                    d2UyI  += B_[faceI].transpose() * vrScheme::Column(d2Uy[j][0], d2Uy[j][1], d2Uy[j][2],
                                                                    d2Uy[j][3], d2Uy[j][4], d2Uy[j][5]);
                    d2UzI  += B_[faceI].transpose() * vrScheme::Column(d2Uz[j][0], d2Uz[j][1], d2Uz[j][2],
                                                                        d2Uz[j][3], d2Uz[j][4], d2Uz[j][5]);
                    d2TI += B_[faceI].transpose() * vrScheme::Column(d2T[j][0], d2T[j][1], d2T[j][2],
                                                                    d2T[j][3], d2T[j][4], d2T[j][5]);
                }
            }
        }
        d2RhoI = rA_[cellI] * d2RhoI;
        d2UxI = rA_[cellI] * d2UxI;
        d2UyI = rA_[cellI] * d2UyI;
        d2UzI = rA_[cellI] * d2UzI;
        d2TI = rA_[cellI] * d2TI;
        d2Rho[cellI] = symmTensor(d2RhoI(0), d2RhoI(1), d2RhoI(2), d2RhoI(3), d2RhoI(4), d2RhoI(5));
        d2Ux[cellI] = symmTensor(d2UxI(0), d2UxI(1), d2UxI(2), d2UxI(3), d2UxI(4), d2UxI(5));
        d2Uy[cellI] = symmTensor(d2UyI(0), d2UyI(1), d2UyI(2), d2UyI(3), d2UyI(4), d2UyI(5));
        d2Uz[cellI] = symmTensor(d2UzI(0), d2UzI(1), d2UzI(2), d2UzI(3), d2UzI(4), d2UzI(5));
        d2T[cellI] = symmTensor(d2TI(0), d2TI(1), d2TI(2), d2TI(3), d2TI(4), d2TI(5));
    }
    limitWBAP
    (
        rho, U, T,
        rhoGrad, UGrad, TGrad,
        d2Rho, d2Ux, d2Uy, d2Uz, d2T, gamma
    );
    if (positive_)
    {
        positiveCorrect
        (
            rho, U, T,
            rhoGrad, UGrad, TGrad,
            d2Rho, d2Ux, d2Uy, d2Uz, d2T
        );
    }
    rhoGrad.correctBoundaryConditions();
    UGrad.correctBoundaryConditions();
    TGrad.correctBoundaryConditions();
    d2Rho.correctBoundaryConditions();
    d2Ux.correctBoundaryConditions();
    d2Uy.correctBoundaryConditions();
    d2Uz.correctBoundaryConditions();
    d2T.correctBoundaryConditions();
}