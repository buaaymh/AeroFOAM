void Foam::euler3rdSolver::updateCoefficients()
{
    const Switch adaptive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("adaptive", false);
    if (adaptive)
    {
        label count = 0;
        p0_ = std::vector<bool>(mesh_.nCells(), false);
        forAll(mesh_.C(), cellI)
        {
            if (cmptSum(cmptMag(rhoGrad_[cellI])) < 1e-6)
            {
                p0_[cellI] = true;
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))
                               /scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << "P0 cell [%] = " << percent << ", ";
    }
    const vector vrWeight = mesh_.schemesDict().subDict("vrSchemes").lookup<vector>("weightList");
    const scalarField delta = mag(mesh_.delta());
    std::vector<Column> bRho = std::vector<Column>(mesh_.nCells(), Column::Zero());
    std::vector<Column> bUx  = std::vector<Column>(mesh_.nCells(), Column::Zero());
    std::vector<Column> bUy  = std::vector<Column>(mesh_.nCells(), Column::Zero());
    std::vector<Column> bUz  = std::vector<Column>(mesh_.nCells(), Column::Zero());
    std::vector<Column> bT   = std::vector<Column>(mesh_.nCells(), Column::Zero());

    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        if (p0_[i] && p0_[j]) continue;
        Column bRho_i = Column::Zero();
        Column bUx_i = Column::Zero();
        Column bUy_i = Column::Zero();
        Column bUz_i = Column::Zero();
        Column bT_i = Column::Zero();
        Column bRho_j = Column::Zero();
        Column bUx_j = Column::Zero();
        Column bUy_j = Column::Zero();
        Column bUz_j = Column::Zero();
        Column bT_j = Column::Zero();
        const scalar dRho0 = rho_[j] - rho_[i];
        const vector dU0   = U_[j] - U_[i];
        const scalar dT0   = T_[j] - T_[i];
        scalar temp = 0.25*mesh_.magSf()[faceI]*sqr(vrWeight[0])/delta[faceI];
        for (label gaussI = 0; gaussI != 4; ++gaussI)
        {
            const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
            const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
            const scalar wTemp = Quad4::w[gaussI] * temp;
            const scalar dRho = (dRho0 + (rhoGrad_[j]&delta_j) - (rhoGrad_[i]&delta_i))*wTemp;
            const vector dU   = (dU0 + (UGrad_[j]&delta_j)     - (UGrad_[i]&delta_i))  *wTemp;
            const scalar dT   = (dT0 + (TGrad_[j]&delta_j)     - (TGrad_[i]&delta_i))  *wTemp;
            //cell i
            Column Dn0 = Foam::basisPolyEigen(delta_i, rLengthScale_[i], basisConst_[i]);
            bRho_i += Dn0 * dRho;
            bUx_i  += Dn0 * dU.x();
            bUy_i  += Dn0 * dU.y();
            bUx_i  += Dn0 * dU.x();
            bT_i   += Dn0 * dT;
            //cell j
            Dn0 = Foam::basisPolyEigen(delta_j, rLengthScale_[j], basisConst_[j]);
            bRho_j -= Dn0 * dRho;
            bUx_j  -= Dn0 * dU.x();
            bUy_j  -= Dn0 * dU.y();
            bUx_j  -= Dn0 * dU.x();
            bT_j   -= Dn0 * dT;
        }
        bRho[i] += bRho_i;
        bUx[i]  += bUx_i;
        bUy[i]  += bUy_i;
        bUz[i]  += bUz_i;
        bT[i]   += bT_i;
        bRho[j] += bRho_j;
        bUx[j]  += bUx_j;
        bUy[j]  += bUy_j;
        bUz[j]  += bUz_j;
        bT[j]   += bT_j;
        
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        temp = delta[faceI]*mesh_.magSf()[faceI]*sqr(vrWeight[1]);
        const scalar dRho1 = (rhoGrad_[j]-rhoGrad_[i])&normal_[faceI] *temp;
        const vector dU1   = (UGrad_[j]-UGrad_[i])&normal_[faceI]     *temp;
        const scalar dT1   = (TGrad_[j]-TGrad_[i])&normal_[faceI]     *temp;
        //cell i
        Column Dn1 = basisPolyDn1Eigen(delta_i, rLengthScale_[i], normal_[faceI]);
        bRho[i] += Dn1 * dRho1;
        bUx[i]  += Dn1 * dU1.x();
        bUy[i]  += Dn1 * dU1.y();
        bUz[i]  += Dn1 * dU1.z();
        bT[i]   += Dn1 * dT1;
        //cell j
        Dn1 = basisPolyDn1Eigen(delta_j, rLengthScale_[j], normal_[faceI]);
        bRho[j] -= Dn1 * dRho1;
        bUx[j]  -= Dn1 * dU1.x();
        bUy[j]  -= Dn1 * dU1.y();
        bUz[j]  -= Dn1 * dU1.z();
        bT[j]   -= Dn1 * dT1;
    }
    // forAll(mesh_.boundary(), patchI)
    // {
    //     if (mesh_.boundary()[patchI].coupled()) continue;
    //     const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
    //     const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
    //     const scalarField  delta = mag(mesh_.boundary()[patchI].delta());
    //     const vectorField& normal = normal_.boundaryField()[patchI];
    //     const vectorField& Cf = mesh_.boundary()[patchI].Cf();
    //     const label start = mesh_.boundary()[patchI].start();
    //     const word name = mesh_.boundary()[patchI].name();
    //     const word type = mesh_.boundary()[patchI].type();
    //     const fvPatchScalarField& rhoB = rho_.boundaryField()[patchI];
    //     const fvPatchVectorField& UB = U_.boundaryField()[patchI];
    //     const fvPatchScalarField& TB = T_.boundaryField()[patchI];
    //     std::vector<vector> quadToints(4, vector::zero);
    //     forAll(bfaceCells, faceI)
    //     {
    //         const label i = bfaceCells[faceI];
    //         gaussQuad4(start+faceI, quadToints);
    //         Column bRho_i = Column::Zero();
    //         Column bUx_i = Column::Zero();
    //         Column bUy_i = Column::Zero();
    //         Column bUz_i = Column::Zero();
    //         Column bT_i = Column::Zero();
    //         const scalar dRho0 = rhoB[faceI] - rho_[i];
    //         const vector dU0   = UB[faceI]   - U_[i];
    //         const scalar dT0   = TB[faceI]   - T_[i];
    //         scalar temp = 0.25*magSf[faceI]*sqr(vrWeight[0])/delta[faceI];
    //         if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
    //         {
    //             for (label gaussI = 0; gaussI != 4; ++gaussI)
    //             {
    //                 const vector delta_i = quadToints[gaussI] - mesh_.C()[i];
    //                 const scalar wTemp = Quad4::w[gaussI] * temp;
    //                 const scalar dRho = (dRho0 - (rhoGrad_[i]&delta_i))*wTemp;
    //                 const vector dU   = (dU0   - (UGrad_[i]&delta_i))  *wTemp;
    //                 const scalar dT   = (dT0   - (TGrad_[i]&delta_i))  *wTemp;
    //                 //cell i
    //                 Column Dn0 = Foam::basisPolyEigen(delta_i, rLengthScale_[i], basisConst_[i]);
    //                 bRho_i += Dn0 * dRho;
    //                 bUx_i  += Dn0 * dU.x();
    //                 bUy_i  += Dn0 * dU.y();
    //                 bUx_i  += Dn0 * dU.x();
    //                 bT_i   += Dn0 * dT;
    //             }
    //             const vector delta_i = Cf[faceI] - mesh_.C()[i];
    //             temp = delta[faceI]*magSf[faceI]*sqr(vrWeight[1]);
    //             const scalar dRho1 = -(rhoGrad_[i]&normal[faceI])*temp;
    //             const vector dU1   = -(UGrad_[i]&normal[faceI])  *temp;
    //             const scalar dT1   = -(TGrad_[i]&normal[faceI])  *temp;
    //             //cell i
    //             Column Dn1 = basisPolyDn1Eigen(delta_i, rLengthScale_[i], normal[faceI]);
    //             bRho_i += Dn1 * dRho1;
    //             bUx_i  += Dn1 * dU1.x();
    //             bUy_i  += Dn1 * dU1.y();
    //             bUz_i  += Dn1 * dU1.z();
    //             bT_i   += Dn1 * dT1;
    //         }
    //         if (name == "inlet" || name == "outlet" || name == "farField")
    //         {
    //             for (label gaussI = 0; gaussI != 4; ++gaussI)
    //             {
    //                 const vector delta_i = quadToints[gaussI] - mesh_.C()[i];
    //                 const scalar wTemp = Quad4::w[gaussI] * temp;
    //                 const scalar dRho = (dRho0 - (rhoGrad_[i]&delta_i))*wTemp;
    //                 const vector dU   = (dU0   - (UGrad_[i]&delta_i))  *wTemp;
    //                 const scalar dT   = (dT0   - (TGrad_[i]&delta_i))  *wTemp;
    //                 //cell i
    //                 Column Dn0 = Foam::basisPolyEigen(delta_i, rLengthScale_[i], basisConst_[i]);
    //                 bRho_i += Dn0 * dRho;
    //                 bUx_i  += Dn0 * dU.x();
    //                 bUy_i  += Dn0 * dU.y();
    //                 bUx_i  += Dn0 * dU.x();
    //                 bT_i   += Dn0 * dT;
    //             }
    //         }
    //         bRho[i] += bRho_i;
    //         bUx[i]  += bUx_i;
    //         bUy[i]  += bUy_i;
    //         bUz[i]  += bUz_i;
    //         bT[i]   += bT_i;
    //     }
    // }
    // forAll(mesh_.boundary(), patchI)
    // {
    //     if (mesh_.boundary()[patchI].coupled()) continue;
    //     const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
    //     const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
    //     const scalarField  delta = mag(mesh_.boundary()[patchI].delta());
    //     const label start = mesh_.boundary()[patchI].start();
    //     const fvPatchScalarField& rhoB = rho_.boundaryField()[patchI];
    //     const fvPatchVectorField& UB = U_.boundaryField()[patchI];
    //     const fvPatchScalarField& TB = T_.boundaryField()[patchI];
    //     std::vector<vector> quadToints(4, vector::zero);
    //     forAll(bfaceCells, faceI)
    //     {
    //         const label i = bfaceCells[faceI];
    //         gaussQuad4(start+faceI, quadToints);
    //         Column bRho_i = Column::Zero();
    //         Column bUx_i = Column::Zero();
    //         Column bUy_i = Column::Zero();
    //         Column bUz_i = Column::Zero();
    //         Column bT_i = Column::Zero();
    //         const scalar dRho0 = rhoB[faceI] - rho_[i];
    //         const vector dU0   = UB[faceI]   - U_[i];
    //         const scalar dT0   = TB[faceI]   - T_[i];
    //         scalar temp = 0.25*magSf[faceI]*vrWeight[0]/delta[faceI];
    //         for (label gaussI = 0; gaussI != 4; ++gaussI)
    //         {
    //             const vector delta_i = quadToints[gaussI] - mesh_.C()[i];
    //             const scalar wTemp = Quad4::w[gaussI] * temp;
    //             const scalar dRho = (dRho0 - (rhoGrad_[i]&delta_i))*wTemp;
    //             const vector dU   = (dU0   - (UGrad_[i]&delta_i))  *wTemp;
    //             const scalar dT   = (dT0   - (TGrad_[i]&delta_i))  *wTemp;
    //             //cell i
    //             Column Dn0 = Foam::basisPolyEigen(delta_i, rLengthScale_[i], basisConst_[i]);
    //             bRho_i += Dn0 * dRho;
    //             bUx_i  += Dn0 * dU.x();
    //             bUy_i  += Dn0 * dU.y();
    //             bUx_i  += Dn0 * dU.x();
    //             bT_i   += Dn0 * dT;
    //         }
    //         bRho[i] += bRho_i;
    //         bUx[i]  += bUx_i;
    //         bUy[i]  += bUy_i;
    //         bUz[i]  += bUz_i;
    //         bT[i]   += bT_i;
    //     }
    // }
    forAll(mesh_.cells(), cellI)
    {
        if (p0_[cellI])
        {
            d2Rho_[cellI] = symmTensor::zero;
            d2Ux_[cellI]  = symmTensor::zero;
            d2Uy_[cellI]  = symmTensor::zero;
            d2Uz_[cellI]  = symmTensor::zero;
            d2T_[cellI]   = symmTensor::zero;
            continue;
        }
        const labelList& cellFaces = mesh_.cells()[cellI];
        Column d2RhoI = bRho[cellI];
        Column d2UxI  = bUx[cellI];
        Column d2UyI  = bUy[cellI];
        Column d2UzI  = bUz[cellI];
        Column d2TI   = bT[cellI];
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j != cellI) /* j < cellI */
                {
                    d2RhoI += B_[faceI] * Column(d2Rho_[j][0], d2Rho_[j][1], d2Rho_[j][2],
                                                 d2Rho_[j][3], d2Rho_[j][4], d2Rho_[j][5]);
                    d2UxI  += B_[faceI] * Column(d2Ux_[j][0], d2Ux_[j][1], d2Ux_[j][2],
                                                 d2Ux_[j][3], d2Ux_[j][4], d2Ux_[j][5]);
                    d2UyI  += B_[faceI] * Column(d2Uy_[j][0], d2Uy_[j][1], d2Uy_[j][2],
                                                 d2Uy_[j][3], d2Uy_[j][4], d2Uy_[j][5]);
                    d2UzI  += B_[faceI] * Column(d2Uz_[j][0], d2Uz_[j][1], d2Uz_[j][2],
                                                 d2Uz_[j][3], d2Uz_[j][4], d2Uz_[j][5]);
                    d2TI   += B_[faceI] * Column(d2T_[j][0], d2T_[j][1], d2T_[j][2],
                                                 d2T_[j][3], d2T_[j][4], d2T_[j][5]);
                }
                else /* j > cellI */
                {
                    j = neighbour[faceI];
                    d2RhoI += B_[faceI].transpose() * Column(d2Rho_[j][0], d2Rho_[j][1], d2Rho_[j][2],
                                                             d2Rho_[j][3], d2Rho_[j][4], d2Rho_[j][5]);
                    d2UxI  += B_[faceI].transpose() * Column(d2Ux_[j][0], d2Ux_[j][1], d2Ux_[j][2],
                                                             d2Ux_[j][3], d2Ux_[j][4], d2Ux_[j][5]);
                    d2UyI  += B_[faceI].transpose() * Column(d2Uy_[j][0], d2Uy_[j][1], d2Uy_[j][2],
                                                             d2Uy_[j][3], d2Uy_[j][4], d2Uy_[j][5]);
                    d2UzI  += B_[faceI].transpose() * Column(d2Uz_[j][0], d2Uz_[j][1], d2Uz_[j][2],
                                                             d2Uz_[j][3], d2Uz_[j][4], d2Uz_[j][5]);
                    d2TI   += B_[faceI].transpose() * Column(d2T_[j][0], d2T_[j][1], d2T_[j][2],
                                                             d2T_[j][3], d2T_[j][4], d2T_[j][5]);
                }
            }
        }
        d2RhoI = rA_[cellI] * d2RhoI;
        d2UxI  = rA_[cellI] * d2UxI;
        d2UyI  = rA_[cellI] * d2UyI;
        d2UzI  = rA_[cellI] * d2UzI;
        d2TI   = rA_[cellI] * d2TI;
        d2Rho_[cellI] = symmTensor(d2RhoI(0), d2RhoI(1), d2RhoI(2), d2RhoI(3), d2RhoI(4), d2RhoI(5));
        d2Ux_[cellI]  = symmTensor(d2UxI(0), d2UxI(1), d2UxI(2), d2UxI(3), d2UxI(4), d2UxI(5));
        d2Uy_[cellI]  = symmTensor(d2UyI(0), d2UyI(1), d2UyI(2), d2UyI(3), d2UyI(4), d2UyI(5));
        d2Uz_[cellI]  = symmTensor(d2UzI(0), d2UzI(1), d2UzI(2), d2UzI(3), d2UzI(4), d2UzI(5));
        d2T_[cellI]   = symmTensor(d2TI(0), d2TI(1), d2TI(2), d2TI(3), d2TI(4), d2TI(5));
    }
    d2Rho_.correctBoundaryConditions();
    d2Ux_.correctBoundaryConditions();
    d2Uy_.correctBoundaryConditions();
    d2Uz_.correctBoundaryConditions();
    d2T_.correctBoundaryConditions();
}