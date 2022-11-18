void Foam::eulerPrimVar3rdSolver::reconstructionIter()
{
    /* Detect adaptive P0 cells */
    const Switch adaptive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("adaptive", false);
    if (adaptive)
    {
        label count = 0;
        isP0Cell_ = std::vector<bool>(mesh_.nCells(), false);
        forAll(mesh_.C(), cellI)
        {
            scalar indicator = cmptSum(cmptMag(rhoGrad_[cellI]))
                             + cmptSum(cmptMag(UGrad_[cellI]))
                             + cmptSum(cmptMag(pGrad_[cellI]));
            if ( indicator < 1e-6)
            {
                isP0Cell_[cellI] = true;
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << "P0 cell [%] = " << percent;
    }

    /* Calculate b of vr linear system */
    std::vector<Mat6X5> bW(mesh_.nCells(), Mat6X5::Zero());
    const vector vrWeight = mesh_.schemesDict().subDict("vrSchemes").lookup<vector>("weightList");
    const scalarField delta = mag(mesh_.delta());
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        if (isP0Cell_[i] && isP0Cell_[j]) continue;
        const scalar dRho0 = rho_[j] - rho_[i];
        const vector dU0   = U_[j]   - U_[i];
        const scalar dP0   = p_[j]   - p_[i];
        scalar temp = mesh_.magSf()[faceI]*sqr(vrWeight[0])/delta[faceI];
        for (label gaussI = 0; gaussI != 4; ++gaussI)
        {
            const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
            const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
            const scalar wTemp = Quad4::w[gaussI] * temp;
            const scalar dRho = wTemp*(dRho0 + (rhoGrad_[j]&delta_j) - (rhoGrad_[i]&delta_i));
            const vector dU   = wTemp*(dU0   + (UGrad_[j]&delta_j)   - (UGrad_[i]&delta_i));
            const scalar dP   = wTemp*(dP0   + (pGrad_[j]&delta_j)   - (pGrad_[i]&delta_i));
            Col5X1 dVars(dRho, dU[0], dU[1], dU[2], dP);
            //cell i
            bW[i] += basisPoly<Col6X1>(delta_i, rDeltaXYZ_[i], basisMean_[i])*dVars.transpose();
            //cell j
            bW[j] -= basisPoly<Col6X1>(delta_j, rDeltaXYZ_[j], basisMean_[j])*dVars.transpose();
        }
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        temp = delta[faceI]*mesh_.magSf()[faceI]*sqr(vrWeight[1]);
        const scalar dRho1 = temp*(normal_[faceI]&(rhoGrad_[j]-rhoGrad_[i]));
        const vector dU1   = temp*(normal_[faceI]&(UGrad_[j]  -UGrad_[i]));
        const scalar dP1   = temp*(normal_[faceI]&(pGrad_[j]  -pGrad_[i]));
        Col5X1 dGrad(dRho1, dU1[0], dU1[1], dU1[2], dP1);
        //cell i
        bW[i] += basisPolyDn1<Col6X1>(delta_i, rDeltaXYZ_[i], normal_[faceI])*dGrad.transpose();
        //cell j
        bW[j] -= basisPolyDn1<Col6X1>(delta_j, rDeltaXYZ_[j], normal_[faceI])*dGrad.transpose();
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField delta = mag(mesh_.boundary()[patchI].delta());
        const vectorField delta_i = mesh_.boundary()[patchI].Cf()-mesh_.boundary()[patchI].Cn();
        const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
        const word type = mesh_.boundary()[patchI].type();
        const label start = mesh_.boundary()[patchI].start();
        if (mesh_.boundary()[patchI].coupled())
        {
            const vectorField neighCn = ownerCn + mesh_.boundary()[patchI].delta();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Ux_neigh = d2Ux_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uy_neigh = d2Uy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uz_neigh = d2Uz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2P_neigh = d2P_.boundaryField()[patchI].patchNeighbourField();
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                const scalar dRho0 = rho_neigh[j] - rho_[i];
                const vector dU0   = U_neigh[j]   - U_[i];
                const scalar dP0   = p_neigh[j]   - p_[i];
                scalar temp = magSf[j]*sqr(vrWeight[0])/delta[j];
                gaussQuad4(mesh_, start+j, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];
                    const scalar wTemp = Quad4::w[gaussI] * temp;
                    const scalar dRho = wTemp*(dRho0 + (rhoGrad_neigh[j]&delta_j) - (rhoGrad_[i]&delta_i));
                    const vector dU   = wTemp*(dU0   + (UGrad_neigh[j]&delta_j)   - (UGrad_[i]&delta_i));
                    const scalar dP   = wTemp*(dP0   + (pGrad_neigh[j]&delta_j)   - (pGrad_[i]&delta_i));
                    // cell i
                    bW[i] += basisPoly<Col6X1>(delta_i, rDeltaXYZ_[i], basisMean_[i]) *
                             Col5X1(dRho, dU[0], dU[1], dU[2], dP).transpose();
                }
                temp = delta[j]*magSf[j]*sqr(vrWeight[1]);
                const scalar dRho1 = temp*(normal[j]&(rhoGrad_neigh[j]-rhoGrad_[i]));
                const vector dU1   = temp*(normal[j]&(UGrad_neigh[j]  -UGrad_[i]));
                const scalar dP1   = temp*(normal[j]&(pGrad_neigh[j]  -pGrad_[i]));
                //cell i
                bW[i] += basisPolyDn1<Col6X1>(delta_i[j], rDeltaXYZ_[i], normal[j]) *
                         Col5X1(dRho1, dU1[0], dU1[1], dU1[2], dP1).transpose();
                bW[i] += (Mat5X6
                         {
                            {d2Rho_neigh[j][0], d2Rho_neigh[j][1], d2Rho_neigh[j][2], d2Rho_neigh[j][3], d2Rho_neigh[j][4], d2Rho_neigh[j][5]},
                            {d2Ux_neigh[j][0],  d2Ux_neigh[j][1],  d2Ux_neigh[j][2],  d2Ux_neigh[j][3],  d2Ux_neigh[j][4],  d2Ux_neigh[j][5]},
                            {d2Uy_neigh[j][0],  d2Uy_neigh[j][1],  d2Uy_neigh[j][2],  d2Uy_neigh[j][3],  d2Uy_neigh[j][4],  d2Uy_neigh[j][5]},
                            {d2Uz_neigh[j][0],  d2Uz_neigh[j][1],  d2Uz_neigh[j][2],  d2Uz_neigh[j][3],  d2Uz_neigh[j][4],  d2Uz_neigh[j][5]},
                            {d2P_neigh[j][0],   d2P_neigh[j][1],   d2P_neigh[j][2],   d2P_neigh[j][3],   d2P_neigh[j][4],   d2P_neigh[j][5]}
                         } * coupledB_[patchI][j]).transpose();
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            std::vector<vector> quadPoints(9, vector::zero);
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                gaussQuad9(mesh_, start+faceI, quadPoints);
                const vector dp = magSf[faceI]*cmptMultiply(vector(1/delta[faceI], delta[faceI], pow3(delta[faceI])),
                                                            cmptMultiply(vrWeight, vrWeight));
                Mat6X5 b = Mat6X5::Zero();
                for (label gaussI = 0; gaussI != 9; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[faceI];
                    const vector dpwTemp = Quad9::w[gaussI]*dp;
                    Col6X1 Dn0_i = basisPoly<Col6X1>(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    Col6X1 Dn1_i = basisPolyDn1<Col6X1>(delta_i, rDeltaXYZ_[i], normal[faceI]);
                    // Zero Order
                    const vector U_i = U_[i] + (UGrad_[i]&delta_i) + vector(d2Var_[i].row(1)*Dn0_i,
                                                                            d2Var_[i].row(2)*Dn0_i,
                                                                            d2Var_[i].row(3)*Dn0_i);
                    const vector U_j = U_i - 2*(U_i&normal[faceI])*normal[faceI];
                    const vector dU = dpwTemp[0]*(U_j-(U_[i]+(UGrad_[i]&delta_i)));
                    b.col(1) += Dn0_i*dU[0];
                    b.col(2) += Dn0_i*dU[1];
                    b.col(3) += Dn0_i*dU[2];
                    // First Order
                    Col5X1 Dn1Var_i = d2Var_[i]*Dn1_i;
                    scalar rhoGrad = (rhoGrad_[i]&normal[faceI]) + Dn1Var_i(0);
                    vector UGrad   = (UGrad_[i]&normal[faceI])   + vector(Dn1Var_i(1), Dn1Var_i(2), Dn1Var_i(3));
                    scalar pGrad   = (pGrad_[i]&normal[faceI])   + Dn1Var_i(4);
                    rhoGrad = -dpwTemp[1]*(rhoGrad+(rhoGrad_[i]&normal[faceI]));
                    UGrad   = -dpwTemp[1]*(UGrad  +(UGrad_[i]&normal[faceI])-2*(UGrad&normal[faceI])*normal[faceI]);
                    pGrad   = -dpwTemp[1]*(pGrad  +(pGrad_[i]&normal[faceI]));
                    b += Dn1_i*Col5X1(rhoGrad, UGrad[0], UGrad[1], UGrad[2], pGrad).transpose();
                }
                // Second Order
                Col6X1 Dn2_i = basisPolyDn2<Col6X1>(rDeltaXYZ_[i], normal[faceI]);
                vector d2U = vector(d2Var_[i].row(1)*Dn2_i, d2Var_[i].row(2)*Dn2_i, d2Var_[i].row(3)*Dn2_i);
                d2U = dp[2]*(d2U-2*(d2U&normal[faceI])*normal[faceI]);
                b.col(1) += Dn2_i*d2U[0];
                b.col(2) += Dn2_i*d2U[1];
                b.col(3) += Dn2_i*d2U[2];
                bW[i] += b;
            }
        }
    }

    /* Block Gauss–Seidel iteration */
    forAll(mesh_.cells(), cellI)
    {
        if (isP0Cell_[cellI])
        {
            d2Var_[cellI] = Mat5X6::Zero();
            continue;
        }
        const labelList& cellFaces = mesh_.cells()[cellI];
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j != cellI) /* j < cellI */
                {
                    bW[cellI] += B_[faceI] * d2Var_[j].transpose();
                }
                else /* j > cellI */
                {
                    j = neighbour[faceI];
                    bW[cellI] += (d2Var_[j] * B_[faceI]).transpose();
                }
            }
        }
        d2Var_[cellI] = (rA_[cellI]*bW[cellI]).transpose();
    }

    /* Limiting coefficients of polynomials */
    limitCoefficients();
    rhoGrad_.correctBoundaryConditions();
    UGrad_.correctBoundaryConditions();
    pGrad_.correctBoundaryConditions();
    d2Rho_.correctBoundaryConditions();
    d2Ux_.correctBoundaryConditions();
    d2Uy_.correctBoundaryConditions();
    d2Uz_.correctBoundaryConditions();
    d2P_.correctBoundaryConditions();
}