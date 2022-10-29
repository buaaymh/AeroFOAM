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
            if (cmptSum(cmptMag(rhoGrad_[cellI])) < 1e-6)
            {
                isP0Cell_[cellI] = true;
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << "P0 cell [%] = " << percent << ", ";
    }

    /* Update d1Var */
    forAll(mesh_.C(), cellI)
    {
        d1Var_[cellI] << rhoGrad_[cellI][0], rhoGrad_[cellI][1], rhoGrad_[cellI][2],
                         UGrad_[cellI][0],   UGrad_[cellI][1],   UGrad_[cellI][2],
                         UGrad_[cellI][3],   UGrad_[cellI][4],   UGrad_[cellI][5],
                         UGrad_[cellI][6],   UGrad_[cellI][7],   UGrad_[cellI][8],
                         TGrad_[cellI][0],   TGrad_[cellI][1],   TGrad_[cellI][2];
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
        Col5X1 dW0 = Col5X1(rho_[j], U_[j][0], U_[j][1], U_[j][2], T_[j]) -
                     Col5X1(rho_[i], U_[i][0], U_[i][1], U_[i][2], T_[i]);
        scalar temp = 0.25*mesh_.magSf()[faceI]*sqr(vrWeight[0])/delta[faceI];
        for (label gaussI = 0; gaussI != 4; ++gaussI)
        {
            const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
            const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
            const scalar wTemp = Quad4::w[gaussI] * temp;
            Col5X1 dW = (d1Var_[j]*Col3X1(delta_j[0], delta_j[1], delta_j[2]) -
                         d1Var_[i]*Col3X1(delta_i[0], delta_i[1], delta_i[2]) + dW0) * wTemp;
            //cell i
            bW[i] += basisPolyEigen(delta_i, rDeltaXYZ_[i], basisMean_[i])*dW.transpose();
            //cell j
            bW[j] -= basisPolyEigen(delta_j, rDeltaXYZ_[j], basisMean_[j])*dW.transpose();
        }
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        temp = delta[faceI]*mesh_.magSf()[faceI]*sqr(vrWeight[1]);
        Col5X1 dW1 = (d1Var_[j]-d1Var_[i])*Col3X1(normal_[faceI][0], normal_[faceI][1], normal_[faceI][2])*temp;
        //cell i
        bW[i] += basisPolyDn1Eigen(delta_i, rDeltaXYZ_[i], normal_[faceI])*dW1.transpose();
        //cell j
        bW[j] -= basisPolyDn1Eigen(delta_j, rDeltaXYZ_[j], normal_[faceI])*dW1.transpose();
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField delta = mag(mesh_.boundary()[patchI].delta());
        const vectorField delta_i = mesh_.boundary()[patchI].Cf()-mesh_.boundary()[patchI].Cn();
        const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
        const label start = mesh_.boundary()[patchI].start();
        if (mesh_.boundary()[patchI].coupled())
        {
            const vectorField neighCn = ownerCn + mesh_.boundary()[patchI].delta();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh   = T_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField TGrad_neigh   = TGrad_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Ux_neigh = d2Ux_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uy_neigh = d2Uy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uz_neigh = d2Uz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2T_neigh = d2T_.boundaryField()[patchI].patchNeighbourField();
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                const scalar dRho0 = rho_neigh[j] - rho_[i];
                const vector dU0   = U_neigh[j] - U_[i];
                const scalar dT0   = T_neigh[j] - T_[i];
                scalar temp = 0.25*magSf[j]*sqr(vrWeight[0])/delta[j];
                gaussQuad4(start+j, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];
                    const scalar wTemp = Quad4::w[gaussI] * temp;
                    const scalar dRho = (dRho0 + (rhoGrad_neigh[j]&delta_j) - (rhoGrad_[i]&delta_i))*wTemp;
                    const vector dU   = (dU0   + (UGrad_neigh[j]&delta_j)   - (UGrad_[i]&delta_i))  *wTemp;
                    const scalar dT   = (dT0   + (TGrad_neigh[j]&delta_j)   - (TGrad_[i]&delta_i))  *wTemp;
                    //cell i
                    bW[i] += basisPolyEigen(delta_i, rDeltaXYZ_[i], basisMean_[i]) *
                             Col5X1(dRho, dU.x(), dU.y(), dU.z(), dT).transpose();
                }
                temp = delta[j]*magSf[j]*sqr(vrWeight[1]);
                const scalar dRho1 = temp*(normal[j]&(rhoGrad_neigh[j]-rhoGrad_[i]));
                const vector dU1   = temp*(normal[j]&(UGrad_neigh[j]  -UGrad_[i]));
                const scalar dT1   = temp*(normal[j]&(TGrad_neigh[j]  -TGrad_[i]));
                //cell i
                bW[i] += basisPolyDn1Eigen(delta_i[j], rDeltaXYZ_[i], normal[j]) *
                         Col5X1(dRho1, dU1.x(), dU1.y(), dU1.z(), dT1).transpose();
                bW[i] += (Mat5X6
                        {
                            {d2Rho_neigh[j][0], d2Rho_neigh[j][1], d2Rho_neigh[j][2], d2Rho_neigh[j][3], d2Rho_neigh[j][4], d2Rho_neigh[j][5]},
                            {d2Ux_neigh[j][0],  d2Ux_neigh[j][1],  d2Ux_neigh[j][2],  d2Ux_neigh[j][3],  d2Ux_neigh[j][4],  d2Ux_neigh[j][5]},
                            {d2Uy_neigh[j][0],  d2Uy_neigh[j][1],  d2Uy_neigh[j][2],  d2Uy_neigh[j][3],  d2Uy_neigh[j][4],  d2Uy_neigh[j][5]},
                            {d2Uz_neigh[j][0],  d2Uz_neigh[j][1],  d2Uz_neigh[j][2],  d2Uz_neigh[j][3],  d2Uz_neigh[j][4],  d2Uz_neigh[j][5]},
                            {d2T_neigh[j][0],   d2T_neigh[j][1],   d2T_neigh[j][2],   d2T_neigh[j][3],   d2T_neigh[j][4],   d2T_neigh[j][5]}
                        } * coupledB_[patchI][j]).transpose();
            }
        }
        else
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& TBound = T_.boundaryField()[patchI];
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                Col5X1 dW0 = Col5X1(rhoBound[faceI], UBound[faceI][0], UBound[faceI][1], UBound[faceI][2], TBound[faceI]) -
                             Col5X1(rho_[i], U_[i][0], U_[i][1], U_[i][2], T_[i]);
                scalar temp = 0.25*magSf[faceI]*sqr(vrWeight[0])/delta[faceI];
                gaussQuad4(start+faceI, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[faceI];
                    const scalar wTemp = Quad4::w[gaussI] * temp;
                    Col5X1 dW = (dW0-d1Var_[i]*Col3X1(delta_i[0], delta_i[1], delta_i[2]))*wTemp;
                    //cell i
                    bW[i] += basisPolyEigen(delta_i, rDeltaXYZ_[i], basisMean_[i])*dW.transpose();
                }
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
    TGrad_.correctBoundaryConditions();
    d2Rho_.correctBoundaryConditions();
    d2Ux_.correctBoundaryConditions();
    d2Uy_.correctBoundaryConditions();
    d2Uz_.correctBoundaryConditions();
    d2T_.correctBoundaryConditions();
}