void Foam::eulerConsVar3rdSolver::reconstructionIter()
{
    /* Detect adaptive P0 cells */
    const Switch adaptive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("adaptive", false);
    if (adaptive)
    {
        label count = 0;
        isP0Cell_ = std::vector<bool>(mesh_.nCells(), false);
        forAll(mesh_.C(), cellI)
        {
            if (cmptSum(cmptMag(rhoEGrad_[cellI])) < 1e-6)
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
        const scalar dRho0  = rho_[j]  - rho_[i];
        const vector dRhoU0 = rhoU_[j] - rhoU_[i];
        const scalar dRhoE0 = rhoE_[j] - rhoE_[i];
        scalar temp = mesh_.magSf()[faceI]*sqr(vrWeight[0])/delta[faceI];
        for (label gaussI = 0; gaussI != 4; ++gaussI)
        {
            const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
            const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
            const scalar wTemp = Quad4::w[gaussI] * temp;
            const scalar dRho  = wTemp*(dRho0  + (rhoGrad_[j]&delta_j)  - (rhoGrad_[i]&delta_i));
            const vector dRhoU = wTemp*(dRhoU0 + (rhoUGrad_[j]&delta_j) - (rhoUGrad_[i]&delta_i));
            const scalar dRhoE = wTemp*(dRhoE0 + (rhoEGrad_[j]&delta_j) - (rhoEGrad_[i]&delta_i));
            Col5X1 dVars(dRho, dRhoU[0], dRhoU[1], dRhoU[2], dRhoE);
            //cell i
            bW[i] += basisPoly<Col6X1>(delta_i, rDeltaXYZ_[i], basisMean_[i])*dVars.transpose();
            //cell j
            bW[j] -= basisPoly<Col6X1>(delta_j, rDeltaXYZ_[j], basisMean_[j])*dVars.transpose();
        }
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        temp = delta[faceI]*mesh_.magSf()[faceI]*sqr(vrWeight[1]);
        const scalar dRho1  = temp*(normal_[faceI]&(rhoGrad_[j]  - rhoGrad_[i]));
        const vector dRhoU1 = temp*(normal_[faceI]&(rhoUGrad_[j] - rhoUGrad_[i]));
        const scalar dRhoE1 = temp*(normal_[faceI]&(rhoEGrad_[j] - rhoEGrad_[i]));
        Col5X1 dGrad(dRho1, dRhoU1[0], dRhoU1[1], dRhoU1[2], dRhoE1);
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
        const label start = mesh_.boundary()[patchI].start();
        if (mesh_.boundary()[patchI].coupled())
        {
            const vectorField neighCn = ownerCn + mesh_.boundary()[patchI].delta();
            const scalarField rho_neigh  = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoU_neigh = rhoU_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoE_neigh = rhoE_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh  = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField rhoUGrad_neigh = rhoUGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoEGrad_neigh = rhoEGrad_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUx_neigh = d2RhoUx_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUy_neigh = d2RhoUy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUz_neigh = d2RhoUz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoE_neigh  = d2RhoE_.boundaryField()[patchI].patchNeighbourField();
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                const scalar dRho0  = rho_neigh[j]  - rho_[i];
                const vector dRhoU0 = rhoU_neigh[j] - rhoU_[i];
                const scalar dRhoE0 = rhoE_neigh[j] - rhoE_[i];
                scalar temp = magSf[j]*sqr(vrWeight[0])/delta[j];
                gaussQuad4(mesh_, start+j, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];
                    const scalar wTemp = Quad4::w[gaussI] * temp;
                    const scalar dRho  = wTemp*(dRho0  + (rhoGrad_neigh[j]&delta_j)  - (rhoGrad_[i]&delta_i));
                    const vector dRhoU = wTemp*(dRhoU0 + (rhoUGrad_neigh[j]&delta_j) - (rhoUGrad_[i]&delta_i));
                    const scalar dRhoE = wTemp*(dRhoE0 + (rhoEGrad_neigh[j]&delta_j) - (rhoEGrad_[i]&delta_i));
                    // cell i
                    bW[i] += basisPoly<Col6X1>(delta_i, rDeltaXYZ_[i], basisMean_[i]) *
                             Col5X1(dRho, dRhoU[0], dRhoU[1], dRhoU[2], dRhoE).transpose();
                }
                temp = delta[j]*magSf[j]*sqr(vrWeight[1]);
                const scalar dRho1  = temp*(normal[j]&(rhoGrad_neigh[j]  - rhoGrad_[i]));
                const vector dRhoU1 = temp*(normal[j]&(rhoUGrad_neigh[j] - rhoUGrad_[i]));
                const scalar dRhoE1 = temp*(normal[j]&(rhoEGrad_neigh[j] - rhoEGrad_[i]));
                //cell i
                bW[i] += basisPolyDn1<Col6X1>(delta_i[j], rDeltaXYZ_[i], normal[j]) *
                         Col5X1(dRho1, dRhoU1[0], dRhoU1[1], dRhoU1[2], dRhoE1).transpose();
                bW[i] += (Mat5X6
                        {
                            {d2Rho_neigh[j][0],   d2Rho_neigh[j][1],   d2Rho_neigh[j][2],   d2Rho_neigh[j][3],   d2Rho_neigh[j][4],   d2Rho_neigh[j][5]},
                            {d2RhoUx_neigh[j][0], d2RhoUx_neigh[j][1], d2RhoUx_neigh[j][2], d2RhoUx_neigh[j][3], d2RhoUx_neigh[j][4], d2RhoUx_neigh[j][5]},
                            {d2RhoUy_neigh[j][0], d2RhoUy_neigh[j][1], d2RhoUy_neigh[j][2], d2RhoUy_neigh[j][3], d2RhoUy_neigh[j][4], d2RhoUy_neigh[j][5]},
                            {d2RhoUz_neigh[j][0], d2RhoUz_neigh[j][1], d2RhoUz_neigh[j][2], d2RhoUz_neigh[j][3], d2RhoUz_neigh[j][4], d2RhoUz_neigh[j][5]},
                            {d2RhoE_neigh[j][0],  d2RhoE_neigh[j][1],  d2RhoE_neigh[j][2],  d2RhoE_neigh[j][3],  d2RhoE_neigh[j][4],  d2RhoE_neigh[j][5]}
                        } * coupledB_[patchI][j]).transpose();
            }
        }
        else
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& rhoUBound = rhoU_.boundaryField()[patchI];
            const fvPatchScalarField& rhoEBound = rhoE_.boundaryField()[patchI];
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const scalar dRho0  = rhoBound[faceI]  - rho_[i];
                const vector dRhoU0 = rhoUBound[faceI] - rhoU_[i];
                const scalar dRhoE0 = rhoEBound[faceI] - rhoE_[i];
                scalar temp = magSf[faceI]*sqr(vrWeight[0])/delta[faceI];
                gaussQuad4(mesh_, start+faceI, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[faceI];
                    const scalar wTemp = Quad4::w[gaussI] * temp;
                    const scalar dRho  = wTemp*(dRho0  - (rhoGrad_[i]&delta_i));
                    const vector dRhoU = wTemp*(dRhoU0 - (rhoUGrad_[i]&delta_i));
                    const scalar dRhoE = wTemp*(dRhoE0 - (rhoEGrad_[i]&delta_i));
                    //cell i
                    bW[i] += basisPoly<Col6X1>(delta_i, rDeltaXYZ_[i], basisMean_[i]) *
                             Col5X1(dRho, dRhoU[0], dRhoU[1], dRhoU[2], dRhoE).transpose();
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
    rhoUGrad_.correctBoundaryConditions();
    rhoEGrad_.correctBoundaryConditions();
    d2Rho_.correctBoundaryConditions();
    d2RhoUx_.correctBoundaryConditions();
    d2RhoUy_.correctBoundaryConditions();
    d2RhoUz_.correctBoundaryConditions();
    d2RhoE_.correctBoundaryConditions();
}