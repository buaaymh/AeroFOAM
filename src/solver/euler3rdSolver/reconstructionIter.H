void Foam::euler3rdSolver::reconstructionIter()
{
    const Switch adaptive = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<Switch>("adaptive", false);
    isP1Cell_ = std::vector<bool>(mesh_.nCells(), false);
    if (adaptive)
    {
        label count = 0;
        volVectorField rhoGrad = fvc::grad(rhoE_);
        forAll(mesh_.cells(), cellI)
        {
            scalar dRhoE = cmptSum(cmptMag(cmptDivide(rhoGrad[cellI], rDeltaXYZ_[cellI])));
            if (dRhoE < 1e-6)
            {
                isP1Cell_[cellI] = true;
                coefs_[cellI] = Mat9X5::Zero();
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << "P1 cell [%] = " << percent << ", ";
    }
    /* Calculate b of vr linear system */
    std::vector<Mat9X5> bW(mesh_.nCells(), Mat9X5::Zero());
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const scalar dRho0  = rho_[j]  - rho_[i];
        const vector dRhoU0 = rhoU_[j] - rhoU_[i];
        const scalar dRhoE0 = rhoE_[j] - rhoE_[i];
        const Col5X1 dVars(dRho0, dRhoU0[0], dRhoU0[1], dRhoU0[2], dRhoE0);
        bW[i] += lowerb_[faceI]*dVars.transpose();
        bW[j] -= upperb_[faceI]*dVars.transpose();
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const word type = mesh_.boundary()[patchI].type();
        if (mesh_.boundary()[patchI].coupled())
        {
            const scalarField rho_neigh  = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoU_neigh = rhoU_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoE_neigh = rhoE_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1Rho_neigh  = d1Rho_.boundaryField()[patchI].patchNeighbourField();
            const tensorField d1RhoU_neigh = d1RhoU_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1RhoE_neigh = d1RhoE_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUx_neigh = d2RhoUx_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUy_neigh = d2RhoUy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUz_neigh = d2RhoUz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoE_neigh  = d2RhoE_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                const scalar dRho0  = rho_neigh[j]  - rho_[i];
                const vector dRhoU0 = rhoU_neigh[j] - rhoU_[i];
                const scalar dRhoE0 = rhoE_neigh[j] - rhoE_[i];
                const Col5X1 dVars(dRho0, dRhoU0[0], dRhoU0[1], dRhoU0[2], dRhoE0);
                bW[i] += coupledb_[patchI][j]*dVars.transpose();
                bW[i] += (Mat5X9
                {
                    {d1Rho_neigh[j][0], d1Rho_neigh[j][1], d1Rho_neigh[j][2],
                     d2Rho_neigh[j][0], d2Rho_neigh[j][1], d2Rho_neigh[j][2],
                     d2Rho_neigh[j][3], d2Rho_neigh[j][4], d2Rho_neigh[j][5]},
                    {d1RhoU_neigh[j][0],  d1RhoU_neigh[j][1],  d1RhoU_neigh[j][2],
                     d2RhoUx_neigh[j][0], d2RhoUx_neigh[j][1], d2RhoUx_neigh[j][2],
                     d2RhoUx_neigh[j][3], d2RhoUx_neigh[j][4], d2RhoUx_neigh[j][5]},
                    {d1RhoU_neigh[j][3],  d1RhoU_neigh[j][4],  d1RhoU_neigh[j][5],
                     d2RhoUy_neigh[j][0], d2RhoUy_neigh[j][1], d2RhoUy_neigh[j][2],
                     d2RhoUy_neigh[j][3], d2RhoUy_neigh[j][4], d2RhoUy_neigh[j][5]},
                    {d1RhoU_neigh[j][6],  d1RhoU_neigh[j][7],  d1RhoU_neigh[j][8],
                     d2RhoUz_neigh[j][0], d2RhoUz_neigh[j][1], d2RhoUz_neigh[j][2],
                     d2RhoUz_neigh[j][3], d2RhoUz_neigh[j][4], d2RhoUz_neigh[j][5]},
                    {d1RhoE_neigh[j][0], d1RhoE_neigh[j][1], d1RhoE_neigh[j][2],
                     d2RhoE_neigh[j][0], d2RhoE_neigh[j][1], d2RhoE_neigh[j][2],
                     d2RhoE_neigh[j][3], d2RhoE_neigh[j][4], d2RhoE_neigh[j][5]}
                } * coupledB_[patchI][j]).transpose();
            }
        }
        if (type == "symmetryPlane" || type == "symmetry" || (type == "wall" && fluidProps_.simulationType == "Euler"))
        {
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const scalarField d_ij = mag(mesh_.boundary()[patchI].delta());
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch4stQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const vector dp = vector(vrWeightSqr_[0]/d_ij[faceI],
                                         vrWeightSqr_[1]*d_ij[faceI],
                                         vrWeightSqr_[2]*pow3(d_ij[faceI]));
                Mat9X5 b = Mat9X5::Zero();
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - ownerCn[faceI];
                    const vector dpwTemp = boundaryQuad[faceI]->weight(gaussI)*dp;
                    Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    Col5X1 varDn0 = coefs_[i].transpose()*Dn0_i;
                    // Dn0
                    const vector rhoUDn0_i = rhoU_[i] + vector(varDn0(1), varDn0(2), varDn0(3));
                    const vector rhoUDn0_j = rhoUDn0_i - 2*(rhoUDn0_i&boundaryQuad[faceI]->normal(gaussI))
                                           * boundaryQuad[faceI]->normal(gaussI);
                    // * * * * * * * * //
                    // vector N = vector(-boundaryQuad[faceI]->at(gaussI)[0], -boundaryQuad[faceI]->at(gaussI)[1], 0);
                    // N /= mag(N);
                    // const vector rhoUDn0_j = rhoUDn0_i - 2*(rhoUDn0_i&N)*N;
                    // * * * * * * * * //
                    const scalar dRhoDn0  = dpwTemp[0]*varDn0(0);
                    const vector dRhoUDn0 = dpwTemp[0]*(rhoUDn0_j - rhoU_[i]);
                    const scalar dRhoEDn0 = dpwTemp[0]*varDn0(4);
                    b += Dn0_i*Col5X1(dRhoDn0, dRhoUDn0[0], dRhoUDn0[1], dRhoUDn0[2], dRhoEDn0).transpose();
                }
                bW[i] += b;
            }
        }
    }

    /* Block Gauss–Seidel iteration */
    forAll(mesh_.cells(), cellI)
    {
        if (isP1Cell_[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j != cellI) /* j < cellI */
                {
                    bW[cellI] += B_[faceI] * coefs_[j];
                }
                else /* j > cellI */
                {
                    j = neighbour[faceI];
                    bW[cellI] += B_[faceI].transpose() * coefs_[j];
                }
            }
        }
        coefs_[cellI] = rA_[cellI] * bW[cellI];
    }
    /* Limiting coefficients of polynomials */
    limitCoefficients();

    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                d1Rho_[i]  = vector(coefs_[i](0,0), coefs_[i](1,0), coefs_[i](2,0));
                d1RhoU_[i] = tensor(coefs_[i](0,1), coefs_[i](1,1), coefs_[i](2,1),
                                    coefs_[i](0,2), coefs_[i](1,2), coefs_[i](2,2),
                                    coefs_[i](0,3), coefs_[i](1,3), coefs_[i](2,3));
                d1RhoE_[i] = vector(coefs_[i](0,4), coefs_[i](1,4), coefs_[i](2,4));
                d2Rho_[i]   = symmTensor(coefs_[i](3,0), coefs_[i](4,0), coefs_[i](5,0), coefs_[i](6,0), coefs_[i](7,0), coefs_[i](8,0));
                d2RhoUx_[i] = symmTensor(coefs_[i](3,1), coefs_[i](4,1), coefs_[i](5,1), coefs_[i](6,1), coefs_[i](7,1), coefs_[i](8,1));
                d2RhoUy_[i] = symmTensor(coefs_[i](3,2), coefs_[i](4,2), coefs_[i](5,2), coefs_[i](6,2), coefs_[i](7,2), coefs_[i](8,2));
                d2RhoUz_[i] = symmTensor(coefs_[i](3,3), coefs_[i](4,3), coefs_[i](5,3), coefs_[i](6,3), coefs_[i](7,3), coefs_[i](8,3));
                d2RhoE_[i]  = symmTensor(coefs_[i](3,4), coefs_[i](4,4), coefs_[i](5,4), coefs_[i](6,4), coefs_[i](7,4), coefs_[i](8,4));
            }
        }
    }
}