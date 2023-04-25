void Foam::euler3rdSolver::initVrLinearSystem()
{
    //rDeltaXYZ and basisMean Initialize
    forAll(mesh_.C(), cellI)
    {
        evaluateBasisFunctionConst(mesh_, cellI, rDeltaXYZ_[cellI], basisMean_[cellI], polyNorm_[cellI]);
    }
    rDeltaXYZ_.correctBoundaryConditions();
    basisMean_.correctBoundaryConditions();

    const vector vrWeight = mesh_.schemes().subDict("vrSchemes").lookup<vector>("weightList");
    vrWeightSqr_ = cmptMultiply(vrWeight, vrWeight);
    std::unique_ptr<Face> faceQuad;
    forAll(mesh_.owner(), faceI)
    {
        // rA, B, b Initialize
        build4stFace(mesh_, faceI, faceQuad);
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        Mat9X9 A_i  = Mat9X9::Zero();
        Mat9X9 A_j  = Mat9X9::Zero();
        Mat9X9 B_ij = Mat9X9::Zero();
        Col9X1 b_i  = Col9X1::Zero();
        Col9X1 b_j  = Col9X1::Zero();
        const vector dp = vector(vrWeightSqr_[0]/d_ij_[faceI],
                                 vrWeightSqr_[1]*d_ij_[faceI],
                                 vrWeightSqr_[2]*pow3(d_ij_[faceI]));
        for (label gaussI = 0; gaussI != faceQuad->size(); ++gaussI)
        {
            const vector delta_i = faceQuad->at(gaussI) - mesh_.C()[i];
            const vector delta_j = faceQuad->at(gaussI) - mesh_.C()[j];
            const vector dpwTemp = faceQuad->weight(gaussI) * dp;
            Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
            Col9X1 Dn0_j = polynomialDn0(delta_j, rDeltaXYZ_[j], basisMean_[j]);
            Col9X1 Dn1_i = polynomialDn1(delta_i, rDeltaXYZ_[i], normal_[faceI]);
            Col9X1 Dn1_j = polynomialDn1(delta_j, rDeltaXYZ_[j], normal_[faceI]);
            Col9X1 Dn2_i = polynomialDn2(rDeltaXYZ_[i], normal_[faceI]);
            Col9X1 Dn2_j = polynomialDn2(rDeltaXYZ_[j], normal_[faceI]);
            for (label m = 0; m != 9; ++m) {
                for (label n = 0; n != 9; ++n)
                {
                    const vector n_i(Dn0_i[n], Dn1_i[n], Dn2_i[n]);
                    const vector m_i(Dn0_i[m], Dn1_i[m], Dn2_i[m]);
                    const vector n_j(Dn0_j[n], Dn1_j[n], Dn2_j[n]);
                    const vector m_j(Dn0_j[m], Dn1_j[m], Dn2_j[m]);
                    A_i(m, n)  += dpwTemp&vector(cmptMultiply(n_i, m_i));
                    A_j(m, n)  += dpwTemp&vector(cmptMultiply(n_j, m_j));
                    B_ij(m, n) += dpwTemp&vector(cmptMultiply(n_i, m_j));
                }    
            }
            b_i += dpwTemp[0] * Dn0_i;
            b_j += dpwTemp[0] * Dn0_j;
        }
        rA_[i] += A_i;
        rA_[j] += A_j;
        B_[faceI] = B_ij;
        lowerb_[faceI] = b_i;
        upperb_[faceI] = b_j;
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField d_ij = mag(mesh_.boundary()[patchI].delta());
        const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
        const label start = mesh_.boundary()[patchI].start();
        const word type = mesh_.boundary()[patchI].type();
        if (mesh_.boundary()[patchI].coupled())
        {
            coupledB_[patchI] = std::vector<Mat9X9>(mesh_.boundary()[patchI].size(), Mat9X9::Zero());
            coupledb_[patchI] = std::vector<Col9X1>(mesh_.boundary()[patchI].size(), Col9X1::Zero());
            const vectorField rDeltaXYZ_neigh = rDeltaXYZ_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisMean_neigh = basisMean_.boundaryField()[patchI].patchNeighbourField();
            const vectorField neighCn = ownerCn + mesh_.boundary()[patchI].delta();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                build4stFace(mesh_, start+j, faceQuad);
                Mat9X9 A_i  = Mat9X9::Zero();
                Mat9X9 B_ij = Mat9X9::Zero();
                Col9X1 b_i  = Col9X1::Zero();
                const vector dp = vector(vrWeightSqr_[0]/d_ij[j],
                                         vrWeightSqr_[1]*d_ij[j],
                                         vrWeightSqr_[2]*pow3(d_ij[j]));
                for (label gaussI = 0; gaussI != faceQuad->size(); ++gaussI)
                {
                    const vector delta_i = faceQuad->at(gaussI) - ownerCn[j];
                    const vector delta_j = faceQuad->at(gaussI) - neighCn[j];
                    const vector dpwTemp = faceQuad->weight(gaussI) * dp;
                    Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    Col9X1 Dn0_j = polynomialDn0(delta_j, rDeltaXYZ_neigh[j], basisMean_neigh[j]);
                    Col9X1 Dn1_i = polynomialDn1(delta_i, rDeltaXYZ_[i], normal[j]);
                    Col9X1 Dn1_j = polynomialDn1(delta_j, rDeltaXYZ_neigh[j], normal[j]);
                    Col9X1 Dn2_i = polynomialDn2(rDeltaXYZ_[i], normal[j]);
                    Col9X1 Dn2_j = polynomialDn2(rDeltaXYZ_neigh[j], normal[j]);
                    for (label m = 0; m != 9; ++m) {
                        for (label n = 0; n != 9; ++n)
                        {
                            const vector n_i(Dn0_i[n], Dn1_i[n], Dn2_i[n]);
                            const vector m_i(Dn0_i[m], Dn1_i[m], Dn2_i[m]);
                            const vector m_j(Dn0_j[m], Dn1_j[m], Dn2_j[m]);
                            A_i(m, n)  += dpwTemp&vector(cmptMultiply(n_i, m_i));
                            B_ij(m, n) += dpwTemp&vector(cmptMultiply(n_i, m_j));
                        }    
                    }
                    b_i += dpwTemp[0] * Dn0_i;
                }
                rA_[i] += A_i;
                coupledB_[patchI][j] = B_ij;
                coupledb_[patchI][j] = b_i;
            }
        }
        if (type == "symmetryPlane" || type == "symmetry" || (type == "wall" && fluidProps_.simulationType == "Euler"))
        {
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                build4stFace(mesh_, start+faceI, faceQuad);
                Mat9X9 A_i = Mat9X9::Zero();
                const vector dp = vector(vrWeightSqr_[0]/d_ij[faceI],
                                         vrWeightSqr_[1]*d_ij[faceI],
                                         vrWeightSqr_[2]*pow3(d_ij[faceI]));
                for (label gaussI = 0; gaussI != faceQuad->size(); ++gaussI)
                {
                    const vector delta_i = faceQuad->at(gaussI) - ownerCn[faceI];
                    const vector dpwTemp = faceQuad->weight(gaussI) * dp;
                    Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    for (label m = 0; m != 9; ++m)
                    {
                        for (label n = 0; n != 9; ++n)
                        {
                            const vector n_i(Dn0_i[n], 0, 0);
                            const vector m_i(Dn0_i[m], 0, 0);
                            A_i(m, n) += dpwTemp&vector(cmptMultiply(n_i, m_i));
                        }
                    }
                }
                rA_[i] += A_i;
            }
        }
    }
    forAll(mesh_.C(), cellI) rA_[cellI] = rA_[cellI].inverse();

    // Quadrature points Init
    faceQuad_.reserve(mesh_.nInternalFaces());
    forAll(mesh_.owner(), faceI)
    {
        build2ndFace(mesh_, faceI, faceQuad);
        faceQuad_.emplace_back(std::move(faceQuad));
    }

    // Initialize N*h^((k+1)/2)
    forAll(mesh_.owner(), faceI)
    {
        const scalar h_k = pow(d_ij_[faceI], 1.5);
        N_h_[mesh_.owner()[faceI]] += h_k;
        N_h_[mesh_.neighbour()[faceI]] += h_k;
    }
}