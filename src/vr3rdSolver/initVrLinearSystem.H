void Foam::vr3rdSolver::initVrLinearSystem()
{
    //rDeltaXYZ and basisMean Initialize
    forAll(mesh_.C(), cellI)
    {
        evaluateBasisFunctionConst(mesh_, cellI, rDeltaXYZ_[cellI], basisMean_[cellI]);
    }
    rDeltaXYZ_.correctBoundaryConditions();
    basisMean_.correctBoundaryConditions();

    // rA, B Initialize
    const vector vrWeight = mesh_.schemesDict().subDict("vrSchemes").lookup<vector>("weightList");
    const scalarField delta = mag(mesh_.delta());
    forAll(mesh_.owner(), faceI)
    {
        std::vector<vector> quadPoints(9, vector::zero);
        gaussQuad9(mesh_, faceI, quadPoints);
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        Mat6X6 A_i  = Mat6X6::Zero();
        Mat6X6 A_j  = Mat6X6::Zero();
        Mat6X6 B_ij = Mat6X6::Zero();
        const vector dp = mesh_.magSf()[faceI]*cmptMultiply(vector(1/delta[faceI], delta[faceI], pow3(delta[faceI])),
                                                                 cmptMultiply(vrWeight, vrWeight));
        for (label gaussI = 0; gaussI != 9; ++gaussI)
        {
            const vector delta_i = quadPoints[gaussI] - mesh_.C()[i];
            const vector delta_j = quadPoints[gaussI] - mesh_.C()[j];
            const vector dpwTemp = Quad9::w[gaussI] * dp;
            symmTensor Dn0_i = basisPoly<symmTensor>(delta_i, rDeltaXYZ_[i], basisMean_[i]);
            symmTensor Dn0_j = basisPoly<symmTensor>(delta_j, rDeltaXYZ_[j], basisMean_[j]);
            symmTensor Dn1_i = basisPolyDn1<symmTensor>(delta_i, rDeltaXYZ_[i], normal_[faceI]);
            symmTensor Dn1_j = basisPolyDn1<symmTensor>(delta_j, rDeltaXYZ_[j], normal_[faceI]);
            symmTensor Dn2_i = basisPolyDn2<symmTensor>(rDeltaXYZ_[i], normal_[faceI]);
            symmTensor Dn2_j = basisPolyDn2<symmTensor>(rDeltaXYZ_[j], normal_[faceI]);
            for (label m = 0; m != 6; ++m) {
                for (label n = 0; n != 6; ++n)
                {
                    const vector n_i(Dn0_i[n], Dn1_i[n], Dn2_i[n]);
                    const vector m_i(Dn0_i[m], Dn1_i[m], Dn2_i[m]);
                    const vector n_j(Dn0_j[n], Dn1_j[n], Dn2_j[n]);
                    const vector m_j(Dn0_j[m], Dn1_j[m], Dn2_j[m]);
                    A_i(m, n) += dpwTemp&vector(cmptMultiply(n_i, m_i));
                    A_j(m, n) += dpwTemp&vector(cmptMultiply(n_j, m_j));
                    B_ij(m, n) += dpwTemp&vector(cmptMultiply(n_i, m_j));
                }    
            }
        }
        rA_[i] += A_i;
        rA_[j] += A_j;
        B_[faceI] = B_ij;
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField delta = mag(mesh_.boundary()[patchI].delta());
        const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
        const label start = mesh_.boundary()[patchI].start();
        const word type = mesh_.boundary()[patchI].type();
        std::vector<vector> quadPoints(9, vector::zero);
        if (mesh_.boundary()[patchI].coupled())
        {
            coupledB_[patchI] = std::vector<Mat6X6>(mesh_.boundary()[patchI].size(), Mat6X6::Zero());
            const vectorField rDeltaXYZ_neigh = rDeltaXYZ_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisMean_neigh = basisMean_.boundaryField()[patchI].patchNeighbourField();
            const vectorField neighCn = ownerCn + mesh_.boundary()[patchI].delta();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                gaussQuad9(mesh_, start+j, quadPoints);
                Mat6X6 A_i  = Mat6X6::Zero();
                Mat6X6 B_ij = Mat6X6::Zero();
                const vector dp = magSf[j]*cmptMultiply(vector(1/delta[j], delta[j], pow3(delta[j])),
                                                             cmptMultiply(vrWeight, vrWeight));
                for (label gaussI = 0; gaussI != 9; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];
                    const vector dpwTemp = Quad9::w[gaussI] * dp;
                    symmTensor Dn0_i = basisPoly<symmTensor>(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    symmTensor Dn0_j = basisPoly<symmTensor>(delta_j, rDeltaXYZ_neigh[j], basisMean_neigh[j]);
                    symmTensor Dn1_i = basisPolyDn1<symmTensor>(delta_i, rDeltaXYZ_[i], normal[j]);
                    symmTensor Dn1_j = basisPolyDn1<symmTensor>(delta_j, rDeltaXYZ_neigh[j], normal[j]);
                    symmTensor Dn2_i = basisPolyDn2<symmTensor>(rDeltaXYZ_[i], normal[j]);
                    symmTensor Dn2_j = basisPolyDn2<symmTensor>(rDeltaXYZ_neigh[j], normal[j]);
                    for (label m = 0; m != 6; ++m) {
                        for (label n = 0; n != 6; ++n)
                        {
                            const vector n_i(Dn0_i[n], Dn1_i[n], Dn2_i[n]);
                            const vector m_i(Dn0_i[m], Dn1_i[m], Dn2_i[m]);
                            const vector m_j(Dn0_j[m], Dn1_j[m], Dn2_j[m]);
                            A_i(m, n) += dpwTemp&vector(cmptMultiply(n_i, m_i));
                            B_ij(m, n) += dpwTemp&vector(cmptMultiply(n_i, m_j));
                        }    
                    }
                }
                rA_[i] += A_i;
                coupledB_[patchI][j] = B_ij;
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                gaussQuad9(mesh_, start+faceI, quadPoints);
                Mat6X6 A_i  = Mat6X6::Zero();
                const vector dp = magSf[faceI]*cmptMultiply(vector(1/delta[faceI], delta[faceI], pow3(delta[faceI])),
                                                            cmptMultiply(vrWeight, vrWeight));
                for (label gaussI = 0; gaussI != 9; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[faceI];
                    const vector dpwTemp = Quad9::w[gaussI] * dp;
                    symmTensor Dn0_i = basisPoly<symmTensor>(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    symmTensor Dn1_i = basisPolyDn1<symmTensor>(delta_i, rDeltaXYZ_[i], normal[faceI]);
                    symmTensor Dn2_i = basisPolyDn2<symmTensor>(rDeltaXYZ_[i], normal[faceI]);
                    for (label m = 0; m != 6; ++m) {
                        for (label n = 0; n != 6; ++n)
                        {
                            const vector n_i(Dn0_i[n], Dn1_i[n], Dn2_i[n]);
                            const vector m_i(Dn0_i[m], Dn1_i[m], Dn2_i[m]);
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
    forAll(mesh_.owner(), faceI) { gaussQuad4(mesh_, faceI, quad_[faceI]); }

    // Initialize N*h^((k+1)/2)
    forAll(mesh_.owner(), faceI)
    {
        const scalar h_k = pow(delta[faceI], 1.5);
        N_h_[mesh_.owner()[faceI]] += h_k;
        N_h_[mesh_.neighbour()[faceI]] += h_k;
    }
}