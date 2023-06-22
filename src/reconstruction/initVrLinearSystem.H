void Foam::Reconstruction::initVrLinearSystem()
{
    //rDeltaXYZ and basisMean Initialize
    forAll(mesh_.C(), cellI)
    {
        evaluateBasisFunctionConst(mesh_, cellI, rDeltaXYZ_[cellI], basisMean_[cellI], polyNorm_[cellI]);
    }
    rDeltaXYZ_.correctBoundaryConditions();
    basisMean_.correctBoundaryConditions();
    
    std::unique_ptr<Face> faceQuad;
    // Interior face quadrature points Init
    faceQuad_.reserve(mesh_.nInternalFaces());
    forAll(mesh_.owner(), faceI)
    {
        build2ndFace(mesh_, faceI, faceQuad);
        faceQuad_.emplace_back(std::move(faceQuad));
    }
    // Boundary face quadrature points Init
    const Switch curvedWall = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<Switch>("curvedWall", false);
    forAll(mesh_.boundary(), patchI)
    {
        const label start = mesh_.boundary()[patchI].start();
        const word type = mesh_.boundary()[patchI].type();
        if (mesh_.boundary()[patchI].coupled())
        {
            patch2ndQuad_[patchI] = std::vector<std::unique_ptr<Face>>();
            patch2ndQuad_[patchI].reserve(mesh_.boundary()[patchI].size());
            forAll(mesh_.boundary()[patchI], faceI)
            {
                build2ndFace(mesh_, start+faceI, faceQuad);
                patch2ndQuad_[patchI].emplace_back(std::move(faceQuad));
            }
        }
        if (type == "symmetryPlane" || type == "symmetry" || (type == "wall" && !curvedWall))
        {
            patch2ndQuad_[patchI] = std::vector<std::unique_ptr<Face>>();
            patch4stQuad_[patchI] = std::vector<std::unique_ptr<Face>>();
            patch2ndQuad_[patchI].reserve(mesh_.boundary()[patchI].size());
            patch4stQuad_[patchI].reserve(mesh_.boundary()[patchI].size());
            forAll(mesh_.boundary()[patchI], faceI)
            {
                build2ndFace(mesh_, start+faceI, faceQuad);
                patch2ndQuad_[patchI].emplace_back(std::move(faceQuad));
                build4stFace(mesh_, start+faceI, faceQuad);
                patch4stQuad_[patchI].emplace_back(std::move(faceQuad));
            }
        }
        if (type == "wall" && curvedWall)
        {
            patch2ndQuad_[patchI] = std::vector<std::unique_ptr<Face>>();
            patch4stQuad_[patchI] = std::vector<std::unique_ptr<Face>>();
            patch2ndQuad_[patchI].reserve(mesh_.boundary()[patchI].size());
            patch4stQuad_[patchI].reserve(mesh_.boundary()[patchI].size());
            forAll(mesh_.boundary()[patchI], faceI)
            {

                build2ndCurvedBoundary(mesh_, patchI, start+faceI, faceQuad);
                patch2ndQuad_[patchI].emplace_back(std::move(faceQuad));
                build4stCurvedBoundary(mesh_, patchI, start+faceI, faceQuad);
                patch4stQuad_[patchI].emplace_back(std::move(faceQuad));
            }
        }
    }

    const vector vrWeight = mesh_.schemes().subDict("vrSchemes").lookup<vector>("weightList");
    vrWeightSqr_ = cmptMultiply(vrWeight, vrWeight);
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
        const vector dp = (1.0/d_ij_[faceI])*vector(vrWeightSqr_[0],
                                                    vrWeightSqr_[1]*sqr(d_ij_[faceI]),
                                                    vrWeightSqr_[2]*pow4(d_ij_[faceI]));
        for (label gaussI = 0; gaussI != faceQuad->size(); ++gaussI)
        {
            const vector delta_i = faceQuad->at(gaussI) - mesh_.C()[i];
            const vector delta_j = faceQuad->at(gaussI) - mesh_.C()[j];
            const vector dpwTemp = faceQuad->weight(gaussI) * dp;
            Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
            Col9X1 Dn0_j = polynomialDn0(delta_j, rDeltaXYZ_[j], basisMean_[j]);
            Mat9X3 Dn1_i = polynomialGrad(delta_i, rDeltaXYZ_[i]);
            Mat9X3 Dn1_j = polynomialGrad(delta_j, rDeltaXYZ_[j]);
            Col9X1 Dn2_i = Col9X1(0, 0, 0,
                                rDeltaXYZ_[i].x()*rDeltaXYZ_[i].x(),
                                rDeltaXYZ_[i].x()*rDeltaXYZ_[i].y(),
                                rDeltaXYZ_[i].x()*rDeltaXYZ_[i].z(),
                                rDeltaXYZ_[i].y()*rDeltaXYZ_[i].y(),
                                rDeltaXYZ_[i].y()*rDeltaXYZ_[i].z(),
                                rDeltaXYZ_[i].z()*rDeltaXYZ_[i].z());
            Col9X1 Dn2_j = Col9X1(0, 0, 0,
                                rDeltaXYZ_[j].x()*rDeltaXYZ_[j].x(),
                                rDeltaXYZ_[j].x()*rDeltaXYZ_[j].y(),
                                rDeltaXYZ_[j].x()*rDeltaXYZ_[j].z(),
                                rDeltaXYZ_[j].y()*rDeltaXYZ_[j].y(),
                                rDeltaXYZ_[j].y()*rDeltaXYZ_[j].z(),
                                rDeltaXYZ_[j].z()*rDeltaXYZ_[j].z());
            for (label m = 0; m != 9; ++m) {
                for (label n = 0; n != 9; ++n)
                {
                    A_i(m, n)  += dpwTemp&vector(Dn0_i(n)*Dn0_i(m), Dn1_i.row(n).dot(Dn1_i.row(m)), Dn2_i(n)*Dn2_i(m));
                    A_j(m, n)  += dpwTemp&vector(Dn0_j(n)*Dn0_j(m), Dn1_j.row(n).dot(Dn1_j.row(m)), Dn2_j(n)*Dn2_j(m));
                    B_ij(m, n) += dpwTemp&vector(Dn0_i(n)*Dn0_j(m), Dn1_i.row(n).dot(Dn1_j.row(m)), Dn2_i(n)*Dn2_j(m));
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
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
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
                const vector dp = (1.0/d_ij[j])*vector(vrWeightSqr_[0],
                                                       vrWeightSqr_[1]*sqr(d_ij[j]),
                                                       vrWeightSqr_[2]*pow4(d_ij[j]));
                for (label gaussI = 0; gaussI != faceQuad->size(); ++gaussI)
                {
                    const vector delta_i = faceQuad->at(gaussI) - ownerCn[j];
                    const vector delta_j = faceQuad->at(gaussI) - neighCn[j];
                    const vector dpwTemp = faceQuad->weight(gaussI) * dp;
                    Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    Col9X1 Dn0_j = polynomialDn0(delta_j, rDeltaXYZ_neigh[j], basisMean_neigh[j]);
                    Mat9X3 Dn1_i = polynomialGrad(delta_i, rDeltaXYZ_[i]);
                    Mat9X3 Dn1_j = polynomialGrad(delta_j, rDeltaXYZ_neigh[j]);
                    Col9X1 Dn2_i = Col9X1(0, 0, 0,
                                        rDeltaXYZ_[i].x()*rDeltaXYZ_[i].x(),
                                        rDeltaXYZ_[i].x()*rDeltaXYZ_[i].y(),
                                        rDeltaXYZ_[i].x()*rDeltaXYZ_[i].z(),
                                        rDeltaXYZ_[i].y()*rDeltaXYZ_[i].y(),
                                        rDeltaXYZ_[i].y()*rDeltaXYZ_[i].z(),
                                        rDeltaXYZ_[i].z()*rDeltaXYZ_[i].z());
                    Col9X1 Dn2_j = Col9X1(0, 0, 0,
                                        rDeltaXYZ_neigh[j].x()*rDeltaXYZ_neigh[j].x(),
                                        rDeltaXYZ_neigh[j].x()*rDeltaXYZ_neigh[j].y(),
                                        rDeltaXYZ_neigh[j].x()*rDeltaXYZ_neigh[j].z(),
                                        rDeltaXYZ_neigh[j].y()*rDeltaXYZ_neigh[j].y(),
                                        rDeltaXYZ_neigh[j].y()*rDeltaXYZ_neigh[j].z(),
                                        rDeltaXYZ_neigh[j].z()*rDeltaXYZ_neigh[j].z());
                    for (label m = 0; m != 9; ++m) {
                        for (label n = 0; n != 9; ++n)
                        {
                            A_i(m, n)  += dpwTemp&vector(Dn0_i(n)*Dn0_i(m), Dn1_i.row(n).dot(Dn1_i.row(m)), Dn2_i(n)*Dn2_i(m));
                            B_ij(m, n) += dpwTemp&vector(Dn0_i(n)*Dn0_j(m), Dn1_i.row(n).dot(Dn1_j.row(m)), Dn2_i(n)*Dn2_j(m));
                        }    
                    }
                    b_i += dpwTemp[0] * Dn0_i;
                }
                rA_[i] += A_i;
                coupledB_[patchI][j] = B_ij;
                coupledb_[patchI][j] = b_i;
            }
        }
        if (type == "symmetryPlane" || type == "symmetry" || type == "wall")
        {
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch4stQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                Mat9X9 A_i = Mat9X9::Zero();
                const vector dp = (0.5/d_ij[faceI])*vector(vrWeightSqr_[0], 0, 0);
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - ownerCn[faceI];
                    const vector dpwTemp = boundaryQuad[faceI]->weight(gaussI) * dp;
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

    // Initialize N*h^((k+1)/2)
    forAll(mesh_.owner(), faceI)
    {
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        N_h_[i] = max(N_h_[i], d_ij_[faceI]);;
        N_h_[j] = max(N_h_[j], d_ij_[faceI]);;
    }
    N_h_ = pow(N_h_, 1.5);
}