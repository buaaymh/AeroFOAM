void Foam::vrScheme::rLengthScaleInit()
{
    forAll(mesh_.C(), cellI)
    {
        vector pointMin = mesh_.C()[cellI];
        vector pointMax = mesh_.C()[cellI];
        const UList<label>& cellPointsId = mesh_.cellPoints()[cellI];
        forAll(cellPointsId, pointI)
        {
            pointMin = min(pointMin, mesh_.points()[cellPointsId[pointI]]);
            pointMax = max(pointMax, mesh_.points()[cellPointsId[pointI]]);
        }
        rLengthScale_[cellI] = cmptDivide(vector::one, 0.5*(pointMax-pointMin));
    }
    rLengthScale_.correctBoundaryConditions();
}

void Foam::vrScheme::basisConstInit()
{
    forAll(mesh_.C(), cellI)
    {
        std::vector<vector> quadPoints(8, vector::zero);
        gaussHex8(cellI, quadPoints);
        for (label i = 0; i != 8; ++i)
        {
            const vector B = cmptMultiply(quadPoints[i]-mesh_.C()[cellI], rLengthScale_[cellI]);
            symmTensor temp
            (
                sqr(B.x()), B.x()*B.y(), B.x()*B.z(),
                             sqr(B.y()), B.y()*B.z(),
                                           sqr(B.z())
            );
            basisConst_[cellI] += Hex8::w[i] * temp;
        }
        basisConst_[cellI] *= 0.125;
    }
    basisConst_.correctBoundaryConditions();
}

void Foam::vrScheme::matInit()
{
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        std::vector<vector> quadPoints(9, vector::zero);
        gaussQuad9(faceI, quadPoints);
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        vrScheme::Matrix A_i  = vrScheme::Matrix::Zero();
        vrScheme::Matrix A_j  = vrScheme::Matrix::Zero();
        vrScheme::Matrix B_ij = vrScheme::Matrix::Zero();
        const vector normal = mesh_.Sf()[faceI]/mesh_.magSf()[faceI];
        const vector dp = cmptDivide(vector(1/delta_[faceI], delta_[faceI], pow3(delta_[faceI])),
                                     cmptMultiply(vrWeight_, vrWeight_)) * mesh_.magSf()[faceI] * 0.25;
        for (label gaussI = 0; gaussI != 9; ++gaussI)
        {
            vrScheme::Matrix tempA_i  = vrScheme::Matrix::Zero();
            vrScheme::Matrix tempA_j  = vrScheme::Matrix::Zero();
            vrScheme::Matrix tempB_ij = vrScheme::Matrix::Zero();
            const vector delta_i = quadPoints[gaussI] - mesh_.C()[i];
            const vector delta_j = quadPoints[gaussI] - mesh_.C()[j];
            const vector dpwTemp = Quad9::w[gaussI] * dp;
            symmTensor Dn0_i = Foam::basisPoly(delta_i, rLengthScale_[i], basisConst_[i]);
            symmTensor Dn1_i = Foam::basisPolyDn1(delta_i, rLengthScale_[i], normal);
            symmTensor Dn2_i = Foam::basisPolyDn2(rLengthScale_[i], normal);
            symmTensor Dn0_j = Foam::basisPoly(delta_j, rLengthScale_[j], basisConst_[j]);
            symmTensor Dn1_j = Foam::basisPolyDn1(delta_j, rLengthScale_[j], normal);
            symmTensor Dn2_j = Foam::basisPolyDn2(rLengthScale_[j], normal);
            for (label m = 0; m != 6; ++m) {
                for (label n = 0; n != 6; ++n)
                {
                    const vector n_i(Dn0_i[n], Dn1_i[n], Dn2_i[n]);
                    const vector m_i(Dn0_i[m], Dn1_i[m], Dn2_i[m]);
                    const vector n_j(Dn0_j[n], Dn1_j[n], Dn2_j[n]);
                    const vector m_j(Dn0_j[m], Dn1_j[m], Dn2_j[m]);
                    tempA_i(m, n)  = dpwTemp&vector(cmptMultiply(n_i, m_i));
                    tempA_j(m, n)  = dpwTemp&vector(cmptMultiply(n_j, m_j));
                    tempB_ij(m, n) = dpwTemp&vector(cmptMultiply(n_i, m_j));
                }    
            }
            A_i  += tempA_i;
            A_j  += tempA_j;
            B_ij += tempB_ij;
        }
        rA_[i] += A_i;
        rA_[j] += A_j;
        B_[faceI] = B_ij;
    }
    forAll(mesh_.C(), cellI) rA_[cellI] = rA_[cellI].inverse();
}