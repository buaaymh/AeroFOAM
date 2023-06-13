void Foam::evaluateBasisFunctionConst
(
    const fvMesh& mesh,
    const label& cellI,
    vector& rDeltaXYZ,
    symmTensor& basisMean,
    Col9X1& polyNorm
)
{
    //rDeltaXYZ Initialize
    vector pointMin = mesh.C()[cellI];
    vector pointMax = mesh.C()[cellI];
    const UList<label>& cellPointsId = mesh.cellPoints()[cellI];
    forAll(cellPointsId, pointI)
    {
        pointMin = min(pointMin, mesh.points()[cellPointsId[pointI]]);
        pointMax = max(pointMax, mesh.points()[cellPointsId[pointI]]);
    }
    rDeltaXYZ = cmptDivide(vector::one, 0.5*(pointMax-pointMin));
    // basisMean Initialize
    std::vector<scalar> weights;
    std::vector<vector> quadPoints;
    build2ndCell(mesh, cellI, weights, quadPoints);
    for (size_t i = 0; i != weights.size(); ++i)
    {
        const vector B = cmptMultiply(quadPoints[i]-mesh.C()[cellI], rDeltaXYZ);
        symmTensor temp
        (
            sqr(B.x()), B.x()*B.y(), B.x()*B.z(),
            sqr(B.y()), B.y()*B.z(),
            sqr(B.z())
        );
        basisMean += weights[i] * temp;
    }
    basisMean /= mesh.V()[cellI];
    // polyNorm Initialize
    vector B = cmptMultiply(max(mesh.C()[cellI]-pointMin, pointMax-mesh.C()[cellI]), rDeltaXYZ);
    polyNorm = Col9X1
    (
        B.x(), B.y(), B.z(),
        0.5*max(sqr(B.x())-mag(basisMean.xx()), mag(basisMean.xx())), B.x()*B.y()+mag(basisMean.xy()), B.x()*B.z()+mag(basisMean.xz()),
        0.5*max(sqr(B.y())-mag(basisMean.yy()), mag(basisMean.yy())), B.y()*B.z()+mag(basisMean.yz()),
        0.5*max(sqr(B.z())-mag(basisMean.zz()), mag(basisMean.zz()))
    );
}

Col9X1 Foam::polynomialDn0
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    return Col9X1
    (
        B.x(), B.y(), B.z(),
        0.5*(sqr(B.x())-basisMean.xx()), B.x()*B.y()-basisMean.xy(), B.x()*B.z()-basisMean.xz(),
        0.5*(sqr(B.y())-basisMean.yy()), B.y()*B.z()-basisMean.yz(),
        0.5*(sqr(B.z())-basisMean.zz())
    );
}

Col9X1 Foam::polynomialDn1
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const vector& normal
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    const vector rLn = cmptMultiply(rDeltaXYZ, normal);
    return Col9X1
    (
        rLn.x(), rLn.y(), rLn.z(),
        B.x()*rLn.x(), B.y()*rLn.x()+B.x()*rLn.y(), B.z()*rLn.x()+B.x()*rLn.z(),
        B.y()*rLn.y(), B.z()*rLn.y()+B.y()*rLn.z(),
        B.z()*rLn.z()
    );
}

Col9X1 Foam::polynomialDn2
(
    const vector& rDeltaXYZ,
    const vector& normal
)
{
    const vector Bn2 = cmptMultiply(rDeltaXYZ, normal);
    return Col9X1
    (
        0, 0, 0,
        sqr(Bn2.x()), 2*Bn2.x()*Bn2.y(), 2*Bn2.x()*Bn2.z(),
        sqr(Bn2.y()), 2*Bn2.y()*Bn2.z(),
        sqr(Bn2.z())
    );
}

Mat9X3 Foam::polynomialGrad
(
    const vector& delta,
    const vector& rDeltaXYZ
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    Mat9X3 polyGrad;
    polyGrad.col(0) << rDeltaXYZ.x(), 0, 0,
                       B.x()*rDeltaXYZ.x(), B.y()*rDeltaXYZ.x(), B.z()*rDeltaXYZ.x(),
                       0, 0,
                       0;
    polyGrad.col(1) << 0, rDeltaXYZ.y(), 0,
                       0, B.x()*rDeltaXYZ.y(), 0,
                       B.y()*rDeltaXYZ.y(), B.z()*rDeltaXYZ.y(),
                       0;
    polyGrad.col(2) << 0, 0, rDeltaXYZ.z(),
                       0, 0, B.x()*rDeltaXYZ.z(),
                       0, B.y()*rDeltaXYZ.z(),
                       B.z()*rDeltaXYZ.z();
    return polyGrad;
}

Mat5X3 Foam::transformD1Coefs
(
    const Mat5X9& coefs,
    const vector& delta,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
)
{
    const vector rDeltaXYZ_ij = cmptDivide(rDeltaXYZ_j, rDeltaXYZ_i);
    const vector grad = vector(cmptMultiply(delta, rDeltaXYZ_j));
    Mat5X3 d1Coefs_ij = coefs.block<5,3>(0,0);
    d1Coefs_ij.col(0) += grad[0]*coefs.col(3) + grad[1]*coefs.col(4) + grad[2]*coefs.col(5);
    d1Coefs_ij.col(1) += grad[0]*coefs.col(4) + grad[1]*coefs.col(5) + grad[2]*coefs.col(7);
    d1Coefs_ij.col(2) += grad[0]*coefs.col(5) + grad[1]*coefs.col(7) + grad[2]*coefs.col(8);
    d1Coefs_ij.col(0) *= rDeltaXYZ_ij.x();
    d1Coefs_ij.col(1) *= rDeltaXYZ_ij.y();
    d1Coefs_ij.col(2) *= rDeltaXYZ_ij.z();
    return d1Coefs_ij;
}

Mat5X6 Foam::transformD2Coefs
(
    const Mat5X6& d2Coefs,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
)
{
    const vector rDeltaXYZ_ij = cmptDivide(rDeltaXYZ_j, rDeltaXYZ_i);
    Mat5X6 d2Coefs_ij = d2Coefs;
    d2Coefs_ij.col(0) *= sqr(rDeltaXYZ_ij.x());
    d2Coefs_ij.col(1) *= rDeltaXYZ_ij.x()*rDeltaXYZ_ij.y();
    d2Coefs_ij.col(2) *= rDeltaXYZ_ij.x()*rDeltaXYZ_ij.z();
    d2Coefs_ij.col(3) *= sqr(rDeltaXYZ_ij.y());
    d2Coefs_ij.col(4) *= rDeltaXYZ_ij.y()*rDeltaXYZ_ij.z();
    d2Coefs_ij.col(5) *= sqr(rDeltaXYZ_ij.z());
    return d2Coefs_ij;
}

Mat5X3 Foam::slope_Limiter
(
    const std::vector<vector>& normalList,
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ
)
{
    Mat5X3 d1_limited(d1_Own);
    Mat5X5 L, R;
    const scalar rho  = d0_Own(0);
    const vector rhoU = vector(d0_Own(1), d0_Own(2), d0_Own(3));
    const scalar rhoE = d0_Own(4);
    for (const auto& normal : normalList)
    {
        evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal);
        std::vector<Col5X1> d0_AdjListTemp;
        d0_AdjListTemp.reserve(d0_AdjList.size());
        for (const auto& d0_Adj : d0_AdjList)
        {
            d0_AdjListTemp.emplace_back(L*d0_Adj);
        }
        const Col5X1 d0_OwnTemp = L*d0_Own;
        const Mat5X3 d1_OwnTemp = L*d1_Own;
        Mat5X3 d1Temp = R*slope_Limiter(deltaList, d0_AdjListTemp, d0_OwnTemp, d1_OwnTemp, rDeltaXYZ);
        for (label i = 0; i < 5; i++)
        {
            for (label j = 0; j < 3; j++)
            {
                d1_limited(i,j) = MinMod(d1_limited(i,j), d1Temp(i,j));
            }
        }
    }
    return d1_limited;
}

Mat5X3 Foam::slope_Limiter
(
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ
)
{
    Mat5X3 d1_limited = d1_Own;
    Col5X1 pj_Min(Col5X1::Zero());
    Col5X1 pj_Max(Col5X1::Zero());
    Col5X1 d0_Min(d0_Own);
    Col5X1 d0_Max(d0_Own);
    for (const auto& delta : deltaList)
    {
        const Col5X1 project = d1_Own*Col3X1(delta[0]*rDeltaXYZ[0],
                                             delta[1]*rDeltaXYZ[1],
                                             delta[2]*rDeltaXYZ[2]);
        pj_Min = pj_Min.array().min(project.array());
        pj_Max = pj_Max.array().max(project.array());
    }
    for (const auto& d0_Adj : d0_AdjList)
    {
        d0_Min = d0_Min.array().min(d0_Adj.array());
        d0_Max = d0_Max.array().max(d0_Adj.array());
    }
    d0_Min -= d0_Own;
    d0_Max -= d0_Own;
    for (label i = 0; i < 5; i++)
    {
        scalar limiter = 1.0;
        limiter = min(limiter, Barth(pj_Min(i), d0_Min(i), d0_Max(i)));
        limiter = min(limiter, Barth(pj_Max(i), d0_Min(i), d0_Max(i)));
        d1_limited.row(i) *= limiter;
    }
    return d1_limited;
}

void Foam::New_WBAP_Limiter
(
    Mat5X9& coefs_i,
    Mat5X9& coefs_j,
    const vector& delta,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j,
    const scalar& value
)
{   
    std::vector<Mat5X3> d1tempList(2);
    std::vector<Mat5X6> d2tempList(2);
    // limit d2
    d2tempList[0] = coefs_i.block<5,6>(0,3);
    d2tempList[1] = transformD2Coefs(coefs_j.block<5,6>(0,3), rDeltaXYZ_i, rDeltaXYZ_j);
    Mat5X6 d2Temp = WBAP_L2_Limiter<6>(d2tempList, value);
    d2tempList[0] = coefs_j.block<5,6>(0,3);
    d2tempList[1] = transformD2Coefs(coefs_i.block<5,6>(0,3), rDeltaXYZ_j, rDeltaXYZ_i);
    coefs_i.block<5,6>(0,3) = d2Temp;
    coefs_j.block<5,6>(0,3) = WBAP_L2_Limiter<6>(d2tempList, value);
    // limit d1
    d1tempList[0] = coefs_i.block<5,3>(0,0);
    d1tempList[1] = transformD1Coefs(coefs_j, -delta, rDeltaXYZ_i, rDeltaXYZ_j);
    Mat5X3 d1Temp = WBAP_L2_Limiter<3>(d1tempList, value);
    d1tempList[0] = coefs_j.block<5,3>(0,0);
    d1tempList[1] = transformD1Coefs(coefs_i,  delta, rDeltaXYZ_j, rDeltaXYZ_i);
    coefs_i.block<5,3>(0,0) = d1Temp;
    coefs_j.block<5,3>(0,0) = WBAP_L2_Limiter<3>(d1tempList, value);    
}