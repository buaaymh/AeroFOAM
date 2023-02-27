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
    for (label i = 0; i != weights.size(); ++i)
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
        0.5*(sqr(B.x())+mag(basisMean.xx())), B.x()*B.y()+mag(basisMean.xy()), B.x()*B.z()+mag(basisMean.xz()),
        0.5*(sqr(B.y())+mag(basisMean.yy())), B.y()*B.z()+mag(basisMean.yz()),
        0.5*(sqr(B.z())+mag(basisMean.zz()))
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

Mat3X5 Foam::transformD1Coefs
(
    const Mat9X5& coefs,
    const vector& delta,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
)
{
    const vector rDeltaXYZ_ij = cmptDivide(rDeltaXYZ_j, rDeltaXYZ_i);
    const tensor gradTensor = rDeltaXYZ_ij * vector(cmptMultiply(delta, rDeltaXYZ_j));
    Mat3X5 d1Coefs_ij = coefs.block<3,5>(0,0);
    d1Coefs_ij.row(0) *= rDeltaXYZ_ij.x();
    d1Coefs_ij.row(1) *= rDeltaXYZ_ij.y();
    d1Coefs_ij.row(2) *= rDeltaXYZ_ij.z();
    d1Coefs_ij.col(0) += Col3X1(gradTensor.x()&vector(coefs(3,0), coefs(4,0), coefs(5,0)),
                                gradTensor.y()&vector(coefs(4,0), coefs(6,0), coefs(7,0)),
                                gradTensor.z()&vector(coefs(5,0), coefs(7,0), coefs(8,0)));
    d1Coefs_ij.col(1) += Col3X1(gradTensor.x()&vector(coefs(3,1), coefs(4,1), coefs(5,1)),
                                gradTensor.y()&vector(coefs(4,1), coefs(6,1), coefs(7,1)),
                                gradTensor.z()&vector(coefs(5,1), coefs(7,1), coefs(8,1)));
    d1Coefs_ij.col(2) += Col3X1(gradTensor.x()&vector(coefs(3,2), coefs(4,2), coefs(5,2)),
                                gradTensor.y()&vector(coefs(4,2), coefs(6,2), coefs(7,2)),
                                gradTensor.z()&vector(coefs(5,2), coefs(7,2), coefs(8,2)));
    d1Coefs_ij.col(3) += Col3X1(gradTensor.x()&vector(coefs(3,3), coefs(4,3), coefs(5,3)),
                                gradTensor.y()&vector(coefs(4,3), coefs(6,3), coefs(7,3)),
                                gradTensor.z()&vector(coefs(5,3), coefs(7,3), coefs(8,3)));
    d1Coefs_ij.col(4) += Col3X1(gradTensor.x()&vector(coefs(3,4), coefs(4,4), coefs(5,4)),
                                gradTensor.y()&vector(coefs(4,4), coefs(6,4), coefs(7,4)),
                                gradTensor.z()&vector(coefs(5,4), coefs(7,4), coefs(8,4)));
    return d1Coefs_ij;
}

Mat6X5 Foam::transformD2Coefs
(
    const Mat6X5& d2Coefs,
    const vector& rDeltaXYZ_i,
    const vector& rDeltaXYZ_j
)
{
    const vector rDeltaXYZ_ij = cmptDivide(rDeltaXYZ_j, rDeltaXYZ_i);
    Mat6X5 d2Coefs_ij = d2Coefs;
    d2Coefs_ij.row(0) *= sqr(rDeltaXYZ_ij.x());
    d2Coefs_ij.row(1) *= rDeltaXYZ_ij.x()*rDeltaXYZ_ij.y();
    d2Coefs_ij.row(2) *= rDeltaXYZ_ij.x()*rDeltaXYZ_ij.z();
    d2Coefs_ij.row(3) *= sqr(rDeltaXYZ_ij.y());
    d2Coefs_ij.row(4) *= rDeltaXYZ_ij.y()*rDeltaXYZ_ij.z();
    d2Coefs_ij.row(5) *= sqr(rDeltaXYZ_ij.z());
    return d2Coefs_ij;
}

void Foam::evaluateEigenMatrix
(
    Mat5X5& L,
    Mat5X5& R,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& normal,
    const scalar& gamma
)
{
    scalar rRho = 1.0/rho;
    vector U = rhoU*rRho;
    scalar ek = 0.5*magSqr(U);
    scalar p = (rhoE-rho*ek)*(gamma-1);
    scalar a = sqrt(p*gamma*rRho);
    scalar h = (rhoE+p)*rRho;
    scalar Vn = U&normal;
    // build Mat5x5 R
    vector an = a*normal;
    R.col(0) << 1, U.x()-an.x(), U.y()-an.y(), U.z()-an.z(), h-a*Vn;
    R.col(1) << 1, U.x()       , U.y()       , U.z()       , ek;
    R.col(2) << 1, U.x()+an.x(), U.y()+an.y(), U.z()+an.z(), h+a*Vn;
    // build Mat5x5 L
    scalar b1 = (gamma-1.0)/sqr(a), b2 = b1*ek, ra = 1.0/a;
    vector n_a = normal*ra, b1xU = b1*U;
    L.row(0) << (b2+Vn*ra), -(b1xU.x()+n_a.x()), -(b1xU.y()+n_a.y()), -(b1xU.z()+n_a.z()), b1;
    L.row(1) << -b2+1     ,   b1xU.x()         ,   b1xU.y()         ,   b1xU.z()     ,    -b1;
    L.row(2) << (b2-Vn*ra), -(b1xU.x()-n_a.x()), -(b1xU.y()-n_a.y()), -(b1xU.z()-n_a.z()), b1;
    L.row(0) *= 0.5;
    L.row(2) *= 0.5;
    if (mag(normal.x()) > 1e-8)
    {
        R.col(3) << 0,  normal.y(), -normal.x(), 0, U.x()*normal.y()-U.y()*normal.x();
        R.col(4) << 0, -normal.z(),  0, normal.x(), U.z()*normal.x()-U.x()*normal.z();
        scalar rnx = 1.0/normal.x();
        L.row(3) <<  (U.y()-Vn*normal.y())*rnx,  normal.y(),  (sqr(normal.y())-1.0)*rnx, normal.y()*normal.z()*rnx, 0;
        L.row(4) << -(U.z()-Vn*normal.z())*rnx, -normal.z(), -normal.y()*normal.z()*rnx, (1.0-sqr(normal.z()))*rnx, 0;
    }
    else if (mag(normal.y()) > 1e-8)
    {
        R.col(3) << 0, normal.y(), -normal.x(), 0, U.x()*normal.y()-U.y()*normal.x();
        R.col(4) << 0, 0, normal.z(), -normal.y(), U.y()*normal.z()-U.z()*normal.y();
        scalar rny = 1.0/normal.y();
        L.row(3) << -(U.x()-Vn*normal.x())*rny, (1.0-sqr(normal.x()))*rny, -normal.x(), normal.x()*normal.z()*rny, 0;
        L.row(4) <<  (U.z()-Vn*normal.z())*rny, normal.x()*normal.z()*rny,  normal.z(), (sqr(normal.z())-1.0)*rny, 0;
    }
    else
    {
        R.col(3) << 0, -normal.z(), 0, normal.x(), U.z()*normal.x()-U.x()*normal.z();
        R.col(4) << 0, 0, normal.z(), -normal.y(), U.y()*normal.z()-U.z()*normal.y();
        scalar rnz = 1.0/normal.z();
        L.row(3) <<  (U.x()-Vn*normal.x())*rnz,  (sqr(normal.x())-1.0)*rnz,  normal.x()*normal.y()*rnz,  normal.x(), 0;
        L.row(4) << -(U.y()-Vn*normal.y())*rnz, -normal.x()*normal.y()*rnz,  (1.0-sqr(normal.y()))*rnz, -normal.y(), 0;
    }
}

void Foam::euler3rdSolver::evaluateVars
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const Mat9X5& coefs,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
) const
{
    Col5X1 varProj = coefs.transpose()*polynomialDn0(delta, rDeltaXYZ, basisMean);
    rho  = rho_0  + varProj(0);
    rhoU = rhoU_0 + vector(varProj(1), varProj(2), varProj(3));
    rhoE = rhoE_0 + varProj(4);
}

void Foam::euler3rdSolver::evaluateVars
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const vector& d1Rho,
    const tensor& d1RhoU,
    const vector& d1RhoE,
    const symmTensor& d2Rho,
    const symmTensor& d2RhoUx,
    const symmTensor& d2RhoUy,
    const symmTensor& d2RhoUz,
    const symmTensor& d2RhoE,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
) const
{
    Col5X1 varProj = Mat5X9
    {
        {d1Rho[0],  d1Rho[1],  d1Rho[2],  d2Rho[0],   d2Rho[1],   d2Rho[2],   d2Rho[3],   d2Rho[4],   d2Rho[5]},
        {d1RhoU[0], d1RhoU[1], d1RhoU[2], d2RhoUx[0], d2RhoUx[1], d2RhoUx[2], d2RhoUx[3], d2RhoUx[4], d2RhoUx[5]},
        {d1RhoU[3], d1RhoU[4], d1RhoU[5], d2RhoUy[0], d2RhoUy[1], d2RhoUy[2], d2RhoUy[3], d2RhoUy[4], d2RhoUy[5]},
        {d1RhoU[6], d1RhoU[7], d1RhoU[8], d2RhoUz[0], d2RhoUz[1], d2RhoUz[2], d2RhoUz[3], d2RhoUz[4], d2RhoUz[5]},
        {d1RhoE[0], d1RhoE[1], d1RhoE[2], d2RhoE[0],  d2RhoE[1],  d2RhoE[2],  d2RhoE[3],  d2RhoE[4],  d2RhoE[5]}
    } * polynomialDn0(delta, rDeltaXYZ, basisMean);
    rho  = rho_0  + varProj(0);
    rhoU = rhoU_0 + vector(varProj(1), varProj(2), varProj(3));
    rhoE = rhoE_0 + varProj(4);
}

Mat5X3 Foam::slope1_Limiter
(
    const std::vector<vector>& normalList,
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ,
    const scalar& gamma
)
{
    Mat5X3 d1_limited(d1_Own);
    Mat5X5 L, R;
    const scalar rho  = d0_Own(0);
    const vector rhoU = vector(d0_Own(1), d0_Own(2), d0_Own(3));
    const scalar rhoE = d0_Own(4);
    for (const auto& normal : normalList)
    {
        evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal, gamma);
        std::vector<Col5X1> d0_AdjListTemp;
        d0_AdjListTemp.reserve(d0_AdjList.size());
        for (const auto& d0_Adj : d0_AdjList)
        {
            d0_AdjListTemp.emplace_back(L*d0_Adj);
        }
        const Col5X1 d0_OwnTemp = L*d0_Own;
        const Mat5X3 d1_OwnTemp = L*d1_Own;
        Mat5X3 d1Temp = R*slope1_Limiter(deltaList, d0_AdjListTemp, d0_OwnTemp, d1_OwnTemp, rDeltaXYZ);
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

Mat5X3 Foam::slope1_Limiter
(
    const std::vector<vector>& deltaList,
    const std::vector<Col5X1>& d0_AdjList,
    const Col5X1& d0_Own,
    const Mat5X3& d1_Own,
    const vector& rDeltaXYZ
)
{
    Mat5X3 d1_limited = d1_Own;
    Col5X1 pj_Min( GREAT*Col5X1::Ones());
    Col5X1 pj_Max(-GREAT*Col5X1::Ones());
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