void Foam::evaluateVars
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
)
{
    Col5X1 varProj = coefs.transpose()*polynomialDn0(delta, rDeltaXYZ, basisMean);
    rho  = rho_0  + varProj(0);
    rhoU = rhoU_0 + vector(varProj(1), varProj(2), varProj(3));
    rhoE = rhoE_0 + varProj(4);
}

void Foam::evaluateVars
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
)
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

void Foam::evaluateVarsAndGrads
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    vector& rhoGrad,
    tensor& rhoUGrad,
    vector& rhoEGrad,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const Mat9X5& coefs,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
)
{
    evaluateVars(rho, rhoU, rhoE, rho_0, rhoU_0, rhoE_0, coefs, basisMean, rDeltaXYZ, delta);
    Mat5X3 consGrad = coefs.transpose()*polynomialGrad(delta, rDeltaXYZ);
    rhoGrad  = vector(consGrad(0,0), consGrad(0,1), consGrad(0,2));
    rhoUGrad = tensor(consGrad(1,0), consGrad(1,1), consGrad(1,2),
                      consGrad(2,0), consGrad(2,1), consGrad(2,2),
                      consGrad(3,0), consGrad(3,1), consGrad(3,2));
    rhoEGrad = vector(consGrad(4,0), consGrad(4,1), consGrad(4,2));
}

void Foam::evaluateVarsAndGrads
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    vector& rhoGrad,
    tensor& rhoUGrad,
    vector& rhoEGrad,
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
)
{
    evaluateVars(rho, rhoU, rhoE, rho_0, rhoU_0, rhoE_0, d1Rho, d1RhoU, d1RhoE,
                 d2Rho, d2RhoUx, d2RhoUy, d2RhoUz, d2RhoE, basisMean, rDeltaXYZ, delta);
    Mat5X3 consGrad = Mat5X9
    {
        {d1Rho[0],  d1Rho[1],  d1Rho[2],  d2Rho[0],   d2Rho[1],   d2Rho[2],   d2Rho[3],   d2Rho[4],   d2Rho[5]},
        {d1RhoU[0], d1RhoU[1], d1RhoU[2], d2RhoUx[0], d2RhoUx[1], d2RhoUx[2], d2RhoUx[3], d2RhoUx[4], d2RhoUx[5]},
        {d1RhoU[3], d1RhoU[4], d1RhoU[5], d2RhoUy[0], d2RhoUy[1], d2RhoUy[2], d2RhoUy[3], d2RhoUy[4], d2RhoUy[5]},
        {d1RhoU[6], d1RhoU[7], d1RhoU[8], d2RhoUz[0], d2RhoUz[1], d2RhoUz[2], d2RhoUz[3], d2RhoUz[4], d2RhoUz[5]},
        {d1RhoE[0], d1RhoE[1], d1RhoE[2], d2RhoE[0],  d2RhoE[1],  d2RhoE[2],  d2RhoE[3],  d2RhoE[4],  d2RhoE[5]}
    } * polynomialGrad(delta, rDeltaXYZ);
    rhoGrad  = vector(consGrad(0,0), consGrad(0,1), consGrad(0,2));
    rhoUGrad = tensor(consGrad(1,0), consGrad(1,1), consGrad(1,2),
                      consGrad(2,0), consGrad(2,1), consGrad(2,2),
                      consGrad(3,0), consGrad(3,1), consGrad(3,2));
    rhoEGrad = vector(consGrad(4,0), consGrad(4,1), consGrad(4,2));
}