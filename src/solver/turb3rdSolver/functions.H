// scalar Foam::eddyViscosityFunc
// (
//     const scalar nuTilde,
//     const scalar laminarViscosity,
//     const scalar rho
// ) const
// {
//     if (nuTilde <= 0.0) return 0.0;
//     const scalar nu  = laminarViscosity/rho;
//     const scalar X3  = pow3(nuTilde/nu);
//     const scalar fv1 = X3/(X3+Cv1_3);
//     return rho*nuTilde*fv1;
// }

// scalar Foam::turb3rdSolver::effectiveNuFunc
// (
//     const scalar nu,
//     const scalar nuTilde
// ) const
// {
//     scalar nu_e;
//     if (nuTilde > 0.0)
//     {
//         nu_e = nu + nuTilde;
//     }
//     else
//     {
//         const scalar X3 = pow3(nuTilde/nu);
//         const scalar fn = (Cn1 + X3)/(Cn1 - X3);
//         nu_e = nu + fn*nuTilde;
//     }
//     return nu_e;
// }

// void Foam::turb3rdSolver::computeProjectedGradient
// (
//     const vector& normal,
//     const vector& edgeVec,
//     const scalar& var_i,
//     const scalar& var_j,
//     const vector& grad_i,
//     const vector& grad_j,
//     scalar& proj_vector_ij,
//     scalar& correctedMeanGrad
// ) const
// {
//     proj_vector_ij = (edgeVec&normal)/magSqr(edgeVec);
//     const vector meanGrad = 0.5*(grad_i+grad_j); 
//     const scalar edgeProj = meanGrad&edgeVec;
//     correctedMeanGrad = (meanGrad&normal) - (edgeProj-(var_j-var_i))*proj_vector_ij;
// }

// void Foam::turb3rdSolver::calculateDelta()
// {
//     delta_ = wallDist::New(mesh_).y();
// }

void Foam::turb3rdSolver::evaluateVarsAndGrads
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
) const
{
    Col5X1 varProj = coefs.transpose()*polynomialDn0(delta, rDeltaXYZ, basisMean);
    rho  = rho_0  + varProj(0);
    rhoU = rhoU_0 + vector(varProj(1), varProj(2), varProj(3));
    rhoE = rhoE_0 + varProj(4);
}

void Foam::turb3rdSolver::evaluateVarsAndGrads
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
) const
{

}