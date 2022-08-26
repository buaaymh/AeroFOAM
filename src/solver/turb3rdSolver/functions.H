scalar Foam::turb3rdSolver::eddyViscosityFunc
(
    const scalar nuTilde,
    const scalar laminarViscosity,
    const scalar rho
) const
{
    if (nuTilde <= 0.0) return 0.0;
    const scalar nu  = laminarViscosity/rho;
    const scalar X3  = pow3(nuTilde/nu);
    const scalar fv1 = X3/(X3+Cv1_3);
    return rho*nuTilde*fv1;
}

scalar Foam::turb3rdSolver::effectiveNuFunc
(
    const scalar nu,
    const scalar nuTilde
) const
{
    scalar nu_e;
    if (nuTilde > 0.0)
    {
        nu_e = nu + nuTilde;
    }
    else
    {
        const scalar X3 = pow3(nuTilde/nu);
        const scalar fn = (Cn1 + X3)/(Cn1 - X3);
        nu_e = nu + fn*nuTilde;
    }
    return nu_e;
}

void Foam::turb3rdSolver::computeProjectedGradient
(
    const vector& normal,
    const vector& edgeVec,
    const scalar& var_i,
    const scalar& var_j,
    const vector& grad_i,
    const vector& grad_j,
    scalar& proj_vector_ij,
    scalar& correctedMeanGrad
) const
{
    proj_vector_ij = (edgeVec&normal)/magSqr(edgeVec);
    const vector meanGrad = 0.5*(grad_i+grad_j); 
    const scalar edgeProj = meanGrad&edgeVec;
    correctedMeanGrad = (meanGrad&normal) - (edgeProj-(var_j-var_i))*proj_vector_ij;
}

void Foam::turb3rdSolver::calculateDelta()
{
    delta_ = wallDist::New(mesh_).y();
}