void Foam::eulerConsVar3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    rhoGrad_  = fvc::grad(rho_);
    rhoUGrad_ = fvc::grad(rhoU_);
    rhoEGrad_ = fvc::grad(rhoE_);
    reconstructionIter();
    limitCoefficients();
}

void Foam::eulerConsVar3rdSolver::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& rhoU_L,
    const vector& rhoU_R,
    const scalar& rhoE_L,
    const scalar& rhoE_R,
    const vector& normal
) const
{
    scalar p_L, p_R, T_L, T_R;
    vector U_L, U_R;
    consToPrim(rho_L, rhoU_L, rhoE_L, fluidProps_.gamma, U_L, p_L, T_L);
    consToPrim(rho_R, rhoU_R, rhoE_R, fluidProps_.gamma, U_R, p_R, T_R);
    riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
}
