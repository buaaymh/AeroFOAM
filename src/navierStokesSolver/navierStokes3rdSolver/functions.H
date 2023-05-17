void Foam::navierStokes3rdSolver::evaluateFlux
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
    const vector& normal,
    const scalar& magSf
) const
{
    scalar rhoFluxTemp;
    vector rhoUFluxTemp;
    scalar rhoEFluxTemp;

    scalar p_L, p_R, T_L, T_R;
    vector U_L, U_R;
    consToPrim(rho_L, rhoU_L, rhoE_L, U_L, p_L, T_L);
    consToPrim(rho_R, rhoU_R, rhoE_R, U_R, p_R, T_R);
    p_L = max(p_L, SMALL);
    p_R = max(p_R, SMALL);
    riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal);
    rhoFlux  += magSf * rhoFluxTemp;
    rhoUFlux += magSf * rhoUFluxTemp;
    rhoEFlux += magSf * rhoEFluxTemp;
}

void Foam::navierStokes3rdSolver::evaluateFlux
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
    const vector& rhoGrad_L,
    const vector& rhoGrad_R,
    const tensor& rhoUGrad_L,
    const tensor& rhoUGrad_R,
    const vector& rhoEGrad_L,
    const vector& rhoEGrad_R,
    const scalar& deltaTilde,
    const vector& normal,
    const scalar& magSf
) const
{

}

void Foam::navierStokes3rdSolver::evaluateWallFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& rhoGrad,
    const tensor& rhoUGrad,
    const vector& rhoEGrad,
    const vector& normal,
    const scalar& magSf
) const
{

}

void Foam::navierStokes3rdSolver::evaluateInletOutletFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rhoBound,
    const vector& UBound,
    const scalar& pBound,
    const scalar& rhoInter,
    const vector& rhoUInter,
    const scalar& rhoEInter,
    const vector& rhoGrad,
    const tensor& rhoUGrad,
    const vector& rhoEGrad,
    const vector& normal,
    const scalar& magSf
) const
{

}