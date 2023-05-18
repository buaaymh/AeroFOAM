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
    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_L, rho_R, rhoU_L, rhoU_R, rhoE_L, rhoE_R, normal, magSf);
    scalar rho  = 0.5*(rho_L +rho_R);
    vector rhoU = 0.5*(rhoU_L+rhoU_R);
    scalar rhoE = 0.5*(rhoE_L+rhoE_R);
    vector U    = 0.5*(rhoU_L/rho_L + rhoU_R/rho_R);
    scalar T    = 0.5*Gamma*(Gamma-1.0)*((rhoE_L-0.5*magSqr(rhoU_L)/rho_L)/rho_L
                                       + (rhoE_R-0.5*magSqr(rhoU_R)/rho_R)/rho_R);
    vector rhoGrad  = 0.5*(rhoGrad_L + rhoGrad_R + normal*(rho_R -rho_L) /deltaTilde);
    tensor rhoUGrad = 0.5*(rhoUGrad_L+ rhoUGrad_R+ normal*(rhoU_R-rhoU_L)/deltaTilde);
    vector rhoEGrad = 0.5*(rhoEGrad_L+ rhoEGrad_R+ normal*(rhoE_R-rhoE_L)/deltaTilde);
    tensor UGrad = (rho*rhoUGrad-rhoU*rhoGrad)/sqr(rho);
    vector TGrad = Gamma*(Gamma-1)*(rhoEGrad-rhoE*rhoGrad/rho-(rhoU&UGrad))/rho;
    scalar laminarViscosity = pow(T, 1.5)*(1.0+S_T_)/(T+S_T_);
    scalar kappa  = (laminarViscosity/Pr_Lam)/(Gamma-1.0);
    vector StrainRateU = (UGrad + UGrad.T())&normal;
    scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
    vector rhoUViscousFlux = laminarViscosity*StrainRateU - 2.0/3.0*laminarViscosity*divU*normal;
    rhoUFlux -= magSf*Ma_Re_ * rhoUViscousFlux;
    rhoEFlux -= magSf*Ma_Re_ *((rhoUViscousFlux&U)+kappa*(TGrad&normal));
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
    vector rhoU_dummy = rhoU - 2*(rhoU&normal)*normal;
    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho, rho, rhoU, rhoU_dummy, rhoE, rhoE, normal, magSf);
    scalar p = (rhoE-0.5*magSqr(rhoU)/rho)*(Gamma-1.0);
    scalar T = Gamma*p/rho;
    tensor UGrad = (rho*rhoUGrad-rhoU*rhoGrad)/sqr(rho);
    scalar laminarViscosity = pow(T, 1.5)*(1.0+S_T_)/(T+S_T_);
    vector StrainRateU = (UGrad + UGrad.T())&normal;
    scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
    vector rhoUViscousFlux = laminarViscosity*StrainRateU - 2.0/3.0*laminarViscosity*divU*normal;
    rhoUFlux -= magSf*Ma_Re_ * rhoUViscousFlux;
}

void Foam::navierStokes3rdSolver::evaluateInletOutletFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rhoBound,
    const vector& UBound,
    const scalar& pBound,
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
    scalar rhoFluxTemp;
    vector rhoUFluxTemp;
    scalar rhoEFluxTemp;
    Foam::evaluateConvectiveFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                 rhoBound, UBound, pBound, normal);
    rhoFlux  += rhoFluxTemp*magSf;
    rhoUFlux += rhoUFluxTemp*magSf;
    rhoEFlux += rhoEFluxTemp*magSf;
    scalar T = Gamma*pBound/rhoBound;
    tensor UGrad = (rhoUGrad*rho-rhoU*rhoGrad)/sqr(rho);
    vector TGrad = Gamma*(Gamma-1)*(rhoEGrad-rhoE*rhoGrad/rho-(rhoU&UGrad))/rho;
    scalar laminarViscosity = pow(T, 1.5)*(1.0+S_T_)/(T+S_T_);
    scalar kappa  = (laminarViscosity/Pr_Lam)/(Gamma-1.0);
    vector StrainRateU = (UGrad + UGrad.T())&normal;
    scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
    vector rhoUViscousFlux = laminarViscosity*StrainRateU - 2.0/3.0*laminarViscosity*divU*normal;
    rhoUFlux -= magSf*Ma_Re_ * rhoUViscousFlux;
    rhoEFlux -= magSf*Ma_Re_ *((rhoUViscousFlux&UBound)+kappa*(TGrad&normal));
}