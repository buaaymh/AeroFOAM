void Foam::turb3rdSolver::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& U_L,
    const vector& U_R,
    const scalar& p_L,
    const scalar& p_R,
    const vector& rhoGrad_L,
    const vector& rhoGrad_R,
    const tensor& UGrad_L,
    const tensor& UGrad_R,
    const vector& pGrad_L,
    const vector& pGrad_R,
    const scalar& nuTilde_L,
    const scalar& nuTilde_R,
    const scalar& dxTilde,
    const vector& normal
) const
{
    riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
    const scalar rho = 0.5*(rho_L+rho_R);
    const vector U   = 0.5*(U_L+U_R);
    const scalar p   = 0.5*(p_L+p_R);
    const scalar T   = fluidProps_.gamma*p/rho;
    const scalar nuTilde = 0.5*(nuTilde_L+nuTilde_R);
    const vector rhoGrad = 0.5*(rhoGrad_L + rhoGrad_R + normal*(rho_R-rho_L)/dxTilde);
    const tensor UGrad   = 0.5*(UGrad_L   + UGrad_R   + normal*(U_R-U_L)/dxTilde);
    const vector pGrad   = 0.5*(pGrad_L   + pGrad_R   + normal*(p_R-p_L)/dxTilde);
    const vector TGrad   = fluidProps_.gamma/sqr(rho) * (rho*pGrad - p*rhoGrad);
    const scalar laminarViscosity = pow(T, 1.5)*(1.0+ST_inf)/(T+ST_inf);
    const scalar eddyViscosity    = eddyViscosityFunc(nuTilde, laminarViscosity, rho);
    const scalar mu = laminarViscosity + eddyViscosity;
    const scalar lambda = -2.0/3.0*mu;
    const scalar kappa  = (laminarViscosity/Pr_Lam+eddyViscosity/Pr_Turb)/(fluidProps_.gamma-1.0);
    const vector StrainRateU
    (
        (UGrad.xx()+UGrad.xx())*normal.x() + (UGrad.xy()+UGrad.yx())*normal.y() + (UGrad.xz()+UGrad.zx())*normal.z(),
        (UGrad.xy()+UGrad.yx())*normal.x() + (UGrad.yy()+UGrad.yy())*normal.y() + (UGrad.yz()+UGrad.zy())*normal.z(),
        (UGrad.xz()+UGrad.zx())*normal.x() + (UGrad.zy()+UGrad.yz())*normal.y() + (UGrad.zz()+UGrad.zz())*normal.z()
    );
    const scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
    const vector rhoUViscousFlux = mu*StrainRateU+lambda*divU*normal;
    rhoUFlux -= Ma_Re_inf*rhoUViscousFlux;
    rhoEFlux -= Ma_Re_inf*((rhoUViscousFlux&U)+kappa*(TGrad&normal));
}

void Foam::turb3rdSolver::evaluateWallFlux
(
    vector& rhoUFlux,
    const scalar& rho,
    const scalar& p,
    const tensor& UGrad,
    const vector& normal
) const
{
    rhoUFlux = p*normal;
    const scalar T  = fluidProps_.gamma*p/rho;
    const scalar mu = pow(T, 1.5)*(1.0+ST_inf)/(T+ST_inf);
    const scalar lambda = -2.0/3.0*mu;
    const vector StrainRateU
    (
        (UGrad.xx()+UGrad.xx())*normal.x() + (UGrad.xy()+UGrad.yx())*normal.y() + (UGrad.xz()+UGrad.zx())*normal.z(),
        (UGrad.xy()+UGrad.yx())*normal.x() + (UGrad.yy()+UGrad.yy())*normal.y() + (UGrad.yz()+UGrad.zy())*normal.z(),
        (UGrad.xz()+UGrad.zx())*normal.x() + (UGrad.zy()+UGrad.yz())*normal.y() + (UGrad.zz()+UGrad.zz())*normal.z()
    );
    const scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
    rhoUFlux -= Ma_Re_inf*(mu*StrainRateU+lambda*divU*normal);
}

void Foam::turb3rdSolver::evaluateInletOutletFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& U_L,
    const vector& U_R,
    const scalar& p_L,
    const scalar& p_R,
    const scalar& nuTilde,
    const vector& rhoGrad,
    const tensor& UGrad,
    const vector& pGrad,
    const vector& normal
) const
{
    riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
    const scalar rho = 0.5*(rho_L+rho_R);
    const vector U   = 0.5*(U_L+U_R);
    const scalar p   = 0.5*(p_L+p_R);
    const scalar T   = fluidProps_.gamma*p/rho;
    const vector TGrad   = fluidProps_.gamma/sqr(rho) * (rho*pGrad - p*rhoGrad);
    
    const scalar laminarViscosity = pow(T, 1.5)*(1.0+ST_inf)/(T+ST_inf);
    const scalar eddyViscosity    = eddyViscosityFunc(nuTilde, laminarViscosity, rho);
    const scalar mu = laminarViscosity + eddyViscosity;
    const scalar lambda = -2.0/3.0*mu;
    const scalar kappa  = (laminarViscosity/Pr_Lam+eddyViscosity/Pr_Turb)/(fluidProps_.gamma-1.0);
    const vector StrainRateU
    (
        (UGrad.xx()+UGrad.xx())*normal.x() + (UGrad.xy()+UGrad.yx())*normal.y() + (UGrad.xz()+UGrad.zx())*normal.z(),
        (UGrad.xy()+UGrad.yx())*normal.x() + (UGrad.yy()+UGrad.yy())*normal.y() + (UGrad.yz()+UGrad.zy())*normal.z(),
        (UGrad.xz()+UGrad.zx())*normal.x() + (UGrad.zy()+UGrad.yz())*normal.y() + (UGrad.zz()+UGrad.zz())*normal.z()
    );
    const scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
    const vector rhoUViscousFlux = mu*StrainRateU+lambda*divU*normal;
    rhoUFlux -= Ma_Re_inf*rhoUViscousFlux;
    rhoEFlux -= Ma_Re_inf*((rhoUViscousFlux&U)+kappa*(TGrad&normal));
}