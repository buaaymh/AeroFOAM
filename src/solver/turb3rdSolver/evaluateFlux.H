void Foam::turb3rdSolver::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& U_L,
    const vector& U_R,
    const scalar& T_L,
    const scalar& T_R,
    const tensor& UGrad_L,
    const tensor& UGrad_R,
    const vector& TGrad_L,
    const vector& TGrad_R,
    const scalar& nuTilde_L,
    const scalar& nuTilde_R,
    const scalar& dxTilde,
    const vector& normal
) const
{
    riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                           rho_L, rho_R, U_L, U_R, T_L, T_R,
                           normal, fluidProps_.gamma);
    const scalar two3 = 2.0/3.0;
    const tensor UGrad = 0.5 * (UGrad_L + UGrad_R + normal*(U_R-U_L)/dxTilde);
    const vector TGrad = 0.5 * (TGrad_L + TGrad_R + normal*(T_R-T_L)/dxTilde);
    const scalar rho_avg = 0.5 * (rho_L + rho_R);
    const vector U_avg = 0.5 * (U_L + U_R);
    const scalar T_avg = 0.5 * (T_L + T_R);
    const scalar nuTilde_avg = 0.5 * (nuTilde_L + nuTilde_R);
    const scalar laminarViscosity = pow(T_avg, 1.5)*(1.0+ST_inf)/(T_avg+ST_inf);
    const scalar eddyViscosity    = eddyViscosityFunc(nuTilde_avg, laminarViscosity, rho_avg);
    tensor tau;
    tau.xx() = two3*(2.0*UGrad.xx()-UGrad.yy()-UGrad.zz());
    tau.yy() = two3*(2.0*UGrad.yy()-UGrad.xx()-UGrad.zz());
    tau.zz() = two3*(2.0*UGrad.zz()-UGrad.xx()-UGrad.yy());
    tau.xy() = UGrad.xy() + UGrad.yx(); tau.yx() = tau.xy();
    tau.xz() = UGrad.xz() + UGrad.zx(); tau.zx() = tau.xz();
    tau.yz() = UGrad.yz() + UGrad.zy(); tau.zy() = tau.yz();
    rhoUFlux -= Ma_Re_inf*(laminarViscosity+eddyViscosity)*(tau&normal);
    rhoEFlux -= Ma_Re_inf*((laminarViscosity+eddyViscosity)*((tau&U_avg)&normal) +
                (laminarViscosity/Pr_Lam+eddyViscosity/Pr_Turb)/(fluidProps_.gamma-1.0)*(TGrad&normal));
}

void Foam::turb3rdSolver::evaluateWallFlux
(
    vector& rhoUFlux,
    const scalar& p,
    const scalar& T,
    const tensor& UGrad,
    const vector& normal
) const
{
    rhoUFlux = p*normal;
    const scalar two3 = 2.0/3.0;
    const scalar laminarViscosity = pow(T, 1.5)*(1.0+ST_inf)/(T+ST_inf);
    tensor tau;
    tau.xx() = two3*(2.0*UGrad.xx()-UGrad.yy()-UGrad.zz());
    tau.yy() = two3*(2.0*UGrad.yy()-UGrad.xx()-UGrad.zz());
    tau.zz() = two3*(2.0*UGrad.zz()-UGrad.xx()-UGrad.yy());
    tau.xy() = UGrad.xy() + UGrad.yx(); tau.yx() = tau.xy();
    tau.xz() = UGrad.xz() + UGrad.zx(); tau.zx() = tau.xz();
    tau.yz() = UGrad.yz() + UGrad.zy(); tau.zy() = tau.yz();
    rhoUFlux -= Ma_Re_inf*(laminarViscosity)*(tau&normal);
}

void Foam::turb3rdSolver::evaluateInletOutletFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho,
    const vector& U,
    const scalar& p,
    const scalar& nuTilde,
    const tensor& UGrad,
    const vector& TGrad,
    const vector& normal
) const
{
    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                       rho, U, p, normal, fluidProps_.gamma);
    const scalar two3 = 2.0/3.0;
    const scalar T = p * fluidProps_.gamma / rho;
    const scalar laminarViscosity = pow(T, 1.5)*(1.0+ST_inf)/(T+ST_inf);
    const scalar eddyViscosity    = eddyViscosityFunc(nuTilde, laminarViscosity, rho);
    tensor tau;
    tau.xx() = two3*(2.0*UGrad.xx()-UGrad.yy()-UGrad.zz());
    tau.yy() = two3*(2.0*UGrad.yy()-UGrad.xx()-UGrad.zz());
    tau.zz() = two3*(2.0*UGrad.zz()-UGrad.xx()-UGrad.yy());
    tau.xy() = UGrad.xy() + UGrad.yx(); tau.yx() = tau.xy();
    tau.xz() = UGrad.xz() + UGrad.zx(); tau.zx() = tau.xz();
    tau.yz() = UGrad.yz() + UGrad.zy(); tau.zy() = tau.yz();
    rhoUFlux -= Ma_Re_inf*(laminarViscosity+eddyViscosity)*(tau&normal);
    rhoEFlux -= Ma_Re_inf*((laminarViscosity+eddyViscosity)*((tau&U)&normal) +
                (laminarViscosity/Pr_Lam+eddyViscosity/Pr_Turb)/(fluidProps_.gamma-1.0)*(TGrad&normal));
}