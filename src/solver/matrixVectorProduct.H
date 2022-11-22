void Foam::solver::matrixVectorProduct
(
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
    const scalarField& v_Rho,
    const vectorField& v_RhoU,
    const scalarField& v_RhoE,
    scalarField& y_Rho,
    vectorField& y_RhoU,
    scalarField& y_RhoE
) const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const scalar h = 1e-8;
    y_Rho  = D*v_Rho;
    y_RhoU = D*v_RhoU;
    y_RhoE = D*v_RhoE;
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];

        scalar rho1  = rho_[i]+v_Rho[i]*h;
        vector rhoU1 = rhoU_[i]+v_RhoU[i]*h;
        scalar rhoE1 = rhoE_[i]+v_RhoE[i]*h;
        scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
        scalar phi0 = U_[i]&Sf[faceI];
        scalar phi1 = (rhoU1/rho1)&Sf[faceI];
        y_Rho[j]  -= L[faceI]*v_Rho[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0)/h;
        y_RhoU[j] -= L[faceI]*v_RhoU[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI])/h;
        y_RhoE[j] -= L[faceI]*v_RhoE[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0)/h;

        rho1  = rho_[j]+v_Rho[j]*h;
        rhoU1 = rhoU_[j]+v_RhoU[j]*h;
        rhoE1 = rhoE_[j]+v_RhoE[j]*h;
        p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
        phi0 = U_[j]&Sf[faceI];
        phi1 = (rhoU1/rho1)&Sf[faceI];
        y_Rho[i]  -= U[faceI]*v_Rho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0)/h;
        y_RhoU[i] -= U[faceI]*v_RhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI])/h;
        y_RhoE[i] -= U[faceI]*v_RhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0)/h;
    }
}

void Foam::solver::matrixVectorProduct
(
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE,
    const scalarField& v_Rho,
    const vectorField& v_RhoU,
    const scalarField& v_RhoE,
    scalarField& y_Rho,
    vectorField& y_RhoU,
    scalarField& y_RhoE
)
{
    const scalar eps = 1e-8/Foam::sqrt(gSumSqr(v_Rho)+gSum(magSqr(v_RhoU))+gSumSqr(v_RhoE));
    volScalarField rho  = rho_;
    volVectorField rhoU = rhoU_;
    volScalarField rhoE = rhoE_;
    rho_.primitiveFieldRef()  += v_Rho*eps;
    rhoU_.primitiveFieldRef() += v_RhoU*eps;
    rhoE_.primitiveFieldRef() += v_RhoE*eps;
    correctFields();
    evaluateFlowRes(y_Rho, y_RhoU, y_RhoE);
    y_Rho  = v_Rho /localDtDv_ - (y_Rho  - resRho) /eps;
    y_RhoU = v_RhoU/localDtDv_ - (y_RhoU - resRhoU)/eps;
    y_RhoE = v_RhoE/localDtDv_ - (y_RhoE - resRhoE)/eps;
    rho_  = rho;
    rhoU_ = rhoU;
    rhoE_ = rhoE;
}