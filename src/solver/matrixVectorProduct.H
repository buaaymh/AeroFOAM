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
    const scalar eps = Foam::sqrt(1e-15)/Foam::sqrt(gSumSqr(v_Rho)+gSum(magSqr(v_RhoU))+gSumSqr(v_RhoE));
    volScalarField rho = rho_;
    volVectorField U   = U_;
    volScalarField p   = p_;
    rho.primitiveFieldRef() += v_Rho*eps;
    U.primitiveFieldRef() = (rhoU_.primitiveField()+v_RhoU*eps)/rho.primitiveField();
    p.primitiveFieldRef() = (rhoE_.primitiveField()+v_RhoE*eps-0.5*rho.primitiveField()*magSqr(U.primitiveField()))*(fluidProps_.gamma-1.0);
    volScalarField T = fluidProps_.gamma*p/rho;
    correctPrimitiveBoundary(rho, U, p, T);
    volVectorField rhoGrad = fvc::grad(rho);
    volTensorField UGrad   = fvc::grad(U);
    volVectorField TGrad   = fvc::grad(T);
    evaluateFlowRes(rho, U, p, T, rhoGrad, UGrad, TGrad, y_Rho, y_RhoU, y_RhoE);
    y_Rho  = v_Rho /localDtDv_ - (y_Rho  - resRho) /eps;
    y_RhoU = v_RhoU/localDtDv_ - (y_RhoU - resRhoU)/eps;
    y_RhoE = v_RhoE/localDtDv_ - (y_RhoE - resRhoE)/eps;
}

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
)
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    y_Rho  = D*v_Rho;
    y_RhoU = D*v_RhoU;
    y_RhoE = D*v_RhoE;
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];

        scalar rho1  = rho_[i]+v_Rho[i];
        vector rhoU1 = rhoU_[i]+v_RhoU[i];
        scalar rhoE1 = rhoE_[i]+v_RhoE[i];
        scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
        scalar phi0 = U_[i]&Sf[faceI];
        scalar phi1 = (rhoU1/rho1)&Sf[faceI];
        scalar dRhoFlux  = L[faceI]*v_Rho[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0);
        vector dRhoUFlux = L[faceI]*v_RhoU[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI]);
        scalar dRhoEFlux = L[faceI]*v_RhoE[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0);
        // y_Rho[i]  += dRhoFlux;
        // y_RhoU[i] += dRhoUFlux;
        // y_RhoE[i] += dRhoEFlux;
        y_Rho[j]  -= dRhoFlux;
        y_RhoU[j] -= dRhoUFlux;
        y_RhoE[j] -= dRhoEFlux;

        rho1  = rho_[j]+v_Rho[j];
        rhoU1 = rhoU_[j]+v_RhoU[j];
        rhoE1 = rhoE_[j]+v_RhoE[j];
        p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
        phi0 = U_[j]&Sf[faceI];
        phi1 = (rhoU1/rho1)&Sf[faceI];
        dRhoFlux  = U[faceI]*v_Rho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0);
        dRhoUFlux = U[faceI]*v_RhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI]);
        dRhoEFlux = U[faceI]*v_RhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0);
        y_Rho[i]  -= dRhoFlux;
        y_RhoU[i] -= dRhoUFlux;
        y_RhoE[i] -= dRhoEFlux;
        // y_Rho[j]  += dRhoFlux;
        // y_RhoU[j] += dRhoUFlux;
        // y_RhoE[j] += dRhoEFlux;
    }
}