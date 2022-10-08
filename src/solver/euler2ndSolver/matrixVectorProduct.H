void Foam::euler2ndSolver::matrixVectorProduct
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
    volScalarField rho = rho_;
    volVectorField U   = U_;
    volScalarField p   = p_;
    rho.primitiveFieldRef() += v_Rho*eps;
    U.primitiveFieldRef() = (rhoU_.primitiveField()+v_RhoU*eps)/rho.primitiveField();
    p.primitiveFieldRef() = (rhoE_.primitiveField()+v_RhoE*eps-0.5*rho.primitiveField()*magSqr(U.primitiveField()))*(fluidProps_.gamma-1.0);
    correctPrimitiveBoundary(rho, U, p);
    volVectorField rhoGrad = fvc::grad(rho);
    volTensorField UGrad   = fvc::grad(U);
    volVectorField pGrad   = fvc::grad(p);
    evaluateFlowRes(rho, U, p, rhoGrad, UGrad, pGrad, y_Rho, y_RhoU, y_RhoE);
    y_Rho  = v_Rho /localDtDv_ - (y_Rho  - resRho) /eps;
    y_RhoU = v_RhoU/localDtDv_ - (y_RhoU - resRhoU)/eps;
    y_RhoE = v_RhoE/localDtDv_ - (y_RhoE - resRhoE)/eps;
}