void Foam::solver::solveFlowPseudoTimeSystem
(
    const scalar& dt,
    const scalar& alpha,
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE,
    scalar& resL1
)
{
    // Calculate L, U and D
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_ + mesh_.V().field()/dt);
    evaluateMatrixLDU(L, U, D);
    // Preconditioning step
    precondition(L, U, D, resRho, resRhoU, resRhoE);
    resL1 = gSum(mag(resRho));
    rho_.primitiveFieldRef()  += resRho;
    rhoU_.primitiveFieldRef() += resRhoU;
    rhoE_.primitiveFieldRef() += resRhoE;
}