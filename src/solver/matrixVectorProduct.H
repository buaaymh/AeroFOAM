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