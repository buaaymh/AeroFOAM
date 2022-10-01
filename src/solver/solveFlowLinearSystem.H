// void Foam::solver::solveFlowLinearSystem
// (
//     scalarField& resRho,
//     vectorField& resRhoU,
//     scalarField& resRhoE
// )
// {
//     updateLTS();
//     scalarField L(mesh_.nInternalFaces());
//     scalarField U(mesh_.nInternalFaces());
//     scalarField D(1.0/localDtDv_);
//     evaluateMatrixLDU(L, U, D);
//     // Preconditioning step
//     precondition(L, U, D, resRho, resRhoU, resRhoE);
//     rho_.primitiveFieldRef()  += resRho;
//     rhoU_.primitiveFieldRef() += resRhoU;
//     rhoE_.primitiveFieldRef() += resRhoE;
// }

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
    resL1 = gSum(mag(resRho)*mesh_.V());
    rho_.primitiveFieldRef()  += resRho;
    rhoU_.primitiveFieldRef() += resRhoU;
    rhoE_.primitiveFieldRef() += resRhoE;
}

void Foam::solver::evaluateMatrixLDU
(
    scalarField& L,
    scalarField& U,
    scalarField& D
) const
{
    const scalar omega = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relaxOmega", 1.2);
    // Calculate L, U and D
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        L[faceI] = 0.5*omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        U[faceI] = 0.5*omega*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI]);
        D[i] += L[faceI];
        D[j] += U[faceI];
    }
}

void Foam::solver::precondition
(
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
) const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    // Forward Sweep
    scalarField dRhoStar(resRho*localDtDv_);
    vectorField dRhoUStar(resRhoU*localDtDv_);
    scalarField dRhoEStar(resRhoE*localDtDv_);
    forAll(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rRho = resRho[cellI];
        vector rRhoU = resRhoU[cellI];
        scalar rRhoE = resRhoE[cellI];
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label i = owner[faceI];
                if (i != cellI) /* L(cellI) < cellI */
                {
                    const scalar rho1  = rho_[i]+dRhoStar[i];
                    const vector rhoU1 = rhoU_[i]+dRhoUStar[i];
                    const scalar rhoE1 = rhoE_[i]+dRhoEStar[i];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[i]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += L[faceI]*dRhoStar[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0);
                    rRhoU += L[faceI]*dRhoUStar[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI]);
                    rRhoE += L[faceI]*dRhoEStar[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0);
                }
            }
        }
        dRhoStar[cellI]  = rRho  / D[cellI];
        dRhoUStar[cellI] = rRhoU / D[cellI];
        dRhoEStar[cellI] = rRhoE / D[cellI];
    }
    // Backward Sweep
    scalarField& dRho  = resRho;
    vectorField& dRhoU = resRhoU;
    scalarField& dRhoE = resRhoE;
    forAllReverse(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rRho  = 0;
        vector rRhoU = vector::zero;
        scalar rRhoE = 0;
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label j = neighbour[faceI];
                if (j != cellI) /* U(cellI) > cellI */
                {
                    const scalar rho1  = rho_[j]+dRho[j];
                    const vector rhoU1 = rhoU_[j]+dRhoU[j];
                    const scalar rhoE1 = rhoE_[j]+dRhoE[j];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[j]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += U[faceI]*dRho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0);
                    rRhoU += U[faceI]*dRhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI]);
                    rRhoE += U[faceI]*dRhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0);
                }
            }
        }
        dRho[cellI]  = dRhoStar[cellI]  + rRho  / D[cellI];
        dRhoU[cellI] = dRhoUStar[cellI] + rRhoU / D[cellI];
        dRhoE[cellI] = dRhoEStar[cellI] + rRhoE / D[cellI];
    }
}