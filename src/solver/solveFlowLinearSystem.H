// void Foam::solver::solveFlowLinearSystem
// (
//     const scalarField& resRho,
//     const vectorField& resRhoU,
//     const scalarField& resRhoE
// )
// {
//     updateLTS();
//     scalarField L(mesh_.nInternalFaces());
//     scalarField U(mesh_.nInternalFaces());
//     scalarField D(1.0/localDtDv_);
//     evaluateMatrixLDU(L, U, D);
//     // Preconditioning step
//     scalarField dRho(resRho);
//     vectorField dRhoU(resRhoU);
//     scalarField dRhoE(resRhoE);
//     precondition(L, U, D, dRho, dRhoU, dRhoE);
//     rho_.primitiveFieldRef()  += dRho;
//     rhoU_.primitiveFieldRef() += dRhoU;
//     rhoE_.primitiveFieldRef() += dRhoE;
// }

void Foam::solver::solveFlowPseudoTimeSystem
(
    const scalar& dt,
    const scalar& alpha,
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE,
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
    scalarField dRho(resRho);
    vectorField dRhoU(resRhoU);
    scalarField dRhoE(resRhoE);
    precondition(L, U, D, dRho, dRhoU, dRhoE);
    resL1 = gSum(mag(resRho)*mesh_.V());
    rho_.primitiveFieldRef()  += dRho;
    rhoU_.primitiveFieldRef() += dRhoU;
    rhoE_.primitiveFieldRef() += dRhoE;
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
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& Sf = mesh_.Sf().boundaryField()[patchI];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        forAll(bfaceCells, faceI)
        {
            const label i = bfaceCells[faceI];
            D[i] += 0.5*omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        }
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
    scalarField dRho(resRho*localDtDv_);
    vectorField dRhoU(resRhoU*localDtDv_);
    scalarField dRhoE(resRhoE*localDtDv_);
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
                    const scalar rho1  = rho_[i]+dRho[i];
                    const vector rhoU1 = rhoU_[i]+dRhoU[i];
                    const scalar rhoE1 = rhoE_[i]+dRhoE[i];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[i]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += L[faceI]*dRho[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0);
                    rRhoU += L[faceI]*dRhoU[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI]);
                    rRhoE += L[faceI]*dRhoE[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0);
                }
            }
        }
        dRho[cellI]  = rRho  / D[cellI];
        dRhoU[cellI] = rRhoU / D[cellI];
        dRhoE[cellI] = rRhoE / D[cellI];
    }
    // Backward Sweep
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
                    const scalar rho1  = rho_[j]+resRho[j];
                    const vector rhoU1 = rhoU_[j]+resRhoU[j];
                    const scalar rhoE1 = rhoE_[j]+resRhoE[j];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[j]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += U[faceI]*resRho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0);
                    rRhoU += U[faceI]*resRhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI]);
                    rRhoE += U[faceI]*resRhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0);
                }
            }
        }
        resRho[cellI]  = dRho[cellI]  + rRho  / D[cellI];
        resRhoU[cellI] = dRhoU[cellI] + rRhoU / D[cellI];
        resRhoE[cellI] = dRhoE[cellI] + rRhoE / D[cellI];
    }
}