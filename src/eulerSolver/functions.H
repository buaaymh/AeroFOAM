void Foam::eulerSolver::evaluateMatrixLDU
(
    const scalar& alpha,
    scalarField& L,
    scalarField& U,
    scalarField& D
) const
{
    // Calculate L, U and D
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        L[faceI] = 0.5*alpha*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        U[faceI] = 0.5*alpha*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI]);
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
            D[i] += 0.5*alpha*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        }
    }
}

void Foam::eulerSolver::precondition
(
    const scalar& alpha,
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE,
    scalarField& dRho,
    vectorField& dRhoU,
    scalarField& dRhoE
) const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const scalar h = 1e-8/Foam::sqrt(gSumSqr(resRho)+gSum(magSqr(resRhoU))+gSumSqr(resRhoE)+1e-16);
    // Forward Sweep
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
                    const scalar rho1  = rho_[i]+dRho[i]*h;
                    const vector rhoU1 = rhoU_[i]+dRhoU[i]*h;
                    const scalar rhoE1 = rhoE_[i]+dRhoE[i]*h;
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(Gamma-1.0);
                    const scalar phi0 = U_[i]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += L[faceI]*dRho[i]  + 0.5*alpha*(rho1*phi1       - rho_[i]*phi0)/h;
                    rRhoU += L[faceI]*dRhoU[i] + 0.5*alpha*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI])/h;
                    rRhoE += L[faceI]*dRhoE[i] + 0.5*alpha*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0)/h;
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
                    const scalar rho1  = rho_[j]+dRho[j]*h;
                    const vector rhoU1 = rhoU_[j]+dRhoU[j]*h;
                    const scalar rhoE1 = rhoE_[j]+dRhoE[j]*h;
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(Gamma-1.0);
                    const scalar phi0 = U_[j]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += U[faceI]*dRho[j]  - 0.5*alpha*(rho1*phi1       - rho_[j]*phi0)/h;
                    rRhoU += U[faceI]*dRhoU[j] - 0.5*alpha*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI])/h;
                    rRhoE += U[faceI]*dRhoE[j] - 0.5*alpha*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0)/h;
                }
            }
        }
        dRho[cellI]  += rRho  / D[cellI];
        dRhoU[cellI] += rRhoU / D[cellI];
        dRhoE[cellI] += rRhoE / D[cellI];
    }
}

void Foam::eulerSolver::matrixVectorProduct
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
    const scalar eps = 1e-8/Foam::sqrt(gSumSqr(v_Rho)+gSum(magSqr(v_RhoU))+gSumSqr(v_RhoE)+1e-16);
    volScalarField rho  = rho_;
    volVectorField rhoU = rhoU_;
    volScalarField rhoE = rhoE_;
    rho_.primitiveFieldRef()  += v_Rho*eps;
    rhoU_.primitiveFieldRef() += v_RhoU*eps;
    rhoE_.primitiveFieldRef() += v_RhoE*eps;
    correctFields();
    evaluateFlowResForGMRES(y_Rho, y_RhoU, y_RhoE);
    y_Rho  = v_Rho /localDtDv_ - (y_Rho  - resRho) /eps;
    y_RhoU = v_RhoU/localDtDv_ - (y_RhoU - resRhoU)/eps;
    y_RhoE = v_RhoE/localDtDv_ - (y_RhoE - resRhoE)/eps;
    rho_  = rho;
    rhoU_ = rhoU;
    rhoE_ = rhoE;
}