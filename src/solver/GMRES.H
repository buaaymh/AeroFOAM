void Foam::solver::matrixVectorProduct
(
    const scalarField& L,
    const scalarField& U,
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
    y_Rho  = v_Rho/localDtDv_;
    y_RhoU = v_RhoU/localDtDv_;
    y_RhoE = v_RhoE/localDtDv_;

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
        y_Rho[i]  += dRhoFlux;
        y_RhoU[i] += dRhoUFlux;
        y_RhoE[i] += dRhoEFlux;
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
        y_Rho[j]  += dRhoFlux;
        y_RhoU[j] += dRhoUFlux;
        y_RhoE[j] += dRhoEFlux;
    }
}

void Foam::solver::solveFlowLinearSystem
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    // Read settings from dictionary
    const label nDirs  = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("searchDir", 8);
    const label innerIter = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("innerIter", 20);
    const scalar relTol = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relTol", 0.1);
    // Calculate L, U and D
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateMatrixLDU(L, U, D);
    scalarField dRho(resRho*localDtDv_);
    vectorField dRhoU(resRhoU*localDtDv_);
    scalarField dRhoE(resRhoE*localDtDv_);
    // Initialize V
    List<scalarField> v_Rho(nDirs+1, scalarField(mesh_.nCells()));
    List<vectorField> v_RhoU(nDirs+1, vectorField(mesh_.nCells()));
    List<scalarField> v_RhoE(nDirs+1, scalarField(mesh_.nCells()));
    matrixVectorProduct(resRho, resRhoU, resRhoE, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
    v_Rho[0]  = resRho  - v_Rho[0];
    v_RhoU[0] = resRhoU - v_RhoU[0];
    v_RhoE[0] = resRhoE - v_RhoE[0];
    scalar initialRes = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
    scalar finalRes = initialRes;
    label solverIter = 0;

    for (label m = 0; m < innerIter; m++)
    {
        if (finalRes < relTol*initialRes) break;
        // Create the Hesenberg matrix
        scalarSquareMatrix H(nDirs, 0);
        // Create y and b for Hessenberg matrix
        scalarField yh(nDirs, 0);
        scalarField bh(nDirs + 1, 0);
        // Givens rotation vectors
        scalarField c(nDirs, 0);
        scalarField s(nDirs, 0);
        precondition(L, U, D, v_Rho[0], v_RhoU[0], v_RhoE[0]);
        scalar beta = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
        // Set initial rhs and bh[0] = beta
        bh = 0;
        bh[0] = beta;
        for (label i = 0; i < nDirs; i++)
        {
            v_Rho[i]  /= beta;
            v_RhoU[i] /= beta;
            v_RhoE[i] /= beta;
            matrixVectorProduct(resRho, resRhoU, resRhoE, v_Rho[i], v_RhoU[i], v_RhoE[i], v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
            precondition(L, U, D, v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
            for (label j = 0; j <= i; j++)
            {
                beta = gSumProd(v_Rho[i+1], v_Rho[j]) + gSumProd(v_RhoU[i+1], v_RhoU[j]) + gSumProd(v_RhoE[i+1], v_RhoE[j]);;
                H[j][i] = beta;
                v_Rho[i+1]  -= beta * v_Rho[j];
                v_RhoU[i+1] -= beta * v_RhoU[j];
                v_RhoE[i+1] -= beta * v_RhoE[j];
            }
            beta = Foam::sqrt(gSumSqr(v_Rho[i+1]) + gSum(magSqr(v_RhoU[i+1])) + gSumSqr(v_RhoE[i+1]));
            // Apply previous Givens rotations to new column of H.
            for (label j = 0; j < i; j++)
            {
                const scalar Hji = H[j][i];
                H[j][i] = c[j]*Hji - s[j]*H[j+1][i];
                H[j+1][i] = s[j]*Hji + c[j]*H[j + 1][i];
            }
            // Apply Givens rotation to current row.
            givensRotation(H[i][i], beta, c[i], s[i]);
            const scalar bhi = bh[i];
            bh[i] = c[i]*bhi - s[i]*bh[i+1];
            bh[i+1] = s[i]*bhi + c[i]*bh[i+1];
            H[i][i] = c[i]*H[i][i] - s[i]*beta;
        }
        // Back substitute to solve Hy = b
        for (label i = nDirs-1; i >= 0; i--)
        {
            scalar sum = bh[i];
            for (label j = i+1; j < nDirs; j++)
            {
                sum -= H[i][j]*yh[j];
            }
            yh[i] = sum/H[i][i];
        }
        // Update solution
        for (label i = 0; i < nDirs; i++)
        {
            dRho  += yh[i] * v_Rho[i];
            dRhoU += yh[i] * v_RhoU[i];
            dRhoE += yh[i] * v_RhoE[i];
        }
        // Re-calculate the residual
        matrixVectorProduct(resRho, resRhoU, resRhoE, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
        v_Rho[0]  = resRho  - v_Rho[0];
        v_RhoU[0] = resRhoU - v_RhoU[0];
        v_RhoE[0] = resRhoE - v_RhoE[0];
        finalRes = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
        solverIter++;
    }
    Info << "GMRES iteration: " << solverIter << " Residual: " << finalRes << endl;
    rho_.primitiveFieldRef()  += dRho;
    rhoU_.primitiveFieldRef() += dRhoU;
    rhoE_.primitiveFieldRef() += dRhoE;
}

void Foam::givensRotation
(
    const scalar& h,
    const scalar& beta,
    scalar& c,
    scalar& s
)
{
    if (beta == 0)
    {
        c = 1;
        s = 0;
    }
    else if (mag(beta) > mag(h))
    {
        scalar tau = -h/beta;
        s = 1.0/Foam::sqrt(1.0 + sqr(tau));
        c = s*tau;
    }
    else
    {
        scalar tau = -beta/h;
        c = 1.0/Foam::sqrt(1.0 + sqr(tau));
        s = c*tau;
    }
}