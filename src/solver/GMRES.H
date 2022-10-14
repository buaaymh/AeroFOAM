void Foam::solver::solveFlowLinearSystem
(
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE
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
    matrixVectorProduct(L, U, D, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
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
            matrixVectorProduct(L, U, D, v_Rho[i], v_RhoU[i], v_RhoE[i], v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
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
        matrixVectorProduct(L, U, D, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
        v_Rho[0]  = resRho  - v_Rho[0];
        v_RhoU[0] = resRhoU - v_RhoU[0];
        v_RhoE[0] = resRhoE - v_RhoE[0];
        finalRes = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
        solverIter++;
    }
    Info << "GMRES iterations=" << solverIter << ", initial residual=" << initialRes << ", final residual=" << finalRes << endl;
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