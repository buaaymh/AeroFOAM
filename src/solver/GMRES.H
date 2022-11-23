void Foam::solver::solveFlowLinearSystem
(
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE
)
{
    // Read settings from dictionary
    const label nDirs  = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("searchDir", 10);
    const label innerIter = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("innerIter", 5);
    const scalar relTol = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relTol", 0.1);
    // Calculate L, U and D
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateMatrixLDU(1.0, L, U, D);
    scalarField dRho(mesh_.nCells(), 0);
    vectorField dRhoU(mesh_.nCells(), vector::zero);
    scalarField dRhoE(mesh_.nCells(), 0);

    scalarField resRho0(mesh_.nCells());
    vectorField resRhoU0(mesh_.nCells());
    scalarField resRhoE0(mesh_.nCells());
    solver::evaluateFlowRes(resRho0, resRhoU0, resRhoE0);
    // Initialize W
    List<scalarField> w_Rho(nDirs+1, scalarField(mesh_.nCells()));
    List<vectorField> w_RhoU(nDirs+1, vectorField(mesh_.nCells()));
    List<scalarField> w_RhoE(nDirs+1, scalarField(mesh_.nCells()));
    // Initialize Z
    List<scalarField> z_Rho(nDirs, scalarField(mesh_.nCells()));
    List<vectorField> z_RhoU(nDirs, vectorField(mesh_.nCells()));
    List<scalarField> z_RhoE(nDirs, scalarField(mesh_.nCells()));
    w_Rho[0]  = resRho;
    w_RhoU[0] = resRhoU;
    w_RhoE[0] = resRhoE;
    scalar beta = Foam::sqrt(gSumSqr(w_Rho[0]) + gSum(magSqr(w_RhoU[0])) + gSumSqr(w_RhoE[0]));
    scalar initialRes = beta;
    scalar finalRes = initialRes;
    label solverIter = 0;
    while (solverIter < innerIter)
    {
        scalarSquareMatrix H(nDirs, 0);
        scalarField g(nDirs+1, 0);
        scalarField sn(nDirs, 0);
        scalarField cs(nDirs, 0);
        scalarField y(nDirs, 0);
        g[0] = beta;
        label i;
        for (i = 0; i < nDirs; i++)
        {
            if (finalRes < relTol*initialRes) break;
            w_Rho[i]  /= beta;
            w_RhoU[i] /= beta;
            w_RhoE[i] /= beta;
            precondition(1.0, L, U, D, w_Rho[i], w_RhoU[i], w_RhoE[i], z_Rho[i], z_RhoU[i], z_RhoE[i]);
            matrixVectorProduct(resRho0, resRhoU0, resRhoE0, z_Rho[i], z_RhoU[i], z_RhoE[i], w_Rho[i+1], w_RhoU[i+1], w_RhoE[i+1]);
            /*---  Modified Gram-Schmidt orthogonalization ---*/
            for (label j = 0; j < i+1; j++)
            {
                H[j][i] = gSumProd(w_Rho[i+1], w_Rho[j]) + gSumProd(w_RhoU[i+1], w_RhoU[j]) + gSumProd(w_RhoE[i+1], w_RhoE[j]);
                w_Rho[i+1]  -= H[j][i] * w_Rho[j];
                w_RhoU[i+1] -= H[j][i] * w_RhoU[j];
                w_RhoE[i+1] -= H[j][i] * w_RhoE[j];
            }
            beta = Foam::sqrt(gSumSqr(w_Rho[i+1]) + gSum(magSqr(w_RhoU[i+1])) + gSumSqr(w_RhoE[i+1]));
            // Apply previous Givens rotations to new column of H.
            for (label j = 0; j < i; j++)
            {
                const scalar Hji = H[j][i];
                H[j][i]   = cs[j]*Hji - sn[j]*H[j+1][i];
                H[j+1][i] = sn[j]*Hji + cs[j]*H[j+1][i];
            }
            // Apply Givens rotation to current row.
            givensRotation(H[i][i], beta, cs[i], sn[i]);
            const scalar gi = g[i];
            g[i]   = cs[i]*gi - sn[i]*g[i+1];
            g[i+1] = sn[i]*gi + cs[i]*g[i+1];
            H[i][i] = cs[i]*H[i][i] - sn[i]*beta;

            finalRes = mag(g[i+1]);
            solverIter++;
        }
        /*---  Solve the least-squares system  ---*/
        for (label k = 0; k < i; k++)
            y[k] = g[k];
        for (label k = i-1; k >= 0; k--)
        {
            y[k] /= H[k][k];
            for (label j = k-1; j >= 0; j--)
                y[j] -= H[j][k]*y[k];
        }
        /*---  Update solution ---*/
        for (label k = 0; k < i; k++)
        {
            dRho  += y[k] * z_Rho[k];
            dRhoU += y[k] * z_RhoU[k];
            dRhoE += y[k] * z_RhoE[k];
        }
        if (finalRes < relTol*initialRes || solverIter >= innerIter)
        {
            Info << "# FGMRES final residual:"  << nl
                 << "# Iteration = " << solverIter << ": |res|/|res0| = " << finalRes/initialRes << "." << endl;
                 break;
        }
        // Re-calculate the residual
        matrixVectorProduct(resRho0, resRhoU0, resRhoE0, dRho, dRhoU, dRhoE, w_Rho[0], w_RhoU[0], w_RhoE[0]);
        w_Rho[0]  = resRho  - w_Rho[0];
        w_RhoU[0] = resRhoU - w_RhoU[0];
        w_RhoE[0] = resRhoE - w_RhoE[0];
        beta = Foam::sqrt(gSumSqr(w_Rho[0]) + gSum(magSqr(w_RhoU[0])) + gSumSqr(w_RhoE[0]));
    }
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