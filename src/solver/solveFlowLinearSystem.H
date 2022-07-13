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
        RectangularMatrix<scalar> H(nDirs+1, nDirs, Zero);
        scalarField g(nDirs+1, 0);
        scalarField sn(nDirs+1, 0);
        scalarField cs(nDirs+1, 0);
        scalarField y(nDirs, 0);

        precondition(L, U, D, v_Rho[0], v_RhoU[0], v_RhoE[0]);
        scalar beta = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
        v_Rho[0]  /= beta;
        v_RhoU[0] /= beta;
        v_RhoE[0] /= beta;

        g[0] = beta;
        label i = 0;
        for (i = 0; i < nDirs; i++)
        {
            if (g[i] < relTol*initialRes) break;
            matrixVectorProduct(L, U, D, v_Rho[i], v_RhoU[i], v_RhoE[i], v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
            precondition(L, U, D, v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
            for (label j = 0; j <= i; j++)
            {
                scalar prob = gSumProd(v_Rho[i+1], v_Rho[j]) + gSumProd(v_RhoU[i+1], v_RhoU[j]) + gSumProd(v_RhoE[i+1], v_RhoE[j]);
                H[j][i] = prob;
                v_Rho[i+1]  -= prob * v_Rho[j];
                v_RhoU[i+1] -= prob * v_RhoU[j];
                v_RhoE[i+1] -= prob * v_RhoE[j];
                
            }
            beta = Foam::sqrt(gSumSqr(v_Rho[i+1]) + gSum(magSqr(v_RhoU[i+1])) + gSumSqr(v_RhoE[i+1]));
            H[i+1][i] = beta;
            v_Rho[i+1]  /= beta;
            v_RhoU[i+1] /= beta;
            v_RhoE[i+1] /= beta;

            for (label j = 0; j < i; j++)
                ApplyGivens(sn[j], cs[j], H[j][i], H[j+1][i]);
            GenerateGivens(H[i][i], H[i+1][i], sn[i], cs[i]);
            ApplyGivens(sn[i], cs[i], g[i], g[i+1]);
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
            dRho  += y[k] * v_Rho[k];
            dRhoU += y[k] * v_RhoU[k];
            dRhoE += y[k] * v_RhoE[k];
        }
        // Re-calculate the residual
        matrixVectorProduct(L, U, D, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
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