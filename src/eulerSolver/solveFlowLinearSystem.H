void Foam::eulerSolver::solveFlowLinearSystemByLUSGS
(
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE
)
{
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateMatrixLDU(1.0, L, U, D);
    // Preconditioning step
    precondition(1.0, L, U, D, resRho, resRhoU, resRhoE, dRho_, dRhoU_, dRhoE_);
    rho_.primitiveFieldRef()  += dRho_;
    rhoU_.primitiveFieldRef() += dRhoU_;
    rhoE_.primitiveFieldRef() += dRhoE_;
}

void Foam::eulerSolver::solveFlowPseudoTimeSystem
(
    const scalar& dt,
    const scalar& alpha,
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE
)
{
    // Calculate L, U and D
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_ + mesh_.V().field()/dt);
    evaluateMatrixLDU(alpha, L, U, D);
    // Preconditioning step
    precondition(alpha, L, U, D, resRho, resRhoU, resRhoE, dRho_, dRhoU_, dRhoE_);
    rho_.primitiveFieldRef()  += dRho_;
    rhoU_.primitiveFieldRef() += dRhoU_;
    rhoE_.primitiveFieldRef() += dRhoE_;
}

void Foam::eulerSolver::solveFlowLinearSystemByGMRES
(
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE
)
{
    // Read settings from dictionary
    const label nDirs  = mesh_.solution().subDict("SOLVER").lookupOrDefault<label>("searchDir", 10);
    const label innerIter = mesh_.solution().subDict("SOLVER").lookupOrDefault<label>("innerIter", 5);
    const scalar relTol = mesh_.solution().subDict("SOLVER").lookupOrDefault<scalar>("relTol", 0.1);
    // Calculate L, U and D
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateMatrixLDU(1.0, L, U, D);
    dRho_  = scalarField(mesh_.nCells(), 0);
    dRhoU_ = vectorField(mesh_.nCells(), vector::zero);
    dRhoE_ = scalarField(mesh_.nCells(), 0);

    scalarField resRho0(mesh_.nCells());
    vectorField resRhoU0(mesh_.nCells());
    scalarField resRhoE0(mesh_.nCells());
    evaluateFlowResForGMRES(resRho0, resRhoU0, resRhoE0);
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
    label nIter = 0;
    while (nIter < innerIter)
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
            nIter++;
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
            dRho_  += y[k] * z_Rho[k];
            dRhoU_ += y[k] * z_RhoU[k];
            dRhoE_ += y[k] * z_RhoE[k];
        }
        if (finalRes < relTol*initialRes || nIter >= innerIter)
        {
            Info << "# Newton innerIter       [-] = " << nIter << nl
                 << "# Flow   |res|/|res0|    [-] = " << finalRes/initialRes << nl
                 << "----------------------------------------" << endl;
            break;
        }
        // Re-calculate the residual
        matrixVectorProduct(resRho0, resRhoU0, resRhoE0, dRho_, dRhoU_, dRhoE_, w_Rho[0], w_RhoU[0], w_RhoE[0]);
        w_Rho[0]  = resRho  - w_Rho[0];
        w_RhoU[0] = resRhoU - w_RhoU[0];
        w_RhoE[0] = resRhoE - w_RhoE[0];
        beta = Foam::sqrt(gSumSqr(w_Rho[0]) + gSum(magSqr(w_RhoU[0])) + gSumSqr(w_RhoE[0]));
    }
    forAll(mesh_.cells(), cellI)
    {
        scalar factor  = 1.0;
        if (rho_[cellI] + dRho_[cellI] <= 1e-3)
        {
            factor = min(factor, (1e-3-rho_[cellI])/dRho_[cellI]);
        }
        scalar rhoNew  = rho_[cellI] + factor*dRho_[cellI];
        vector rhoUNew = rhoU_[cellI]+ factor*dRhoU_[cellI];
        scalar rhoENew = rhoE_[cellI]+ factor*dRhoE_[cellI];
        scalar pNew = (Gamma-1.0)*(rhoENew - 0.5*magSqr(rhoUNew)/rhoNew);
        if (pNew <= 1e-3)
        {
            factor = min(factor, solveForPressure(rho_[cellI], rhoU_[cellI], rhoE_[cellI],
                                                  dRho_[cellI], dRhoU_[cellI], dRhoE_[cellI], factor));
        }
        if (factor < 1.0)
        {
            dRho_[cellI]  *= factor;
            dRhoU_[cellI] *= factor;
            dRhoE_[cellI] *= factor;
            rhoNew  = rho_[cellI] + dRho_[cellI];
            rhoUNew = rhoU_[cellI]+ dRhoU_[cellI];
            rhoENew = rhoE_[cellI]+ dRhoE_[cellI];
        }
    }
    rho_.primitiveFieldRef()  += dRho_;
    rhoU_.primitiveFieldRef() += dRhoU_;
    rhoE_.primitiveFieldRef() += dRhoE_;
}