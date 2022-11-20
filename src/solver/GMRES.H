void Foam::solver::solveFlowLinearSystem
(
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE
)
{
    // Read settings from dictionary
    const label nDirs  = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("searchDir", 5);
    const label innerIter = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("innerIter", 10);
    const scalar relTol = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relTol", 0.1);
    const scalar omega = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relaxOmega", 1.2);
    // Calculate L, U and D
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateMatrixLDU(omega, L, U, D);
    scalarField dRho(mesh_.nCells(), 0);
    vectorField dRhoU(mesh_.nCells(), vector::zero);
    scalarField dRhoE(mesh_.nCells(), 0);
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
        RectangularMatrix<scalar> H(nDirs+1, nDirs, 0);
        scalarField g(nDirs+1, 0);
        scalarField sn(nDirs+1, 0);
        scalarField cs(nDirs+1, 0);
        scalarField y(nDirs, 0);
        g[0] = beta;
        label i;
        for (i = 0; i < nDirs; i++)
        {
            if (finalRes < relTol*initialRes) break;
            w_Rho[i]  /= beta;
            w_RhoU[i] /= beta;
            w_RhoE[i] /= beta;
            precondition(omega, L, U, D, w_Rho[i], w_RhoU[i], w_RhoE[i], z_Rho[i], z_RhoU[i], z_RhoE[i]);
            matrixVectorProduct(L, U, D, z_Rho[i], z_RhoU[i], z_RhoE[i], w_Rho[i+1], w_RhoU[i+1], w_RhoE[i+1]);
            // matrixVectorProduct(resRho, resRhoU, resRhoE, z_Rho[i], z_RhoU[i], z_RhoE[i], w_Rho[i+1], w_RhoU[i+1], w_RhoE[i+1]);
            /*---  Modified Gram-Schmidt orthogonalization ---*/
            for (label j = 0; j < i+1; j++)
            {
                H[j][i] = gSumProd(w_Rho[i+1], w_Rho[j]) + gSumProd(w_RhoU[i+1], w_RhoU[j]) + gSumProd(w_RhoE[i+1], w_RhoE[j]);
                w_Rho[i+1]  -= H[j][i] * w_Rho[j];
                w_RhoU[i+1] -= H[j][i] * w_RhoU[j];
                w_RhoE[i+1] -= H[j][i] * w_RhoE[j];
            }
            beta = Foam::sqrt(gSumSqr(w_Rho[i+1]) + gSum(magSqr(w_RhoU[i+1])) + gSumSqr(w_RhoE[i+1]));
            H[i+1][i] = beta;
            // Apply previous Givens rotations to new column of H.
            for (label j = 0; j < i; j++)
                ApplyGivens(sn[j], cs[j], H[j][i], H[j+1][i]);
            GenerateGivens(H[i][i], H[i+1][i], sn[i], cs[i]);
            ApplyGivens(sn[i], cs[i], g[i], g[i+1]);

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
        matrixVectorProduct(L, U, D, dRho, dRhoU, dRhoE, w_Rho[0], w_RhoU[0], w_RhoE[0]);
        // matrixVectorProduct(resRho, resRhoU, resRhoE, dRho, dRhoU, dRhoE, w_Rho[0], w_RhoU[0], w_RhoE[0]);
        w_Rho[0]  = resRho  - w_Rho[0];
        w_RhoU[0] = resRhoU - w_RhoU[0];
        w_RhoE[0] = resRhoE - w_RhoE[0];
        beta = Foam::sqrt(gSumSqr(w_Rho[0]) + gSum(magSqr(w_RhoU[0])) + gSumSqr(w_RhoE[0]));
    }
    rho_.primitiveFieldRef()  += dRho;
    rhoU_.primitiveFieldRef() += dRhoU;
    rhoE_.primitiveFieldRef() += dRhoE;
}

void Foam::ApplyGivens
(
    scalar s,
    scalar c,
    scalar& h1,
    scalar& h2
)
{
    scalar temp = c*h1 + s*h2;
    h2 = c*h2 - s*h1;
    h1 = temp;
}

void Foam::GenerateGivens
(
    scalar& dx,
    scalar& dy,
    scalar& s,
    scalar& c
)
{
    if ( (dx == 0.0) && (dy == 0.0) ) {
        c = 1.0;
        s = 0.0;
    }
    else if ( mag(dy) > mag(dx) ) {
        scalar tmp = dx/dy;
        dx = sqrt(1.0 + tmp*tmp);
        s = mag(1.0/dx) * (dy < 0.0 ? -1.0 : 1.0);
        c = tmp*s;
    }
    else if ( mag(dy) <= mag(dx) ) {
        scalar tmp = dy/dx;
        dy = sqrt(1.0 + tmp*tmp);
        c = mag(1.0/dy) * (dx < 0.0 ? -1.0 : 1.0);
        s = tmp*c;
    }
    else {
        dx = 0.0;
        dy = 0.0;
        c = 1.0;
        s = 0.0;
    }
    dx = mag(dx*dy);
    dy = 0.0;
}

Foam::Mat5X5 Foam::evaluateAcJacobian
(
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& normal,
    const scalar& gamma
)
{
    const vector U = rhoU/rho;
    const scalar V = U&normal;
    const scalar a2 = gamma-1;
    const scalar phi = 0.5*a2*magSqr(U);
    const scalar a1 = gamma*rhoE/rho - phi;
    Mat5X5 jacobian = Mat5X5::Zero();
    for (label iDim = 0; iDim != 3; ++iDim)
    {
        jacobian(0, iDim+1) = normal[iDim];
        jacobian(iDim+1, 0) = normal[iDim]*phi-U[iDim]*V;
        jacobian(iDim+1, iDim+1) += V;
        jacobian(4, iDim+1) += normal[iDim]*a1 - a2*U[iDim]*V;
        jacobian(iDim+1, 4) += a2*normal[iDim];
        for (label jDim = 0; jDim != 3; ++jDim)
        {
            jacobian(iDim+1, jDim+1) += normal[jDim]*U[iDim] - a2*normal[iDim]*U[jDim];
        }
    }
    jacobian(4, 0) = V*(phi-a1);
    jacobian(4, 4) = gamma*V;
    return jacobian;
}