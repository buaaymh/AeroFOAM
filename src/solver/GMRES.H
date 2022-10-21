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
    label solverIter = 0;
    Info  << "GMRES initial residual = " << initialRes << endl;
    while (solverIter < innerIter)
    {
        RectangularMatrix<scalar> H(nDirs+1, nDirs, Zero);
        scalarField g(nDirs+1, 0);
        scalarField sn(nDirs+1, 0);
        scalarField cs(nDirs+1, 0);
        scalarField y(nDirs, 0);

        w_Rho[0]  /= beta;
        w_RhoU[0] /= beta;
        w_RhoE[0] /= beta;

        g[0] = beta;
        label i;
        for (i = 0; i < nDirs; i++)
        {
            if (beta < relTol*initialRes) break;
            precondition(L, U, D, w_Rho[i], w_RhoU[i], w_RhoE[i], z_Rho[i], z_RhoU[i], z_RhoE[i]);
            matrixVectorProduct(L, U, D, z_Rho[i], z_RhoU[i], z_RhoE[i], w_Rho[i+1], w_RhoU[i+1], w_RhoE[i+1]);
            /*---  Modified Gram-Schmidt orthogonalization ---*/
            for (label j = 0; j <= i; j++)
            {
                H[j][i] = gSumProd(w_Rho[i+1], w_Rho[j]) + gSumProd(w_RhoU[i+1], w_RhoU[j]) + gSumProd(w_RhoE[i+1], w_RhoE[j]);
                w_Rho[i+1]  -= H[j][i] * w_Rho[j];
                w_RhoU[i+1] -= H[j][i] * w_RhoU[j];
                w_RhoE[i+1] -= H[j][i] * w_RhoE[j];
            }
            H[i+1][i] = Foam::sqrt(gSumSqr(w_Rho[i+1]) + gSum(magSqr(w_RhoU[i+1])) + gSumSqr(w_RhoE[i+1]));
            w_Rho[i+1]  /= H[i+1][i];
            w_RhoU[i+1] /= H[i+1][i];
            w_RhoE[i+1] /= H[i+1][i];
            // Apply previous Givens rotations to new column of H.
            for (label j = 0; j < i; j++)
            {
                const scalar Hji = H[j][i];
                H[j][i]   = cs[j]*Hji - sn[j]*H[j+1][i];
                H[j+1][i] = sn[j]*Hji + cs[j]*H[j+1][i];
            }
            // Apply Givens rotation to current row.
            givensRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
            const scalar gi = g[i];
            g[i]   = cs[i]*gi - sn[i]*g[i+1];
            g[i+1] = sn[i]*gi + cs[i]*g[i+1];
            H[i][i] = cs[i]*H[i][i] - sn[i]*H[i+1][i];

            // for (label j = 0; j < i; j++)
            //     ApplyGivens(sn[j], cs[j], H[j][i], H[j+1][i]);
            // GenerateGivens(H[i][i], H[i+1][i], sn[i], cs[i]);
            // ApplyGivens(sn[i], cs[i], g[i], g[i+1]);
            
            /*---  Set L2 norm of residual and check if solution has converged ---*/
            beta = mag(g[i+1]);
            Info << "g[" << i+1 << "] = " << beta << endl;
        }
        solverIter += i;
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
        Info << "GMRES final residual = " << beta << endl;
        if (beta < relTol*initialRes) break;
        // Re-calculate the residual
        matrixVectorProduct(L, U, D, dRho, dRhoU, dRhoE, w_Rho[0], w_RhoU[0], w_RhoE[0]);
        w_Rho[0]  = resRho  - w_Rho[0];
        w_RhoU[0] = resRhoU - w_RhoU[0];
        w_RhoE[0] = resRhoE - w_RhoE[0];
        beta = Foam::sqrt(gSumSqr(w_Rho[0]) + gSum(magSqr(w_RhoU[0])) + gSumSqr(w_RhoE[0]));
    }
    rho_.primitiveFieldRef()  += dRho;
    rhoU_.primitiveFieldRef() += dRhoU;
    rhoE_.primitiveFieldRef() += dRhoE;
}

// void Foam::solver::solveFlowLinearSystem
// (
//     const scalarField& resRho,
//     const vectorField& resRhoU,
//     const scalarField& resRhoE
// )
// {
//     // Read settings from dictionary
//     const label nDirs  = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("searchDir", 8);
//     const label innerIter = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<label>("innerIter", 20);
//     const scalar relTol = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relTol", 0.1);
//     // Calculate L, U and D
//     updateLTS();
//     scalarField L(mesh_.nInternalFaces());
//     scalarField U(mesh_.nInternalFaces());
//     scalarField D(1.0/localDtDv_);
//     evaluateMatrixLDU(L, U, D);
//     scalarField dRho(mesh_.nCells(), 0);
//     vectorField dRhoU(mesh_.nCells(), vector::zero);
//     scalarField dRhoE(mesh_.nCells(), 0);
//     // Initialize V
//     List<scalarField> v_Rho(nDirs+1, scalarField(mesh_.nCells()));
//     List<vectorField> v_RhoU(nDirs+1, vectorField(mesh_.nCells()));
//     List<scalarField> v_RhoE(nDirs+1, scalarField(mesh_.nCells()));
//     v_Rho[0]  = resRho;
//     v_RhoU[0] = resRhoU;
//     v_RhoE[0] = resRhoE;
//     scalar initialRes = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
//     scalar finalRes = initialRes;
//     label solverIter = 0;

//     for (label m = 0; m < innerIter; m++)
//     {
//         if (finalRes < relTol*initialRes) break;
//         // Create the Hesenberg matrix
//         scalarSquareMatrix H(nDirs, 0);
//         // Create y and b for Hessenberg matrix
//         scalarField yh(nDirs, 0);
//         scalarField bh(nDirs + 1, 0);
//         // Givens rotation vectors
//         scalarField c(nDirs, 0);
//         scalarField s(nDirs, 0);
//         precondition(L, U, D, v_Rho[0], v_RhoU[0], v_RhoE[0]);
//         scalar beta = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
//         // Set initial rhs and bh[0] = beta
//         bh[0] = beta;
//         for (label i = 0; i < nDirs; i++)
//         {
//             v_Rho[i]  /= beta;
//             v_RhoU[i] /= beta;
//             v_RhoE[i] /= beta;
//             matrixVectorProduct(L, U, D, v_Rho[i], v_RhoU[i], v_RhoE[i], v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
//             precondition(L, U, D, v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
//             for (label j = 0; j <= i; j++)
//             {
//                 H[j][i] = gSumProd(v_Rho[i+1], v_Rho[j]) + gSumProd(v_RhoU[i+1], v_RhoU[j]) + gSumProd(v_RhoE[i+1], v_RhoE[j]);
//                 v_Rho[i+1]  -= H[j][i] * v_Rho[j];
//                 v_RhoU[i+1] -= H[j][i] * v_RhoU[j];
//                 v_RhoE[i+1] -= H[j][i] * v_RhoE[j];
//             }
//             beta = Foam::sqrt(gSumSqr(v_Rho[i+1]) + gSum(magSqr(v_RhoU[i+1])) + gSumSqr(v_RhoE[i+1]));
//             // Apply previous Givens rotations to new column of H.
//             for (label j = 0; j < i; j++)
//             {
//                 const scalar Hji = H[j][i];
//                 H[j][i] = c[j]*Hji - s[j]*H[j+1][i];
//                 H[j+1][i] = s[j]*Hji + c[j]*H[j + 1][i];
//             }
//             // Apply Givens rotation to current row.
//             givensRotation(H[i][i], beta, c[i], s[i]);
//             const scalar bhi = bh[i];
//             bh[i] = c[i]*bhi - s[i]*bh[i+1];
//             bh[i+1] = s[i]*bhi + c[i]*bh[i+1];
//             H[i][i] = c[i]*H[i][i] - s[i]*beta;
//         }
//         // Back substitute to solve Hy = b
//         for (label i = nDirs-1; i >= 0; i--)
//         {
//             scalar sum = bh[i];
//             for (label j = i+1; j < nDirs; j++)
//             {
//                 sum -= H[i][j]*yh[j];
//             }
//             yh[i] = sum/H[i][i];
//         }
//         // Update solution
//         for (label i = 0; i < nDirs; i++)
//         {
//             dRho  += yh[i] * v_Rho[i];
//             dRhoU += yh[i] * v_RhoU[i];
//             dRhoE += yh[i] * v_RhoE[i];
//         }
//         // Re-calculate the residual
//         matrixVectorProduct(L, U, D, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
//         v_Rho[0]  = resRho  - v_Rho[0];
//         v_RhoU[0] = resRhoU - v_RhoU[0];
//         v_RhoE[0] = resRhoE - v_RhoE[0];
//         finalRes = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
//         solverIter++;
//     }
//     Info << "GMRES iterations=" << solverIter << ", initial residual=" << initialRes << ", final residual=" << finalRes << endl;
//     rho_.primitiveFieldRef()  += dRho;
//     rhoU_.primitiveFieldRef() += dRhoU;
//     rhoE_.primitiveFieldRef() += dRhoE;
// }

void Foam::ApplyGivens
(
    const scalar s,
    const scalar c,
    scalar& h1,
    scalar& h2
)
{
    scalar temp = c*h1 + s*h2;
    h2 = c*h2 - s*h1;
    h1 = temp;
}

scalar Foam::Sign(scalar x, scalar y)
{
    if (y == 0.0) return 0.0;
    return mag(x) * (y < 0.0 ? -1.0 : 1.0);
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
        s = Sign(1.0/dx, dy);
        c = tmp*s;
    }
    else if ( mag(dy) <= mag(dx) ) {
        scalar tmp = dy/dx;
        dy = sqrt(1.0 + tmp*tmp);
        c = Sign(1.0/dy, dx);
        s = tmp*c;
    }
    else {
        // dx and/or dy must be invalid
        dx = 0.0;
        dy = 0.0;
        c = 1.0;
        s = 0.0;
    }
    dx = mag(dx*dy);
    dy = 0.0;
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