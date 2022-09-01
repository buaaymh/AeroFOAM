void Foam::solver::matrixVectorProduct
(
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
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
        y_Rho[j]  -= dRhoFlux;
        y_RhoU[j] -= dRhoUFlux;
        y_RhoE[j] -= dRhoEFlux;
        y_Rho[i]  += dRhoFlux;
        y_RhoU[i] += dRhoUFlux;
        y_RhoE[i] += dRhoEFlux;

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
