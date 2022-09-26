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
    std::vector<Jacobian> rD(mesh_.nCells(), Jacobian::Zero());
    evaluateMatrixLDUEigen(L, U, rD);
    scalarField dRho(resRho*localDtDv_);
    vectorField dRhoU(resRhoU*localDtDv_);
    scalarField dRhoE(resRhoE*localDtDv_);
    // Initialize V
    List<scalarField> v_Rho(nDirs+1, scalarField(mesh_.nCells()));
    List<vectorField> v_RhoU(nDirs+1, vectorField(mesh_.nCells()));
    List<scalarField> v_RhoE(nDirs+1, scalarField(mesh_.nCells()));
    matrixVectorProduct(L, U, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
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

        preconditionEigen(L, U, rD, v_Rho[0], v_RhoU[0], v_RhoE[0]);
        scalar beta = Foam::sqrt(gSumSqr(v_Rho[0]) + gSum(magSqr(v_RhoU[0])) + gSumSqr(v_RhoE[0]));
        v_Rho[0]  /= beta;
        v_RhoU[0] /= beta;
        v_RhoE[0] /= beta;

        g[0] = beta;
        label i = 0;
        for (i = 0; i < nDirs; i++)
        {
            if (g[i] < relTol*initialRes) break;
            matrixVectorProduct(L, U, v_Rho[i], v_RhoU[i], v_RhoE[i], v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
            preconditionEigen(L, U, rD, v_Rho[i+1], v_RhoU[i+1], v_RhoE[i+1]);
            scalar nrm = gSumSqr(v_Rho[i+1]) + gSum(magSqr(v_RhoU[i+1])) + gSumSqr(v_RhoE[i+1]);
            for (label j = 0; j <= i; j++)
            {
                scalar prob = gSumProd(v_Rho[i+1], v_Rho[j]) + gSumProd(v_RhoU[i+1], v_RhoU[j]) + gSumProd(v_RhoE[i+1], v_RhoE[j]);
                H[j][i] = prob;
                v_Rho[i+1]  -= prob * v_Rho[j];
                v_RhoU[i+1] -= prob * v_RhoU[j];
                v_RhoE[i+1] -= prob * v_RhoE[j];
                if (prob*prob > 0.98*nrm)
                {
                    prob = gSumProd(v_Rho[i+1], v_Rho[j]) + gSumProd(v_RhoU[i+1], v_RhoU[j]) + gSumProd(v_RhoE[i+1], v_RhoE[j]);
                    H[j][i] += prob;
                    v_Rho[i+1]  -= prob * v_Rho[j];
                    v_RhoU[i+1] -= prob * v_RhoU[j];
                    v_RhoE[i+1] -= prob * v_RhoE[j];
                }
                nrm -= sqr(H[j][i]);
                nrm = max(nrm, 0.0);
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
        matrixVectorProduct(L, U, dRho, dRhoU, dRhoE, v_Rho[0], v_RhoU[0], v_RhoE[0]);
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

void Foam::solver::evaluateMatrixLDUEigen
(
    scalarField& L,
    scalarField& U,
    std::vector<Jacobian>& rD
) const
{
    const scalar omega = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relaxOmega", 1.2);
    // Calculate L, U and D
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const scalar ggm1  = fluidProps_.gamma/(fluidProps_.gamma-1.0);
    scalarField lambdaDiag(1.0/localDtDv_);
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        L[faceI] = 0.5*omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        U[faceI] = 0.5*omega*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI]);
        // const scalar H_L = ggm1*p_[i]/rho_[i] + 0.5*magSqr(U_[i]);
        // const scalar H_R = ggm1*p_[j]/rho_[j] + 0.5*magSqr(U_[j]);
        // const scalar rho_A = sqrt(rho_[i]*rho_[j]);
        // const scalar dd    = rho_A/rho_[i];
        // const scalar dd1   = 1.0/(1.0+dd);
        // const vector U_A   = (U_[i]+dd*U_[j])*dd1;
        // const scalar H_A   = (H_L+dd*H_R)*dd1;
        // const scalar c_A   = sqrt((fluidProps_.gamma-1)*(H_A-0.5*magSqr(U_A)));
        // L[faceI] = 0.5*omega*(mag(U_A&Sf[faceI])+c_A*magSf[faceI]);
        // U[faceI] = L[faceI];
        lambdaDiag[i] += L[faceI];
        lambdaDiag[j] += U[faceI];
        rD[i] += GetInviscidProjJac(rho_[i], U_[i], rhoE_[i], normal_[faceI], 0.5*omega*magSf[faceI]);
        rD[j] -= GetInviscidProjJac(rho_[j], U_[j], rhoE_[j], normal_[faceI], 0.5*omega*magSf[faceI]);
    }
    forAll(mesh_.cells(), cellI)
    {
        rD[cellI] += lambdaDiag[cellI]*Jacobian::Identity();
        rD[cellI] = rD[cellI].inverse();
    }
}

Foam::solver::Jacobian Foam::solver::GetInviscidProjJac
(
    const scalar rho,
    const vector U,
    const scalar rhoE,
    const vector normal,
    const scalar scale
) const
{
    Jacobian jacobian = Jacobian::Zero();
    const scalar sqvel = magSqr(U);
    const scalar proj_vel = U&normal;
    const scalar phi = 0.5*(fluidProps_.gamma-1.0)*sqvel;
    const scalar a1 = fluidProps_.gamma*rhoE/rho - phi;
    const scalar a2 = fluidProps_.gamma-1.0;
    for (label iDim = 0; iDim < 3; iDim++)
    {
        jacobian(0, iDim+1) = normal[iDim];
        jacobian(iDim+1, 0) = normal[iDim]*phi - U[iDim]*proj_vel;
        jacobian(iDim+1, iDim+1) = proj_vel;
        jacobian(iDim+1, 4) = a2*normal[iDim];
        jacobian(4, iDim+1) = a1*normal[iDim]-a2*U[iDim]*proj_vel;
        for (label jDim = 0; jDim < 3; jDim++)
            jacobian(iDim+1, jDim+1) += normal[jDim]*U[iDim]-a2*normal[iDim]*U[jDim];
    }
    jacobian(4, 0) = proj_vel*(phi-a1);
    jacobian(4, 4) = proj_vel*fluidProps_.gamma;
    jacobian *= scale;
    return jacobian;
}

void Foam::solver::preconditionEigen
(
    const scalarField& L,
    const scalarField& U,
    const std::vector<Jacobian>& rD,
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
) const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    // Forward Sweep
    scalarField dRhoStar(resRho*localDtDv_);
    vectorField dRhoUStar(resRhoU*localDtDv_);
    scalarField dRhoEStar(resRhoE*localDtDv_);
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
                    const scalar rho1  = rho_[i]+dRhoStar[i];
                    const vector rhoU1 = rhoU_[i]+dRhoUStar[i];
                    const scalar rhoE1 = rhoE_[i]+dRhoEStar[i];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[i]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += L[faceI]*dRhoStar[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0);
                    rRhoU += L[faceI]*dRhoUStar[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI]);
                    rRhoE += L[faceI]*dRhoEStar[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0);
                }
            }
        }
        VarColum temp = VarColum(rRho, rRhoU[0], rRhoU[1], rRhoU[2], rRhoE);
        temp = rD[cellI]*temp;
        dRhoStar[cellI]  = temp(0);
        dRhoUStar[cellI] = vector(temp(1), temp(2), temp(3));
        dRhoEStar[cellI] = temp(4);
    }
    // Backward Sweep
    scalarField& dRho  = resRho;
    vectorField& dRhoU = resRhoU;
    scalarField& dRhoE = resRhoE;
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
                    const scalar rho1  = rho_[j]+dRho[j];
                    const vector rhoU1 = rhoU_[j]+dRhoU[j];
                    const scalar rhoE1 = rhoE_[j]+dRhoE[j];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[j]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += U[faceI]*dRho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0);
                    rRhoU += U[faceI]*dRhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI]);
                    rRhoE += U[faceI]*dRhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0);
                }
            }
        }
        VarColum temp = VarColum(rRho, rRhoU[0], rRhoU[1], rRhoU[2], rRhoE);
        temp = rD[cellI]*temp;

        dRho[cellI]  = dRhoStar[cellI]  + temp(0);
        dRhoU[cellI] = dRhoUStar[cellI] + vector(temp(1), temp(2), temp(3));
        dRhoE[cellI] = dRhoEStar[cellI] + temp(4);
    }
}