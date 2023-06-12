void Foam::turbulenceSolver::solveTurbLinearSystemByLUSGS
(
    const scalarField& resNuTilda
)
{
    updateLTSforSATurb();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateTurbLDU(1.0, L, U, D);
    precondition(L, U, D, resNuTilda, dNuTilda_);
    nuTilda_.primitiveFieldRef() += dNuTilda_;
    correctTurbulenceFields();
}

void Foam::turbulenceSolver::solveTurbLinearSystemByGMRES
(
    const scalarField& resNuTilda
)
{
    // Read settings from dictionary
    const label nDirs  = mesh_.solution().subDict("SOLVER").lookupOrDefault<label>("searchDir", 10);
    const label innerIter = mesh_.solution().subDict("SOLVER").lookupOrDefault<label>("innerIter", 5);
    const scalar relTol = mesh_.solution().subDict("SOLVER").lookupOrDefault<scalar>("relTol", 0.1);
    // Calculate L, U and D
    updateLTSforSATurb();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateTurbLDU(1.0, L, U, D);
    dNuTilda_ = scalarField(mesh_.nCells(), 0);

    // Initialize W
    List<scalarField> w_NuTilda(nDirs+1, scalarField(mesh_.nCells()));
    // Initialize Z
    List<scalarField> z_NuTilda(nDirs, scalarField(mesh_.nCells()));
    w_NuTilda[0] = resNuTilda;
    scalar beta = Foam::sqrt(gSumSqr(w_NuTilda[0]));
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
            w_NuTilda[i]  /= beta;
            precondition(L, U, D, w_NuTilda[i], z_NuTilda[i]);
            matrixVectorProduct(L, U, D, z_NuTilda[i], w_NuTilda[i+1]);
            /*---  Modified Gram-Schmidt orthogonalization ---*/
            for (label j = 0; j < i+1; j++)
            {
                H[j][i] = gSumProd(w_NuTilda[i+1], w_NuTilda[j]);
                w_NuTilda[i+1]  -= H[j][i] * w_NuTilda[j];
            }
            beta = Foam::sqrt(gSumSqr(w_NuTilda[i+1]));
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
            dNuTilda_  += y[k] * z_NuTilda[k];
        }
        if (finalRes < relTol*initialRes || nIter >= innerIter)
        {
            Info << "# Newton innerIter       [-] = " << nIter << nl
                 << "# Turb   |res|/|res0|    [-] = " << finalRes/initialRes << nl
                 << "----------------------------------------" << endl;
            break;
        }
        // Re-calculate the residual
        matrixVectorProduct(L, U, D, dNuTilda_, w_NuTilda[0]);
        w_NuTilda[0]  = resNuTilda  - w_NuTilda[0];
        beta = Foam::sqrt(gSumSqr(w_NuTilda[0]));
    }
    nuTilda_.primitiveFieldRef() += dNuTilda_;
    correctTurbulenceFields();
}

void Foam::turbulenceSolver::precondition
(
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
    const scalarField& resNuTilda,
    scalarField& dNuTilda
) const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    // Forward Sweep /* L(cellI) < cellI */
    forAll(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rNuTilda = resNuTilda[cellI];
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label i = owner[faceI];
                if (i != cellI) rNuTilda += L[faceI] * dNuTilda[i]; 
            }
        }
        dNuTilda[cellI] = rNuTilda / D[cellI];
    }
    // Backward Sweep /* U(cellI) > cellI */
    forAllReverse(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rNuTilda = 0;
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label j = neighbour[faceI];
                if (j != cellI) rNuTilda -= U[faceI] *dNuTilda[j];
            }
        }
        dNuTilda[cellI] += rNuTilda / D[cellI];
    }
}

void Foam::turbulenceSolver::matrixVectorProduct
(
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
    const scalarField& v_NuTilda,
    scalarField& y_NuTilda
) const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    y_NuTilda  = D*v_NuTilda;
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        y_NuTilda[j] -= L[faceI]*v_NuTilda[i];
        y_NuTilda[i] += U[faceI]*v_NuTilda[j];
    }
}