void Foam::turb3rdSolver::solveTurbLinearSystem
(
    scalarField& resNuTilde,
    const scalarField& L,
    const scalarField& U,
    const scalarField& D
)
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    // Forward Sweep
    scalarField dNuTildeStar(resNuTilde*localDtDv_);
    forAll(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rNuTilde = resNuTilde[cellI];
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label i = owner[faceI];
                if (i != cellI) /* L(cellI) < cellI */
                {
                    rNuTilde += L[faceI]*dNuTildeStar[i];
                }
            }
        }
        dNuTildeStar[cellI] = rNuTilde / D[cellI];
    }
    // Backward Sweep
    scalarField& dNuTilde = resNuTilde;
    forAllReverse(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rNuTilde = 0;
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label j = neighbour[faceI];
                if (j != cellI) /* U(cellI) > cellI */
                {
                    rNuTilde -= U[faceI]*dNuTilde[j];
                }
            }
        }
        dNuTilde[cellI] = dNuTildeStar[cellI] + rNuTilde/D[cellI];
    }
    nuTilde_.primitiveFieldRef() += dNuTilde;
}

void Foam::turb3rdSolver::solveTurbPseudoTimeSystem
(
    const scalar& dt,
    const scalar& alpha,
    scalarField& resNuTilde,
    const scalarField& L,
    const scalarField& U,
    const scalarField& D
)
{
    
}