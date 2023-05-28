void Foam::turbulenceSolver::solveTurbLinearSystemByLUSGS
(
    const scalarField& resNuTilda
)
{
    updateLTS();
    scalarField L(mesh_.nInternalFaces());
    scalarField U(mesh_.nInternalFaces());
    scalarField D(1.0/localDtDv_);
    evaluateTurbLDU(1.0, L, U, D);
    
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    // Forward Sweep
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
                if (i != cellI) /* L(cellI) < cellI */
                {
                    rNuTilda += L[faceI] * dNuTilda_[i];
                }
            }
        }
        dNuTilda_[cellI] = rNuTilda / D[cellI];
    }
    // Backward Sweep
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
                if (j != cellI) /* U(cellI) > cellI */
                {
                    rNuTilda -= U[faceI] *dNuTilda_[j];
                }
            }
        }
        dNuTilda_[cellI] += rNuTilda / D[cellI];
    }
    nuTilda_.primitiveFieldRef() += dNuTilda_;
    correctTurbulenceFields();
}