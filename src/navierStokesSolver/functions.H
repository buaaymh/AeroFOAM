void Foam::navierStokesSolver::evaluateMatrixLDU
(
    const scalar& alpha,
    scalarField& L,
    scalarField& U,
    scalarField& D
) const
{
    scalarField nuMax = (max(4.0/3.0, Gamma)/Pr_Lam)*muLam_.primitiveField()/rho_.primitiveField();
    // Calculate L, U and D
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const scalar magSf_vol = Ma_Re_*sqr(magSf[faceI])/mag(delta_[faceI]&Sf[faceI]);
        L[faceI] = 0.5*alpha*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI] + magSf_vol*nuMax[i]);
        U[faceI] = 0.5*alpha*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI] + magSf_vol*nuMax[j]);
        D[i] += L[faceI];
        D[j] += U[faceI];
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& Sf = mesh_.Sf().boundaryField()[patchI];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        forAll(bfaceCells, faceI)
        {
            const label i = bfaceCells[faceI];
            D[i] += 0.5*alpha*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        }
    }
}