void Foam::turb3rdSolver::evaluateMatrixLDU
(
    scalarField& L,
    scalarField& U,
    scalarField& D
) const
{
    Info << "Solving Turb LUSGS." << endl;
    const scalar omega = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relaxOmega", 1.2);
    // Calculate L, U and D
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const scalar magSf_vol = 2.0*sqr(magSf[faceI])/mag((mesh_.C()[i]-mesh_.C()[j])&Sf[faceI]);
        L[faceI] = 0.5*(omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI])+magSf_vol*nuMax_[i]);
        U[faceI] = 0.5*(omega*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI])+magSf_vol*nuMax_[j]);
        D[i] += L[faceI];
        D[j] += U[faceI];
    }
    forAll(mesh_.boundary(), patchI)
    {
        const word type = mesh_.boundary()[patchI].type();
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& Sf = mesh_.boundary()[patchI].Sf();
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        if (type == "symmetryPlane" || type == "wall") continue;
        forAll(bfaceCells, faceI)
        {
            const label i = bfaceCells[faceI];
            D[i] += 0.5*omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        }
    }
}