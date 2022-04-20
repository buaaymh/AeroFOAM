void Foam::eulerSolver::volProjectionsInit()
{
    volProjections_ = vectorField(mesh_.nCells(), vector::zero);
    // Get face-to-cell addressing: face area point from owner to neighbour
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    forAll(owner, faceI)
    {
        const vector faceVector = 0.5*cmptMag(Sf[faceI]);
        volProjections_[owner[faceI]] += faceVector;
        volProjections_[neighbour[faceI]] += faceVector;
    }
    forAll(mesh_.boundaryMesh(), patchI)
    {   
        scalar temp = 0.25;
        if (mesh_.boundary()[patchI].coupled()) temp = 0.5;
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const auto& Sf = mesh_.boundary()[patchI].Sf();
        forAll(bfaceCells, faceI)
            volProjections_[bfaceCells[faceI]] += temp*cmptMag(Sf[faceI]);
    }
}

void Foam::eulerSolver::updateLTS()
{
    scalar localCFL = 1.0;
    localDtDv_ = localCFL/(volProjections_&(cmptMag(U_.primitiveField())+c_.primitiveField()*vector::one));
}