void Foam::turbulenceSolver::evaluateTurbRes
(
    scalarField& resNuTilda
)
{
    resNuTilda  = scalarField(mesh_.nCells(), 0.0);
    nuTildaGrad_ = fvc::grad(nuTilda_);
    nuLam_ = muLam_/rho_;
    forAll(mesh_.owner(), faceI)
    {
        scalar nuTildaFlux = 0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const scalar q_ij = 0.5*((U_[i]+U_[j])&mesh_.Sf()[faceI]);
        const scalar a_i = 0.5*(q_ij+mag(q_ij));
        const scalar a_j = 0.5*(q_ij-mag(q_ij));
        nuTildaFlux += a_i*nuTilda_[i] + a_j*nuTilda_[j];
        const scalar nuLam_ij = 0.5*(nuLam_[i]+nuLam_[j]);
        const scalar nuTilda_ij = 0.5*(nuTilda_[i]+nuTilda_[j]);
        const scalar nu_e = nu(nuTilda_ij, nuLam_ij);
        nuTildaFlux -= Ma_Re_*nu_e*(0.5*(nuTildaGrad_[i]+nuTildaGrad_[j])&mesh_.Sf()[faceI])/SA::sigma;
        resNuTilda[i] -= nuTildaFlux;
        resNuTilda[j] += nuTildaFlux;
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& Sf = mesh_.Sf().boundaryField()[patchI];
        const word name = mesh_.boundary()[patchI].name();
        const word type = mesh_.boundary()[patchI].type();
        if (mesh_.boundary()[patchI].coupled())
        {
            const scalarField nuLam_neigh = nuLam_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField nuTilda_neigh = nuTilda_.boundaryField()[patchI].patchNeighbourField();
            const vectorField nuTildaGrad_neigh = nuTildaGrad_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar nuTildaFlux = 0;
                const label i = bfaceCells[faceI];
                const scalar q_ij = 0.5*((U_[i]+U_neigh[faceI])&Sf[faceI]);
                const scalar a_i = 0.5*(q_ij+mag(q_ij));
                const scalar a_j = 0.5*(q_ij-mag(q_ij));
                nuTildaFlux += a_i*nuTilda_[i] + a_j*nuTilda_neigh[faceI];
                const scalar nu_ij = 0.5*(nuLam_[i]+nuLam_neigh[faceI]);
                const scalar nuTilda_ij = 0.5*(nuTilda_[i]+nuTilda_neigh[faceI]);
                const scalar nu_e = nu(nuTilda_ij, nu_ij);
                nuTildaFlux -= Ma_Re_*nu_e*(0.5*(nuTildaGrad_[i]+nuTildaGrad_neigh[faceI])&Sf[faceI])/SA::sigma;
                resNuTilda[i] -= nuTildaFlux;

            }
        }
        if (type == "wall")
        {
            const fvPatchScalarField& nuLamBound = nuLam_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar nuTildaFlux = 0;
                const label i = bfaceCells[faceI];
                nuTildaFlux -= Ma_Re_*nuLamBound[faceI]*(nuTildaGrad_[i]&Sf[faceI])/SA::sigma;
                resNuTilda[i] -= nuTildaFlux;
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& nuLamBound = nuLam_.boundaryField()[patchI];
            const fvPatchScalarField& nuTildaBound = nuTilda_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar nuTildaFlux = 0;
                const label i = bfaceCells[faceI];
                nuTildaFlux += nuTildaBound[faceI]*(UBound[faceI]&Sf[faceI]);
                const scalar nu_e = nu(nuTildaBound[faceI], nuLamBound[faceI]);
                nuTildaFlux -= Ma_Re_*nu_e*(nuTildaGrad_[i]&Sf[faceI])/SA::sigma;
                resNuTilda[i] -= nuTildaFlux;
            }
        }
    }
    const volScalarField& y = wallDist::New(mesh_).y();
    forAll(mesh_.C(), cellI)
    {
        const scalar Omega = mag(vorticity(UGrad_[cellI]));
        resNuTilda[cellI] += evaluateSource(y[cellI], Omega, nuLam_[cellI], nuTilda_[cellI],
                                            nuTildaGrad_[cellI], mesh_.V()[cellI]);
    }
}