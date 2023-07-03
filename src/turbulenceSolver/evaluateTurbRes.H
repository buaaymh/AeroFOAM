void Foam::turbulenceSolver::evaluateTurbRes
(
    scalarField& resNuTilda
)
{
    resNuTilda  = scalarField(mesh_.nCells(), 0.0);
    nuTildaGrad_ = fvc::grad(nuTilda_);
    forAll(mesh_.owner(), faceI)
    {
        scalar nuTildaFlux = 0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        scalar dx_L  = mag((mesh_.Cf()[faceI] - mesh_.C()[i])&normal_[faceI]) + SMALL;
        scalar dx_R  = mag((mesh_.Cf()[faceI] - mesh_.C()[j])&normal_[faceI]) + SMALL;
        scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx;
        const scalar q_ij = (dx_L*U_[i]+dx_R*U_[j])&mesh_.Sf()[faceI];
        const scalar a_i = 0.5*(q_ij+mag(q_ij));
        const scalar a_j = 0.5*(q_ij-mag(q_ij));
        nuTildaFlux += a_i*nuTilda_[i] + a_j*nuTilda_[j];
        const scalar nuLam_ij = dx_L*nuLam_[i]+dx_R*nuLam_[j];
        const scalar nuTilda_ij = dx_L*nuTilda_[i]+dx_R*nuTilda_[j];
        const scalar nu_e = nu(nuTilda_ij, nuLam_ij);
        nuTildaFlux -= Ma_Re_*nu_e*((dx_L*nuTildaGrad_[i]+dx_R*nuTildaGrad_[j])&mesh_.Sf()[faceI])/SA::sigma;
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
            const vectorField delta_i = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            const vectorField delta_j = delta_i - mesh_.boundary()[patchI].delta();
            const scalarField nuLam_neigh = nuLam_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField nuTilda_neigh = nuTilda_.boundaryField()[patchI].patchNeighbourField();
            const vectorField nuTildaGrad_neigh = nuTildaGrad_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar nuTildaFlux = 0;
                const label i = bfaceCells[faceI];
                scalar dx_L  = mag(delta_i[faceI]&normal_[faceI]) + SMALL;
                scalar dx_R  = mag(delta_j[faceI]&normal_[faceI]) + SMALL;
                scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx;
                const scalar q_ij = (dx_L*U_[i]+dx_R*U_neigh[faceI])&Sf[faceI];
                const scalar a_i = 0.5*(q_ij+mag(q_ij));
                const scalar a_j = 0.5*(q_ij-mag(q_ij));
                nuTildaFlux += a_i*nuTilda_[i] + a_j*nuTilda_neigh[faceI];
                const scalar nu_ij = dx_L*nuLam_[i]+dx_R*nuLam_neigh[faceI];
                const scalar nuTilda_ij = dx_L*nuTilda_[i]+dx_R*nuTilda_neigh[faceI];
                const scalar nu_e = nu(nuTilda_ij, nu_ij);
                nuTildaFlux -= Ma_Re_*nu_e*((dx_L*nuTildaGrad_[i]+dx_R*nuTildaGrad_neigh[faceI])&Sf[faceI])/SA::sigma;
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
    forAll(mesh_.C(), cellI)
    {
        Omega_[cellI] = vorticity(UGrad_[cellI]);
        const scalar Omega = mag(Omega_[cellI]);
        resNuTilda[cellI] += evaluateSource(dist_[cellI], Omega, nuLam_[cellI], nuTilda_[cellI],
                                            nuTildaGrad_[cellI], mesh_.V()[cellI]);
    }
}