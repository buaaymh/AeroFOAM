void Foam::turb3rdSolver::evaluateTurbRes
(
    scalarField& resNuTilde,
    scalarField& L,
    scalarField& U,
    scalarField& D
)
{
    resNuTilde = scalarField(mesh_.nCells(), 0.0);
    L = scalarField(mesh_.nInternalFaces(), 0.0);
    U = scalarField(mesh_.nInternalFaces(), 0.0);
    D = 1.0/localDtDv_;
    nuTildeGrad_ = fvc::grad(nuTilde_);
    calculateDelta();
    evaluateTurbSource(resNuTilde, D);
    forAll(mesh_.owner(), faceI)
    {
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        // convective
        const scalar q_ij = 0.5*(U_[i]+U_[j])&normal_[faceI];
        const scalar ac_i = 0.5*(q_ij+mag(q_ij));
        const scalar ac_j = 0.5*(q_ij-mag(q_ij));
        scalar nuTildeFlux = ac_i*nuTilde_[i] + ac_j*nuTilde_[j];
        scalar jacobian_i = ac_i;
        scalar jacobian_j = ac_j;
        // diffusion
        scalar proj_vector_ij;
        scalar correctedMeanGrad;
        computeProjectedGradient(normal_[faceI], mesh_.C()[j]-mesh_.C()[i],
                                 nuTilde_[i], nuTilde_[j],
                                 nuTildeGrad_[i], nuTildeGrad_[j],
                                 proj_vector_ij, correctedMeanGrad);
        const scalar nu_i = laminarViscosity_[i]/rho_[i];
        const scalar nu_j = laminarViscosity_[j]/rho_[j];
        const scalar nu_ij = 0.5 * (nu_i + nu_j);
        const scalar nuTilde_ij = 0.5 * (nuTilde_[i] + nuTilde_[j]);
        scalar nu_e;
        if (nuTilde_ij > 0.0) {
            nu_e = nu_ij + nuTilde_ij;
        }
        else {
            const scalar X3 = pow3(nuTilde_ij/nu_ij);
            const scalar fn = (Cn1 + X3)/(Cn1 - X3);
            nu_e = nu_ij + fn*nuTilde_ij;
        }
        nuTildeFlux -= nu_e*correctedMeanGrad/sigma;
        jacobian_i -= (0.5*correctedMeanGrad-nu_e*proj_vector_ij)/sigma;
        jacobian_j -= (0.5*correctedMeanGrad+nu_e*proj_vector_ij)/sigma;
        nuTildeFlux *= mesh_.magSf()[faceI];
        jacobian_i *= mesh_.magSf()[faceI];
        jacobian_j *= mesh_.magSf()[faceI];
        L[faceI] += jacobian_i;
        U[faceI] += jacobian_j;
        D[i] += jacobian_i;
        D[j] -= jacobian_j;
        resNuTilde[i] -= nuTildeFlux;
        resNuTilde[j] += nuTildeFlux;
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        const word name = mesh_.boundary()[patchI].name();
        const word type = mesh_.boundary()[patchI].type();
        if (isA<processorFvPatch>(mesh_.boundary()[patchI]))
        {
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const vectorField neighCn = mesh_.C().boundaryField()[patchI].patchNeighbourField();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh = T_.boundaryField()[patchI].patchNeighbourField();
            const scalarField nuTilde_neigh = nuTilde_.boundaryField()[patchI].patchNeighbourField();
            const vectorField nuTildeGrad_neigh = nuTildeGrad_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                // convective
                const scalar q_ij = 0.5*(U_[i]+U_neigh[j])&normal[j];
                const scalar ac_i = 0.5*(q_ij+mag(q_ij));
                const scalar ac_j = 0.5*(q_ij-mag(q_ij));
                scalar nuTildeFlux = ac_i*nuTilde_[i] + ac_j*nuTilde_neigh[j];
                // diffusion
                scalar proj_vector_ij;
                scalar correctedMeanGrad;
                computeProjectedGradient(normal[j], neighCn[j]-ownerCn[j],
                                         nuTilde_[i], nuTilde_neigh[j],
                                         nuTildeGrad_[i], nuTildeGrad_neigh[j],
                                         proj_vector_ij, correctedMeanGrad);
                const scalar nu_i = laminarViscosity_[i]/rho_[i];
                const scalar nu_j = pow(T_neigh[j], 2.0/3.0)*(1.0+ST_inf)/((T_neigh[j]+ST_inf)*rho_neigh[j]);
                const scalar nu_ij = 0.5 * (nu_i + nu_j);
                const scalar nuTilde_ij = 0.5 * (nuTilde_[i] + nuTilde_neigh[j]);
                scalar nu_e = effectiveNuFunc(nu_ij, nuTilde_ij);
                nuTildeFlux -= nu_e*correctedMeanGrad/sigma;
                nuTildeFlux *= magSf[j];
                D[i] += ac_i * magSf[j];
                resNuTilde[i] -= nuTildeFlux;
            }
        }
        if (type == "symmetryPlane") {}
        if (type == "wall")
        {
            forAll(bfaceCells, j)
            {  
                const label i = bfaceCells[j];
                // diffusion
                const scalar nu = laminarViscosity_[i]/rho_[i];
                scalar nu_e = effectiveNuFunc(nu, nuTilde_[i]);
                resNuTilde[i] -= nu_e*(nuTildeGrad_[i]&normal[j])*magSf[j]/sigma;
            }
        }
        if (name == "inlet")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& Ma = Ma_.boundaryField()[patchI];
            const fvPatchScalarField& TBound = T_.boundaryField()[patchI];
            const fvPatchScalarField& nuTildeBound = nuTilde_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                scalar nuTildeFlux = 0.0;
                scalar nu_ij;
                scalar nuTilde_ij;
                // convective
                if (Ma[faceI] < 1)
                {
                    const scalar q_ij = 0.5*(UBound[faceI]+U_[i])&normal[faceI];
                    const scalar ac_i = 0.5*(q_ij+mag(q_ij));
                    const scalar ac_j = 0.5*(q_ij-mag(q_ij));
                    nuTildeFlux += ac_i*nuTilde_[i] + ac_j*nuTildeBound[faceI];
                    D[i] += ac_i*magSf[faceI];
                    const scalar nu_i = laminarViscosity_[i]/rho_[i];
                    const scalar nu_j = pow(TBound[faceI], 2.0/3.0)*(1.0+ST_inf)/((TBound[faceI]+ST_inf)*rhoBound[faceI]);
                    nu_ij = 0.5 * (nu_i + nu_j);
                    nuTilde_ij = 0.5 * (nuTilde_[i] + nuTildeBound[faceI]);
                }
                else {
                    nuTildeFlux += nuTildeBound[faceI]*UBound[faceI]&normal[faceI];
                    nu_ij = pow(TBound[faceI], 2.0/3.0)*(1.0+ST_inf)/((TBound[faceI]+ST_inf)*TBound[faceI]);
                    nuTilde_ij = nuTildeBound[faceI];
                }
                // diffusion
                scalar nu_e = effectiveNuFunc(nu_ij, nuTilde_ij);
                nuTildeFlux -= nu_e*(nuTildeGrad_[i]&normal[faceI])/sigma;
                resNuTilde[i] -= nuTildeFlux*magSf[faceI];
            }
        }
        if (name == "outlet")
        {
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                scalar nuTildeFlux = 0.0;
                // convective
                const scalar ac = U_[i]&normal[faceI];
                nuTildeFlux += nuTilde_[i]*ac;
                // diffusion
                const scalar nu = laminarViscosity_[i]/rho_[i];
                scalar nu_e = effectiveNuFunc(nu, nuTilde_[i]);
                nuTildeFlux -= nu_e*(nuTildeGrad_[i]&normal[faceI])/sigma;
                D[i] += ac*magSf[faceI];
                resNuTilde[i] -= nuTildeFlux*magSf[faceI];
            }
        }
        if (name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& TBound = T_.boundaryField()[patchI];
            const fvPatchScalarField& Ma = Ma_.boundaryField()[patchI];
            const fvPatchScalarField& nuTildeBound = nuTilde_.boundaryField()[patchI];
            const scalarField Vn = normal&UBound;
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                scalar nuTildeFlux = 0.0;
                scalar nu_ij;
                scalar nuTilde_ij;
                if (Vn[faceI] > 0)
                {
                    const scalar ac = U_[i]&normal[faceI];
                    nuTildeFlux += ac*nuTilde_[i];
                    D[i] += ac*magSf[faceI];
                    nu_ij = laminarViscosity_[i]/rho_[i];
                    nuTilde_ij = nuTilde_[i];
                }
                else
                {
                    if (Ma[faceI] < 1.0)
                    {
                        const scalar q_ij = 0.5*(UBound[faceI]+U_[i])&normal[faceI];
                        const scalar ac_i = 0.5*(q_ij+mag(q_ij));
                        const scalar ac_j = 0.5*(q_ij-mag(q_ij));
                        nuTildeFlux += ac_i*nuTilde_[i] + ac_j*nuTildeBound[faceI];
                        D[i] += ac_i*magSf[faceI];
                        const scalar nu_i = laminarViscosity_[i]/rho_[i];
                        const scalar nu_j = pow(TBound[faceI], 2.0/3.0)*(1.0+ST_inf)/((TBound[faceI]+ST_inf)*rhoBound[faceI]);
                        nu_ij = 0.5 * (nu_i + nu_j);
                        nuTilde_ij = 0.5 * (nuTilde_[i] + nuTildeBound[faceI]);
                    }
                    else
                    {
                        nuTildeFlux += nuTildeBound[faceI]*UBound[faceI]&normal[faceI];
                        nu_ij = pow(TBound[faceI], 2.0/3.0)*(1.0+ST_inf)/((TBound[faceI]+ST_inf)*rhoBound[faceI]);
                        nuTilde_ij = nuTildeBound[faceI];
                    }
                    scalar nu_e = effectiveNuFunc(nu_ij, nuTilde_ij);
                    nuTildeFlux -= nu_e*(nuTildeGrad_[i]&normal[faceI])/sigma;
                    resNuTilde[i] -= nuTildeFlux*magSf[faceI];
                }
            }
        }
    }
}
