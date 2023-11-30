void Foam::turbulenceSolver::evaluateMatrixLDU
(
    const scalar& alpha,
    scalarField& L,
    scalarField& U,
    scalarField& D
) const
{
    scalarField nuMax = max(4.0/3.0, Gamma)*(muLam_.primitiveField()/Pr_Lam+muTurb_.primitiveField()/Pr_Turb)/rho_.primitiveField();
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
        const vectorField delta = mesh_.boundary()[patchI].delta();
        forAll(bfaceCells, faceI)
        {
            const label i = bfaceCells[faceI];
            const scalar magSf_vol = Ma_Re_*sqr(magSf[faceI])/mag(delta[faceI]&Sf[faceI]);
            D[i] += 0.5*alpha*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI] + magSf_vol*nuMax[i]);
        }
    }
}

void Foam::turbulenceSolver::evaluateTurbLDU
(
    const scalar& alpha,
    scalarField& L,
    scalarField& U,
    scalarField& D
) const
{
    // Calculate L, U and D
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        scalar dx_L = mag((mesh_.Cf()[faceI] - mesh_.C()[i])&normal_[faceI]) + SMALL;
        scalar dx_R = mag((mesh_.Cf()[faceI] - mesh_.C()[j])&normal_[faceI]) + SMALL;
        scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx;
        const scalar q_ij = (dx_L*U_[i]+dx_R*U_[j])&mesh_.Sf()[faceI];
        const scalar a_i = 0.5*(q_ij+mag(q_ij));
        const scalar a_j = 0.5*(q_ij-mag(q_ij));
        const scalar nuLam_ij = dx_L*nuLam_[i]+dx_R*nuLam_[j];
        const scalar nuTilda_ij = dx_L*nuTilda_[i]+dx_R*nuTilda_[j];
        const scalar nu_e = nu(nuTilda_ij, nuLam_ij);
        const scalar gradSf = 0.5*((dx_L*nuTildaGrad_[i]+dx_R*nuTildaGrad_[j])&Sf[faceI]);
        const scalar projSf = (delta_[faceI]&Sf[faceI])/magSqr(delta_[faceI]);
        L[faceI] = alpha*(a_i - Ma_Re_*(gradSf - nu_e*projSf)/SA::sigma);
        U[faceI] = alpha*(a_j - Ma_Re_*(gradSf + nu_e*projSf)/SA::sigma);
        D[i] += L[faceI];
        D[j] -= U[faceI];
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& Sf = mesh_.Sf().boundaryField()[patchI];
        const vectorField delta = mesh_.boundary()[patchI].delta();
        const word name = mesh_.boundary()[patchI].name();
        const word type = mesh_.boundary()[patchI].type();
        if (mesh_.boundary()[patchI].coupled())
        {
            const vectorField delta_i = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            const vectorField delta_j = delta_i - delta;
            const scalarField nuLam_neigh = nuLam_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField nuTilda_neigh = nuTilda_.boundaryField()[patchI].patchNeighbourField();
            const vectorField nuTildaGrad_neigh = nuTildaGrad_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                scalar dx_L  = mag(delta_i[faceI]&normal_[faceI]) + SMALL;
                scalar dx_R  = mag(delta_j[faceI]&normal_[faceI]) + SMALL;
                scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx;
                const scalar q_ij = (dx_L*U_[i]+dx_R*U_neigh[faceI])&Sf[faceI];
                const scalar a_i = 0.5*(q_ij+mag(q_ij));
                const scalar nu_ij = dx_L*nuLam_[i]+dx_R*nuLam_neigh[faceI];
                const scalar nuTilda_ij = dx_L*nuTilda_[i]+dx_R*nuTilda_neigh[faceI];
                const scalar nu_e = nu(nuTilda_ij, nu_ij);
                const scalar gradSf = 0.5*((dx_L*nuTildaGrad_[i]+dx_R*nuTildaGrad_neigh[faceI])&Sf[faceI]);
                const scalar projSf = (delta[faceI]&Sf[faceI])/magSqr(delta[faceI]);
                D[i] += alpha*(a_i - Ma_Re_*(gradSf - nu_e*projSf)/SA::sigma);
            }
        }
        if (type == "wall")
        {
            const fvPatchScalarField& nuLamBound = nuLam_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const scalar gradSf = 0.5*(nuTildaGrad_[i]&Sf[faceI]);
                const scalar projSf = (delta[faceI]&Sf[faceI])/magSqr(delta[faceI]);
                D[i] -= alpha*(Ma_Re_*(gradSf - nuLamBound[faceI]*projSf)/SA::sigma);
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& nuLamBound = rho_.boundaryField()[patchI];
            const fvPatchScalarField& nuTildaBound = nuTilda_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const scalar q_ij = UBound[faceI]&Sf[faceI];
                const scalar a_i = 0.5*(q_ij+mag(q_ij));
                const scalar nu_e = nu(nuTildaBound[faceI], nuLamBound[faceI]);
                const scalar gradSf = 0.5*(nuTildaGrad_[i]&Sf[faceI]);
                const scalar projSf = (delta[faceI]&Sf[faceI])/magSqr(delta[faceI]);
                D[i] += alpha*(a_i - Ma_Re_*(gradSf - nu_e*projSf)/SA::sigma);
            }
        }
    }
    forAll(mesh_.C(), cellI)
    {
        scalar source_Jacobian = evaluateSourceJacobian(dist_[cellI], Omega_[cellI], Strain_[cellI], nuLam_[cellI], nuTilda_[cellI], mesh_.V()[cellI]);
        D[cellI] -= 0.5*alpha*(source_Jacobian-mag(source_Jacobian));
    }
}

scalar Foam::turbulenceSolver::evaluateSource
(
    const scalar& d,
    const scalar& Omega,
    const scalar& Strain,
    const scalar& nuLam,
    const scalar& nuTilda,
    const vector& nuTildaGrad,
    const scalar& volume
) const
{
    scalar source = 0.0;
    if (d > 1e-10)
    {
        const scalar d2 = sqr(d);
        const scalar S  = 2.0*min(1e-10, Strain-Omega);
        if (nuTilda > 0.0)
        {
            /*--- Production term ---*/
            scalar X = nuTilda/nuLam;
            scalar fv2 = 1.0 - X/(1.0+X*fv1(X));
            scalar inv_k2_d2 = 1.0/(SA::k2*d2);
            scalar Shat = max(S + Ma_Re_*nuTilda*fv2*inv_k2_d2, 1e-10);
            scalar Production = SA::Cb1*Shat*nuTilda;
            /*--- Destruction term ---*/
            scalar r = min(Ma_Re_*nuTilda/Shat*inv_k2_d2, 10.0);
            scalar g = r + SA::Cw2*(pow(r,6.0)-r);
            scalar fw = g*pow((1.0+SA::Cw3_6)/(pow(g,6.0)+SA::Cw3_6),1.0/6.0);
            scalar Destruction = Ma_Re_*SA::Cw1*fw*sqr(nuTilda)/d2;
            /*--- Diffusion term ---*/
            scalar CrossProduction = Ma_Re_*SA::Cb2*magSqr(nuTildaGrad)/SA::sigma;
            source = (Production - Destruction + CrossProduction) * volume;
        }
        else
        {
            /*--- Production term ---*/;
            scalar Production = SA::Cb1*(1.0-SA::Ct3)*S*nuTilda;
            /*--- Destruction term ---*/
            scalar Destruction = Ma_Re_*SA::Cw1*sqr(nuTilda)/d2;
            /*--- Diffusion term ---*/
            scalar CrossProduction = Ma_Re_*SA::Cb2*magSqr(nuTildaGrad)/SA::sigma;
            source = (Production + Destruction + CrossProduction) * volume;
        }
    }
    return source;
}

scalar Foam::turbulenceSolver::evaluateSourceJacobian
(
    const scalar& d,
    const scalar& Omega,
    const scalar& Strain,
    const scalar& nuLam,
    const scalar& nuTilda,
    const scalar& volume
) const
{
    scalar jacobian = 0.0;
    if (d > 1e-10)
    {
        const scalar d2 = sqr(d);
        const scalar S  = 2.0*min(1e-10, Strain-Omega);
        if (nuTilda > 0.0)
        {
            /*--- Implicit part, production term ---*/
            scalar X = nuTilda/nuLam;
            scalar fv2 = 1.0 - X/(1.0+X*fv1(X));
            scalar inv_k2_d2 = 1.0/(SA::k2*d2);
            scalar Shat = max(S + Ma_Re_*nuTilda*fv2*inv_k2_d2, 1e-10);
            scalar inv_Shat = 1.0/Shat;
            scalar dfv1 = 3.0*sqr(X)*SA::Cv1_3/(nuLam*sqr(pow3(X)+SA::Cv1_3));
            scalar dfv2 = -(1/nuLam-sqr(X)*dfv1)/sqr(1.+X*fv1(X));
            scalar dShat = 0.0;
            if (Shat > 1.0e-10) dShat = Ma_Re_*(fv2+nuTilda*dfv2)*inv_k2_d2;
            jacobian += SA::Cb1*(nuTilda*dShat+Shat) * volume;
            /*--- Implicit part, destruction term ---*/
            scalar r = min(Ma_Re_*nuTilda/Shat*inv_k2_d2, 10.0);
            scalar g = r + SA::Cw2*(pow(r,6.0)-r);
            scalar g_6 =  pow(g,6.0);
            scalar glim = pow((1.0+SA::Cw3_6)/(g_6+SA::Cw3_6),1.0/6.0);
            scalar fw = g*glim;
            scalar dr = Ma_Re_*(Shat-nuTilda*dShat)*sqr(inv_Shat)*inv_k2_d2;
            if (r == 10.0) dr = 0.0;
            scalar dg = dr*(1.+SA::Cw2*(6.0*pow(r,5.0)-1.0));
            scalar dfw = dg*glim*(1.-g_6/(g_6+SA::Cw3_6));
            jacobian -= Ma_Re_*SA::Cw1*(dfw*nuTilda +  2.0*fw)*nuTilda/d2 * volume;
        }
        else
        {
            /*--- Implicit part, production term ---*/
            jacobian += SA::Cb1*(1.0-SA::Ct3)*S * volume;
            /*--- Implicit part, destruction term ---*/
            jacobian += Ma_Re_*2.0*SA::Cw1*nuTilda/d2 * volume;
        }
    }
    return jacobian;
}