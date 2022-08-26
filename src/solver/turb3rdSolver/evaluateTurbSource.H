void Foam::turb3rdSolver::evaluateTurbSource
(
    scalarField& resNuTilde,
    scalarField& D
) const
{
    forAll(mesh_.cells(), cellI)
    {
        scalar Production = 0.0;
        scalar Destruction = 0.0;
        scalar CrossProduction = 0.0;
        scalar Jacobian = 0.0;
        if (delta_[cellI] < 1e-10) continue;
        // evaluate Omega
        const vector vorticity(UGrad_[cellI].zy()-UGrad_[cellI].yz(),
                               UGrad_[cellI].xz()-UGrad_[cellI].zx(),
                               UGrad_[cellI].yx()-UGrad_[cellI].xy());
        const scalar Omega = mag(vorticity);
        if (nuTilde_[cellI] > 0.0)
        {
            const scalar delta_2 = Foam::sqr(delta_[cellI]);
            const scalar nu = laminarViscosity_[cellI]/rho_[cellI];
            const scalar X = nuTilde_[cellI]/nu;
            const scalar X_2 = Foam::sqr(X);
            const scalar X_3 = X*X_2;
            const scalar fv1 = X_3/(X_3+Cv1_3);
            const scalar fv2 = 1.0 - X/(1.0+X/fv1);
            // const scalar ft2 = Ct3*Foam::exp(-Ct4*X_2);
            const scalar S = Omega;
            const scalar inv_k2_d2 = 1.0/(k2*delta_2);
            scalar Shat = S + nuTilde_[cellI]*fv2*inv_k2_d2;
            Shat = max(Shat, 1.0e-10);
            const scalar inv_Shat = 1.0/Shat;

            /*--- Production term ---*/
            Production = Cb1*Shat*nuTilde_[cellI];
            /*--- Destruction term ---*/
            const scalar r = min(nuTilde_[cellI]*inv_Shat*inv_k2_d2, 10.0);
            const scalar g = r + Cw2*(pow(r,6.0)-r);
            const scalar g_6 =  pow(g,6.0);
            const scalar glim = pow((1.0+Cw3_6)/(g_6+Cw3_6), 1.0/6.0);
            const scalar fw = g*glim;
            Destruction = Cw1*fw*sqr(nuTilde_[cellI])/delta_2;
            /*--- Diffusion term ---*/
            const scalar norm2_Grad = magSqr(nuTildeGrad_[cellI]);
            CrossProduction = Cb2_sigma*norm2_Grad;

            resNuTilde += (Production-Destruction+CrossProduction)*mesh_.V()[cellI];

            /*--- Implicit part, production term ---*/
            const scalar dfv1 = 3.0*X_2*Cv1_3/(nu*sqr(X_3+Cv1_3));
            const scalar dfv2 = -(1.0/nu-X_2*dfv1)/sqr(1.0+X*fv1);
            scalar dShat;
            if ( Shat <= 1.0e-10 ) dShat = 0.0;
            else dShat = (fv2+nuTilde_[cellI]*dfv2)*inv_k2_d2;
            Jacobian += Cb1*(nuTilde_[cellI]*dShat+Shat);
            /*--- Implicit part, destruction term ---*/
            scalar dr = (Shat-nuTilde_[cellI]*dShat)*inv_Shat*inv_Shat*inv_k2_d2;
            if (r == 10.0) dr = 0.0;
            const scalar dg = dr*(1.0+Cw2*(6.0*pow(r,5.0)-1.0));
            const scalar dfw = dg*glim*(1.0-g_6/(g_6+Cw3_6));
            Jacobian -= Cw1*(dfw*nuTilde_[cellI] + 2.0*fw)*nuTilde_[cellI]/delta_2;
        }
        else
        {
            const scalar delta_2 = Foam::sqr(delta_[cellI]);
            /*--- Production term ---*/;
            Production = Cb1*(1.0-Ct3)*Omega*nuTilde_[cellI];
            /*--- Destruction term ---*/
            Destruction = Cw1*sqr(nuTilde_[cellI])/delta_2;
            /*--- Diffusion term ---*/
            const scalar norm2_Grad = magSqr(nuTildeGrad_[cellI]);
            CrossProduction = Cb2_sigma*norm2_Grad;

            resNuTilde += (Production+Destruction+CrossProduction)*mesh_.V()[cellI];

            /*--- Implicit part, production term ---*/
            Jacobian += Cb1*(1.0-Ct3)*Omega;
            /*--- Implicit part, destruction term ---*/
            Jacobian += 2.0*Cw1*nuTilde_[cellI]/delta_2;
        }
        D -= Jacobian*mesh_.V()[cellI];
    }
}