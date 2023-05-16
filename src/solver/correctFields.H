void Foam::solver::correctFields()
{
    // Internal fields
    U_.ref() = rhoU_/rho_;
    p_.ref() = (rhoE_-0.5*rho_*magSqr(U_))*(Gamma-1.0);
    const bool rhoBool = Foam::positiveCorrect(rho_);
    const bool pBool   = Foam::positiveCorrect(p_);
    if (rhoBool || pBool)
    {
        rhoU_.ref() = rho_*U_;
        rhoE_.ref() = p_/(Gamma-1.0)+0.5*rho_*magSqr(U_);
    }
    T_.primitiveFieldRef() = Gamma*p_.primitiveField()/rho_.primitiveField();
    c_ = sqrt(T_.primitiveField());
    // Boundary conditions
    correctPrimitiveBoundary(rho_, U_, p_, T_);
    // Update other BCs
    rhoU_.correctBoundaryConditions();
    rhoE_.correctBoundaryConditions();
    if (fluidProps_.simulationType != "Euler") correctViscousFields();
}

void Foam::solver::correctPrimitiveBoundary
(
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p,
    volScalarField& T
)
{
    rho.correctBoundaryConditions();
    U.correctBoundaryConditions();
    p.correctBoundaryConditions();
    forAll(mesh_.boundary(), patchI)
    {
        const word name = mesh_.boundary()[patchI].name();
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        fvPatchScalarField& rhoBound = rho.boundaryFieldRef()[patchI];
        fvPatchVectorField& UBound = U.boundaryFieldRef()[patchI];
        fvPatchScalarField& pBound = p.boundaryFieldRef()[patchI];
        if (name == "farField")
        {
            if (fluidProps_.Mach_inf < 1.0)
            {
                forAll(bfaceCells, faceI)
                {
                    const label i = bfaceCells[faceI];
                    const scalar qn = U[i]&normal[faceI];
                    const scalar rhoc = c_[i]*rho[i]/Gamma;
                    if (qn < 0.0)
                    {
                        const scalar pB = pBound[faceI];
                        pBound[faceI] = 0.5*(pB+p[i]-rhoc*(normal[faceI]&(UBound[faceI]-U[i])));
                        rhoBound[faceI] += (pBound[faceI]-pB)*Gamma/sqr(c_[i]);
                        UBound[faceI] += normal[faceI]*(pBound[faceI]-pB)/rhoc;
                    }
                    else
                    {
                        rhoBound[faceI] = rho[i]+(pBound[faceI]-p[i])*Gamma/sqr(c_[i]);
                        UBound[faceI] = U[i]+normal[faceI]*(p[i]-pBound[faceI])/rhoc;
                    }
                }
            }
            else
            {
                forAll(bfaceCells, faceI)
                {
                    const label i = bfaceCells[faceI];
                    const scalar qn = U[i]&normal[faceI];
                    if (qn > 0.0)
                    {
                        rhoBound[faceI] = rho[i];
                        UBound[faceI]   = U[i];
                        pBound[faceI]   = p[i];
                    }
                }
            }
        }
        if (name == "inlet")
        {
            const scalarField Ma = mag(UBound)/sqrt(pBound*Gamma/rhoBound);
            forAll(bfaceCells, faceI)
            {
                if (Ma[faceI] < 1.0)
                {
                    const label i = bfaceCells[faceI];
                    const scalar rhoc = rho[i]*c_[i]/Gamma;
                    const scalar pB = pBound[faceI];
                    pBound[faceI] = 0.5*(pB+p[i]-rhoc*(normal[faceI]&(UBound[faceI]-U[i])));
                    rhoBound[faceI] += (pBound[faceI]-pB)*Gamma/sqr(c_[i]);
                    UBound[faceI] += normal[faceI]*(pBound[faceI]-pB)/rhoc;
                }
            }
        }
        if (name == "outlet")
        {
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const scalar Ma = mag(U[i])/c_[i];
                const scalar rhoc = rho[i]*c_[i]/Gamma;
                if (Ma >= 1.0)
                {
                    pBound[faceI]   = p[i];
                }
                else
                {
                    rhoBound[faceI] = rho[i]+(pBound[faceI]-p[i])*Gamma/sqr(c_[i]);
                    UBound[faceI] = U[i]+normal[faceI]*(p[i]-pBound[faceI])/rhoc;
                    if ((UBound[faceI]&normal[faceI]) < 0.0)
                    {
                        const vector signU = vector(sign(U[i].x()), sign(U[i].y()), sign(U[i].z()));
                        UBound[faceI] = cmptMultiply(signU, max(cmptMag(UBound[faceI]), cmptMag(U[i])));
                    }
                }
            }
        }
    }
    T.boundaryFieldRef() = p.boundaryFieldRef()*Gamma/rho.boundaryFieldRef();
}


void Foam::solver::correctViscousFields()
{
    laminarViscosity_.ref() = pow(T_, 1.5)*(1.0+S_T_)/(T_+S_T_);
    laminarViscosity_.boundaryFieldRef() = pow(T_.boundaryFieldRef(), 1.5)*(1.0+S_T_)/(T_.boundaryFieldRef()+S_T_);
    if (fluidProps_.simulationType == "SATurb")
    {
        eddyViscosity_.ref() = rho_*nuTilde_;
        eddyViscosity_.boundaryFieldRef() = rho_.boundaryFieldRef()*nuTilde_.boundaryFieldRef();
        nuMax_ = max(4.0/3.0, Gamma)*(laminarViscosity_.primitiveFieldRef()/Pr_Lam +
             eddyViscosity_.primitiveFieldRef()/Pr_Turb)/rho_.primitiveField();
    }
    else
    {
        nuMax_ = (max(4.0/3.0, Gamma)/Pr_Lam)*laminarViscosity_.primitiveFieldRef()/rho_.primitiveField();
    }
    
}