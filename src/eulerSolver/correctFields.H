void Foam::eulerSolver::correctFields()
{
    conservativeToPrimitiveFields();
    c_ = sqrt(T_.primitiveField());
    // Update other BCs
    rhoU_.correctBoundaryConditions();
    rhoE_.correctBoundaryConditions();
}

void Foam::eulerSolver::conservativeToPrimitiveFields()
{
    U_.ref() = rhoU_/rho_;
    p_.ref() = (rhoE_-0.5*rho_*magSqr(U_))*(Gamma-1.0);
    const bool rhoBool = Foam::positiveCorrect(rho_);
    const bool pBool   = Foam::positiveCorrect(p_);
    if (rhoBool || pBool)
    {
        rhoU_.ref() = rho_*U_;
        rhoE_.ref() = p_/(Gamma-1.0)+0.5*rho_*magSqr(U_);
    }
    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    forAll(mesh_.boundary(), patchI)
    {
        const word name = mesh_.boundary()[patchI].name();
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        fvPatchScalarField& rhoBound = rho_.boundaryFieldRef()[patchI];
        fvPatchVectorField& UBound = U_.boundaryFieldRef()[patchI];
        fvPatchScalarField& pBound = p_.boundaryFieldRef()[patchI];
        if (name == "farField")
        {
            if (fluidProps_.Mach_inf < 1.0)
            {
                forAll(bfaceCells, faceI)
                {
                    const label i = bfaceCells[faceI];
                    const scalar c2 = Gamma*p_[i]/rho_[i];
                    const scalar qn = U_[i]&normal[faceI];
                    const scalar rhoc = Foam::sqrt(c2)*rho_[i]/Gamma;
                    if (qn < 0.0)
                    {
                        const scalar pB = pBound[faceI];
                        pBound[faceI] = 0.5*(pB+p_[i]-rhoc*(normal[faceI]&(UBound[faceI]-U_[i])));
                        rhoBound[faceI] += (pBound[faceI]-pB)*Gamma/c2;
                        UBound[faceI] += normal[faceI]*(pBound[faceI]-pB)/rhoc;
                    }
                    else
                    {
                        rhoBound[faceI] = rho_[i]+(pBound[faceI]-p_[i])*Gamma/c2;
                        UBound[faceI] = U_[i]+normal[faceI]*(p_[i]-pBound[faceI])/rhoc;
                    }
                }
            }
            else
            {
                forAll(bfaceCells, faceI)
                {
                    const label i = bfaceCells[faceI];
                    const scalar qn = U_[i]&normal[faceI];
                    if (qn > 0.0)
                    {
                        rhoBound[faceI] = rho_[i];
                        UBound[faceI]   = U_[i];
                        pBound[faceI]   = p_[i];
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
                    const scalar c2 = Gamma*p_[i]/rho_[i];
                    const scalar rhoc = rho_[i]*Foam::sqrt(c2)/Gamma;
                    const scalar pB = pBound[faceI];
                    pBound[faceI] = 0.5*(pB+p_[i]-rhoc*(normal[faceI]&(UBound[faceI]-U_[i])));
                    rhoBound[faceI] += (pBound[faceI]-pB)*Gamma/c2;
                    UBound[faceI] += normal[faceI]*(pBound[faceI]-pB)/rhoc;
                }
            }
        }
        if (name == "outlet")
        {
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const scalar c2 = Gamma*p_[i]/rho_[i];
                const scalar Ma = mag(U_[i])/Foam::sqrt(c2);
                const scalar rhoc = rho_[i]*Foam::sqrt(c2)/Gamma;
                if (Ma >= 1.0)
                {
                    pBound[faceI]   = p_[i];
                }
                else
                {
                    rhoBound[faceI] = rho_[i]+(pBound[faceI]-p_[i])*Gamma/c2;
                    UBound[faceI] = U_[i]+normal[faceI]*(p_[i]-pBound[faceI])/rhoc;
                    if ((UBound[faceI]&normal[faceI]) < 0.0)
                    {
                        const vector signU = vector(sign(U_[i].x()), sign(U_[i].y()), sign(U_[i].z()));
                        UBound[faceI] = cmptMultiply(signU, max(cmptMag(UBound[faceI]), cmptMag(U_[i])));
                    }
                }
            }
        }
    }
    T_ = Gamma*p_/rho_;
}