void Foam::solver::correctFields()
{
    // Internal fields
    U_.ref() = rhoU_ / rho_;
    p_.ref() = (rhoE_ - 0.5*rho_*magSqr(U_)) * (fluidProps_.gamma-1.0);
    const bool rhoBool = Foam::positiveCorrect(rho_);
    const bool pBool   = Foam::positiveCorrect(p_);
    if (rhoBool || pBool)
    {
        rhoU_ = rho_ * U_;
        rhoE_ = p_/(fluidProps_.gamma-1.0) + 0.5*rho_*magSqr(U_);
    }
    T_.ref() = p_*fluidProps_.gamma/rho_;
    c_ = sqrt(T_.primitiveField());
    Ma_.primitiveFieldRef() = mag(U_.primitiveFieldRef())/c_;
    // Boundary conditions
    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    T_.correctBoundaryConditions();
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
                    const scalar qn = U_[i]&normal[faceI];
                    const scalar rhoc = c_[i]*rho_[i]/fluidProps_.gamma;
                    if (qn < 0.0)
                    {
                        const scalar p = pBound[faceI];
                        pBound[faceI] = 0.5*(p+p_[i]-rhoc*(normal[faceI]&(UBound[faceI]-U_[i])));
                        rhoBound[faceI] += (pBound[faceI]-p)*fluidProps_.gamma/sqr(c_[i]);
                        UBound[faceI] += normal[faceI]*(pBound[faceI]-p)/rhoc;
                    }
                    else
                    {
                        rhoBound[faceI] = rho_[i]+(pBound[faceI]-p_[i])*fluidProps_.gamma/sqr(c_[i]);
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
            const scalarField Ma = mag(UBound)/sqrt(pBound*fluidProps_.gamma/rhoBound);
            forAll(bfaceCells, faceI)
            {
                if (Ma[faceI] < 1.0)
                {
                    const label i = bfaceCells[faceI];
                    const scalar rhoc = rho_[i]*c_[i]/fluidProps_.gamma;
                    const scalar p = pBound[faceI];
                    pBound[faceI] = 0.5*(p+p_[i]-rhoc*(normal[faceI]&(UBound[faceI]-U_[i])));
                    rhoBound[faceI] += (pBound[faceI]-p)*fluidProps_.gamma/sqr(c_[i]);
                    UBound[faceI] += normal[faceI]*(pBound[faceI]-p)/rhoc;
                }
            }
        }
        if (name == "outlet")
        {
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const scalar Ma = mag(U_[i])/c_[i];
                const scalar rhoc = rho_[i]*c_[i]/fluidProps_.gamma;
                if (Ma >= 1.0)
                {
                    pBound[faceI]   = p_[i];
                }
                else
                {
                    rhoBound[faceI] = rho_[i]+(pBound[faceI]-p_[i])*fluidProps_.gamma/sqr(c_[i]);
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
    // Update other BCs
    const volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U_.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p_.boundaryFieldRef();
    volVectorField::Boundary& rhoUBf = rhoU_.boundaryFieldRef();
    volScalarField::Boundary& rhoEBf = rhoE_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T_.boundaryFieldRef();
    volScalarField::Boundary& MaBf = Ma_.boundaryFieldRef();
    rhoUBf = rhoBf * UBf;
    rhoEBf = pBf/(fluidProps_.gamma-1.0) + 0.5*rhoBf*magSqr(UBf);
    TBf = pBf*fluidProps_.gamma/rhoBf;
    MaBf = mag(UBf)/sqrt(TBf);
}