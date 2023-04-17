void Foam::solver::correctFields()
{
    // Internal fields
    U_.ref() = rhoU_ / rho_;
    p_.ref() = (rhoE_ - 0.5*rho_*magSqr(U_)) * (Gamma-1.0);
    const bool rhoBool = Foam::positiveCorrect(rho_);
    const bool pBool   = Foam::positiveCorrect(p_);
    if (rhoBool || pBool)
    {
        rhoU_ = rho_ * U_;
        rhoE_ = p_/(Gamma-1.0) + 0.5*rho_*magSqr(U_);
    }
    T_.ref() = p_*Gamma/rho_;
    c_ = sqrt(T_.primitiveField());
    Ma_.primitiveFieldRef() = mag(U_.primitiveFieldRef())/c_;
    // Boundary conditions
    correctPrimitiveBoundary(rho_, U_, p_, T_);
    // Update other BCs
    const volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U_.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p_.boundaryFieldRef();
    const volScalarField::Boundary& TBf = T_.boundaryFieldRef();
    volVectorField::Boundary& rhoUBf = rhoU_.boundaryFieldRef();
    volScalarField::Boundary& rhoEBf = rhoE_.boundaryFieldRef();
    volScalarField::Boundary& MaBf = Ma_.boundaryFieldRef();
    rhoUBf = rhoBf * UBf;
    rhoEBf = pBf/(Gamma-1.0) + 0.5*rhoBf*magSqr(UBf);
    MaBf = mag(UBf)/sqrt(TBf);
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