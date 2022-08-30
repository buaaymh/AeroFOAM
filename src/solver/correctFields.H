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
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const word name = mesh_.boundary()[patchI].name();
        const vectorField& normal = normal_.boundaryField()[patchI];
        fvPatchScalarField& rhoBound = rho_.boundaryFieldRef()[patchI];
        fvPatchVectorField& UBound = U_.boundaryFieldRef()[patchI];
        fvPatchScalarField& pBound = p_.boundaryFieldRef()[patchI];
        fvPatchScalarField& Ma = Ma_.boundaryFieldRef()[patchI];
        const scalarField rhoInter = rhoBound.patchInternalField();
        const vectorField UInter = UBound.patchInternalField();
        const scalarField pInter = pBound.patchInternalField();
        if (name == "inlet")
        {
            Ma = mag(UBound)/sqrt(pBound*fluidProps_.gamma/rhoBound);
            forAll(bfaceCells, faceI)
            {
                if (Ma[faceI] < 1.0)
                {   
                    const label i = bfaceCells[faceI];
                    const scalar rhoc = rhoInter[faceI]*c_[i];
                    const scalar p    = pBound[faceI];
                    pBound[faceI] = 0.5*(p+pInter[faceI]-rhoc*(normal[faceI]&(UBound[faceI]-UInter[faceI])));
                    rhoBound[faceI] += (pBound[faceI]-p)/sqr(c_[i]);
                    UBound[faceI] += normal[faceI]*(pBound[faceI]-p)/rhoc;
                }
            }
        }
        if (name == "outlet")
        {
            Ma = mag(UInter)/sqrt(pInter*fluidProps_.gamma/rhoInter);
            forAll(bfaceCells, faceI)
            {
                if (Ma[faceI] < 1.0)
                {  
                    const label i = bfaceCells[faceI];
                    const scalar rhoc = rhoInter[faceI]*c_[i];
                    rhoBound[faceI] += (pBound[faceI]-pInter[faceI])/sqr(c_[i]);
                    UBound[faceI] += normal[faceI]*(pInter[faceI]-pBound[faceI])/rhoc;
                }
            }
        }
        if (name == "farField")
        {
            const scalarField Vn = normal&UBound;
            Ma = mag(UBound)/sqrt(pBound*fluidProps_.gamma/rhoBound);
            forAll(bfaceCells, faceI)
            {
                if (Ma[faceI] < 1.0)
                {
                    const label i = bfaceCells[faceI];
                    const scalar c2 = T_[i];
                    const scalar rhoc = rhoInter[faceI]*c_[i];
                    if (Vn[faceI] < 0)
                    {
                        const scalar p = pBound[faceI];
                        pBound[faceI] = 0.5*(p+pInter[faceI]-rhoc*(normal[faceI]&(UBound[faceI]-UInter[faceI])));
                        rhoBound[faceI] += (pBound[faceI]-p)/c2;
                        UBound[faceI] += normal[faceI]*(pBound[faceI]-p)/rhoc;
                    }
                    else
                    {
                        rhoBound[faceI] = rhoInter[faceI]+(pBound[faceI]-pInter[faceI])/c2;
                        UBound[faceI] = UInter[faceI]+normal[faceI]*(pInter[faceI]-pBound[faceI])/rhoc;
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