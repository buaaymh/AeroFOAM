void Foam::eulerSolver::correctFields()
{
    const dimensionedScalar boundMin("boundMin", dimless, 0.001);
    const dimensionedScalar boundMax("boundMax", dimless, 100.0);
    U_.ref() = rhoU_ / rho_;
    p_.ref() = (rhoE_ - 0.5*rho_*magSqr(U_))*(fluidProps_.gamma-1.0);
    rho_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    p_.correctBoundaryConditions();
    const bool rhoBool = Foam::boundMinMax(rho_, boundMin, boundMax);
    const bool pBool   = Foam::boundMinMax(p_, boundMin, boundMax);
    if (rhoBool || pBool)
    {
        rhoU_ = rho_ * U_;
        rhoE_ = p_/(fluidProps_.gamma-1.0) + 0.5*rho_*magSqr(U_);
    }
    T_.ref() = p_ * fluidProps_.gamma / rho_;
    c_.ref() = sqrt(T_);
    Ma_.ref() = mag(U_)/c_;
    const volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
    const volVectorField::Boundary& UBf = U_.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p_.boundaryFieldRef();
    volVectorField::Boundary& rhoUBf = rhoU_.boundaryFieldRef();
    volScalarField::Boundary& rhoEBf = rhoE_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T_.boundaryFieldRef();
    volScalarField::Boundary& cBf = c_.boundaryFieldRef();
    volScalarField::Boundary& MaBf = Ma_.boundaryFieldRef();
    rhoUBf = rhoBf * UBf;
    rhoEBf = pBf/(fluidProps_.gamma-1.0) + 0.5*rhoBf*magSqr(UBf);
    TBf = pBf * fluidProps_.gamma / rhoBf;
    cBf = Foam::sqrt(TBf);
    MaBf = mag(UBf)/cBf;
}

