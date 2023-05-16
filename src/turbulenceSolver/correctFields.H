

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