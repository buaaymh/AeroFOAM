void Foam::turbulenceSolver::evaluateFlowResForGMRES
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    eulerSolver::evaluateFlowResForGMRES(resRho, resRhoU, resRhoE);
    UGrad_   = fvc::grad(U_);
    TGrad_   = fvc::grad(T_);
    forAll(mesh_.owner(), faceI)
    {
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux = 0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        scalar dx_L  = mag((mesh_.Cf()[faceI] - mesh_.C()[i])&normal_[faceI]) + SMALL;
        scalar dx_R  = mag((mesh_.Cf()[faceI] - mesh_.C()[j])&normal_[faceI]) + SMALL;
        scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx; 
        evaluateViscousFlux(rhoUFlux, rhoEFlux,
                            dx_L*U_[i]     +dx_R*U_[j],
                            dx_L*TGrad_[i] +dx_R*TGrad_[j],
                            dx_L*UGrad_[i] +dx_R*UGrad_[j],
                            dx_L*muLam_[i] +dx_R*muLam_[j],
                            dx_L*muTurb_[i]+dx_R*muTurb_[j], normal_[faceI]);
        resRhoU[i] += rhoUFlux*mesh_.magSf()[faceI];
        resRhoE[i] += rhoEFlux*mesh_.magSf()[faceI];
        resRhoU[j] -= rhoUFlux*mesh_.magSf()[faceI];
        resRhoE[j] -= rhoEFlux*mesh_.magSf()[faceI];
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        const word name = mesh_.boundary()[patchI].name();
        const word type = mesh_.boundary()[patchI].type();
        if (mesh_.boundary()[patchI].coupled())
        {
            const vectorField delta_i = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            const vectorField delta_j = delta_i - mesh_.boundary()[patchI].delta();
            const vectorField U_neigh = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh = T_.boundaryField()[patchI].patchNeighbourField();
            const scalarField muLam_neigh = muLam_.boundaryField()[patchI].patchNeighbourField();
            const scalarField muTurb_neigh= muTurb_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField TGrad_neigh = TGrad_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                scalar dx_L  = mag(delta_i[faceI]&normal_[faceI]) + SMALL;
                scalar dx_R  = mag(delta_j[faceI]&normal_[faceI]) + SMALL;
                scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx; 
                evaluateViscousFlux(rhoUFlux, rhoEFlux,
                                    dx_L*U_[i]     +dx_R*U_neigh[faceI],
                                    dx_L*TGrad_[i] +dx_R*TGrad_neigh[faceI],
                                    dx_L*UGrad_[i] +dx_R*UGrad_neigh[faceI],
                                    dx_L*muLam_[i] +dx_R*muLam_neigh[faceI],
                                    dx_L*muTurb_[i]+dx_R*muTurb_neigh[faceI], normal[faceI]);
                resRhoU[i] += rhoUFlux*magSf[faceI];
                resRhoE[i] += rhoEFlux*magSf[faceI];
            }
        }
        if (type == "wall")
        {
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                evaluateViscousFlux(rhoUFlux, rhoEFlux, vector::zero, vector::zero, UGrad_[i],
                                    muLam_[i], 0, normal[faceI]);
                resRhoU[i] += rhoUFlux * magSf[faceI];
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& muLamBound = muLam_.boundaryField()[patchI];
            const fvPatchScalarField& muTurbBound = muTurb_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                evaluateViscousFlux(rhoUFlux, rhoEFlux, UBound[faceI], TGrad_[i], UGrad_[i],
                                    muLamBound[faceI], muTurbBound[faceI], normal[faceI]);
                resRhoU[i] += rhoUFlux * magSf[faceI];
                resRhoE[i] += rhoEFlux * magSf[faceI];
            }
        }
    }
}
void Foam::turbulenceSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    evaluateFlowResForGMRES(resRho, resRhoU, resRhoE);
}

void Foam::turbulenceSolver::evaluateViscousFlux
(
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const vector& U,
    const vector& TGrad,
    const tensor& UGrad,
    const scalar& muLam,
    const scalar& muTurb,
    const vector& normal
) const
{
    const scalar mu = muLam + muTurb;
    const scalar lambda = -2.0/3.0*mu;
    const scalar kappa  = (muLam/Pr_Lam + muTurb/Pr_Turb)/(Gamma-1.0);
    const vector StrainRateU = (UGrad + UGrad.T())&normal;
    const scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
    const vector rhoUViscousFlux = mu*StrainRateU + lambda*divU*normal;
    rhoUFlux = Ma_Re_*rhoUViscousFlux;
    rhoEFlux = Ma_Re_*((rhoUViscousFlux&U)+kappa*(TGrad&normal));
}
