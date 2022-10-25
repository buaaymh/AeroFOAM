void Foam::euler2ndSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    volVectorField rhoGrad_ = fvc::grad(rho_);
    volTensorField UGrad_   = fvc::grad(U_);
    volVectorField TGrad_   = fvc::grad(T_);
    limitGrad(rho_, U_, T_, rhoGrad_, UGrad_, TGrad_, rhoLimit_, ULimit_, TLimit_);
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux;
        vector rhoUFlux;
        scalar rhoEFlux;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i)*rhoLimit_[i];
        const scalar rho_j = rho_[j] + (rhoGrad_[j]&delta_j)*rhoLimit_[j];
        const vector U_i = U_[i] + cmptMultiply(UGrad_[i]&delta_i, ULimit_[i]);
        const vector U_j = U_[j] + cmptMultiply(UGrad_[j]&delta_j, ULimit_[j]);
        const scalar T_i = T_[i] + (TGrad_[i]&delta_i)*TLimit_[i];
        const scalar T_j = T_[j] + (TGrad_[j]&delta_j)*TLimit_[j];
        evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                     rho_i, rho_j, U_i, U_j, T_i, T_j,
                     normal_[faceI], mesh_.magSf()[faceI]);
        resRho[i]  -= rhoFlux;
        resRhoU[i] -= rhoUFlux;
        resRhoE[i] -= rhoEFlux;
        resRho[j]  += rhoFlux;
        resRhoU[j] += rhoUFlux;
        resRhoE[j] += rhoEFlux;
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
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh   = T_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField TGrad_neigh   = TGrad_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoLimit_neigh = rhoLimit_.boundaryField()[patchI].patchNeighbourField();
            const vectorField ULimit_neigh   = ULimit_.boundaryField()[patchI].patchNeighbourField();
            const scalarField TLimit_neigh   = TLimit_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i[faceI])*rhoLimit_[i];
                const scalar rho_j = rho_neigh[faceI] + (rhoGrad_neigh[faceI]&delta_j[faceI])*rhoLimit_neigh[faceI];
                const vector U_i = U_[i] + cmptMultiply(UGrad_[i]&delta_i[faceI], ULimit_[i]);
                const vector U_j = U_neigh[faceI] + cmptMultiply(UGrad_neigh[faceI]&delta_j[faceI], ULimit_neigh[faceI]);
                const scalar T_i = T_[i] + (TGrad_[i]&delta_i[faceI])*TLimit_[i];
                const scalar T_j = T_neigh[faceI] + (TGrad_neigh[faceI]&delta_j[faceI])*TLimit_neigh[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, T_i, T_j,
                                       normal[faceI], fluidProps_.gamma);
                evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_i, rho_j, U_i, U_j, T_i, T_j,
                             normal[faceI], magSf[faceI]);
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                vector rhoUFlux = pBound[faceI]*normal[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                   rhoBound[faceI], UBound[faceI], pBound[faceI],
                                   normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
    }
}

void Foam::euler2ndSolver::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& U_L,
    const vector& U_R,
    const scalar& T_L,
    const scalar& T_R,
    const vector& normal,
    const scalar& magSf
) const
{
    const scalar p_L = T_L*rho_L/fluidProps_.gamma;
    const scalar p_R = T_R*rho_R/fluidProps_.gamma;
    riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
    rhoFlux *= magSf;
    rhoUFlux *= magSf;
    rhoEFlux *= magSf;
}