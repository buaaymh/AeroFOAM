void Foam::eulerPrimVar3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    scalar rho_i, rho_j, T_i, T_j;
    vector U_i, U_j;
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    rhoGrad_ = fvc::grad(rho_);
    UGrad_   = fvc::grad(U_);
    TGrad_   = fvc::grad(T_);
    reconstructionIter();
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux = 0.0;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux = 0.0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        if (isP0Cell_[i] && isP0Cell_[j])
        {
            const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
            const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
            const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i);
            const scalar rho_j = rho_[j] + (rhoGrad_[j]&delta_j);
            const vector U_i = U_[i] + (UGrad_[i]&delta_i);
            const vector U_j = U_[j] + (UGrad_[j]&delta_j);
            const scalar T_i = T_[i] + (TGrad_[i]&delta_i);
            const scalar T_j = T_[j] + (TGrad_[j]&delta_j);
            evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                         rho_i, rho_j, U_i, U_j, T_i, T_j,
                         normal_[faceI], mesh_.magSf()[faceI]);
        }
        else
        {
            for (label gaussI = 0; gaussI != 4; ++gaussI)
            {
                const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
                const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
                evaluateVars(rho_i, U_i, T_i, rho_[i], U_[i], T_[i], rhoGrad_[i], UGrad_[i], TGrad_[i],
                             d2Rho_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], d2T_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                evaluateVars(rho_j, U_j, T_j, rho_[j], U_[j], T_[j], rhoGrad_[j], UGrad_[j], TGrad_[j],
                             d2Rho_[j], d2Ux_[j], d2Uy_[j], d2Uz_[j], d2T_[j], basisMean_[j], rDeltaXYZ_[j], delta_j);
                evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_i, rho_j, U_i, U_j, T_i, T_j,
                             normal_[faceI], mesh_.magSf()[faceI]*Quad4::w[gaussI]);
            }
        }
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
        const label start = mesh_.boundary()[patchI].start();
        if (mesh_.boundary()[patchI].coupled())
        {
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const vectorField neighCn = ownerCn + mesh_.boundary()[patchI].delta();
            const vectorField rDeltaXYZ_neigh = rDeltaXYZ_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisMean_neigh = basisMean_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh   = T_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField TGrad_neigh   = TGrad_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Ux_neigh = d2Ux_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uy_neigh = d2Uy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uz_neigh = d2Uz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2T_neigh = d2T_.boundaryField()[patchI].patchNeighbourField();
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, j)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[j];
                gaussQuad4(start+j, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];
                    evaluateVars(rho_i, U_i, T_i, rho_[i], U_[i], T_[i], rhoGrad_[i], UGrad_[i], TGrad_[i],
                                 d2Rho_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], d2T_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateVars(rho_j, U_j, T_j, rho_neigh[j], U_neigh[j], T_neigh[j], rhoGrad_neigh[j], UGrad_neigh[j], TGrad_neigh[j],
                                 d2Rho_neigh[j], d2Ux_neigh[j], d2Uy_neigh[j], d2Uz_neigh[j], d2T_neigh[j], basisMean_neigh[j], rDeltaXYZ_neigh[j], delta_j);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                 rho_i, rho_j, U_i, U_j, T_i, T_j,
                                 normal[j], magSf[j]*Quad4::w[gaussI]);
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                vector rhoUFlux = vector::zero;
                gaussQuad4(start+faceI, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - Cn[faceI];
                    rho_i = rho_[i] + scalarPoly(delta_i, rhoGrad_[i], d2Rho_[i], rDeltaXYZ_[i], basisMean_[i]);
                    T_i   = T_[i]   + scalarPoly(delta_i, TGrad_[i], d2T_[i], rDeltaXYZ_[i], basisMean_[i]);
                    rhoUFlux += (T_i*rho_i/fluidProps_.gamma)*normal[faceI]*Quad4::w[gaussI];
                }
                resRhoU[i] -= rhoUFlux*magSf[faceI];
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& TBound = T_.boundaryField()[patchI];
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[faceI]; 
                gaussQuad4(start+faceI, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - Cn[faceI];
                    evaluateVars(rho_i, U_i, T_i, rho_[i], U_[i], T_[i], rhoGrad_[i], UGrad_[i], TGrad_[i],
                                 d2Rho_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], d2T_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                 rho_i, rhoBound[faceI], U_i, UBound[faceI], T_i, TBound[faceI],
                                 normal[faceI], magSf[faceI]*Quad4::w[gaussI]);
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
    }
}

void Foam::eulerPrimVar3rdSolver::evaluateFlux
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
    scalar rhoFluxTemp;
    vector rhoUFluxTemp;
    scalar rhoEFluxTemp;
    const scalar p_L = T_L*rho_L/fluidProps_.gamma;
    const scalar p_R = T_R*rho_R/fluidProps_.gamma;
    riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
    rhoFlux  += magSf * rhoFluxTemp;
    rhoUFlux += magSf * rhoUFluxTemp;
    rhoEFlux += magSf * rhoEFluxTemp;
}
