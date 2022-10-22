void Foam::euler3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    rhoGrad_ = fvc::grad(rho_);
    UGrad_   = fvc::grad(U_);
    TGrad_   = fvc::grad(T_);
    updateCoefficients();
    limitCoefficients();
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux = 0.0;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux = 0.0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        if (p0_[i] && p0_[j])
        {
            const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
            const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
            const scalar rho_i = rho_[i] + (rhoGradLimited_[i]&delta_i);
            const scalar rho_j = rho_[j] + (rhoGradLimited_[j]&delta_j);
            const vector U_i = U_[i] + (UGradLimited_[i]&delta_i);
            const vector U_j = U_[j] + (UGradLimited_[j]&delta_j);
            const scalar T_i = T_[i] + (TGradLimited_[i]&delta_i);
            const scalar T_j = T_[j] + (TGradLimited_[j]&delta_j);
            riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                   rho_i, rho_j, U_i, U_j, T_i, T_j,
                                   normal_[faceI], fluidProps_.gamma);
            rhoFlux  *= mesh_.magSf()[faceI];
            rhoUFlux *= mesh_.magSf()[faceI];
            rhoEFlux *= mesh_.magSf()[faceI];
        }
        else
        {
            for (label gaussI = 0; gaussI != 4; ++gaussI)
            {
                scalar rhoFluxTemp;
                vector rhoUFluxTemp = vector::zero;
                scalar rhoEFluxTemp;
                const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
                const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
                const scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i], basisConst_[i]);
                const scalar rho_j = rho_[j] + scalarPoly(delta_j, rhoGradLimited_[j], d2RhoLimited_[j], rLengthScale_[j], basisConst_[j]);
                const vector U_i = U_[i] + vectorPoly(delta_i, UGradLimited_[i], d2UxLimited_[i], d2UyLimited_[i], d2UzLimited_[i], rLengthScale_[i], basisConst_[i]);
                const vector U_j = U_[j] + vectorPoly(delta_j, UGradLimited_[j], d2UxLimited_[j], d2UyLimited_[j], d2UzLimited_[j], rLengthScale_[j], basisConst_[j]);
                const scalar T_i = T_[i] + scalarPoly(delta_i, TGradLimited_[i], d2TLimited_[i], rLengthScale_[i], basisConst_[i]);
                const scalar T_j = T_[j] + scalarPoly(delta_j, TGradLimited_[j], d2TLimited_[j], rLengthScale_[j], basisConst_[j]);
                evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                             rho_i, rho_j, U_i, U_j, T_i, T_j,
                             normal_[faceI]);
                rhoFlux  += rhoFluxTemp  * Quad4::w[gaussI];
                rhoUFlux += rhoUFluxTemp * Quad4::w[gaussI];
                rhoEFlux += rhoEFluxTemp * Quad4::w[gaussI];
            }
            const scalar temp = mesh_.magSf()[faceI] * 0.25;
            rhoFlux  *= temp;
            rhoUFlux *= temp;
            rhoEFlux *= temp;
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
            vectorField neighCn(bfaceCells.size());
            if (isA<processorFvPatch>(mesh_.boundary()[patchI]))
            {
                neighCn = mesh_.C().boundaryField()[patchI].patchNeighbourField();
            }
            else if (isA<cyclicFvPatch>(mesh_.boundary()[patchI]))
            {
                const cyclicFvPatch& ownerPatch = refCast<const cyclicFvPatch>(mesh_.boundary()[patchI]);
                const cyclicFvPatch& neighPatch = ownerPatch.neighbFvPatch();
                neighCn = neighPatch.Cn()-(neighPatch.Cf() - ownerPatch.Cf());
            }
            const vectorField rLengthScale_neigh = rLengthScale_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisConst_neigh = basisConst_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh   = T_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGradLimited_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGradLimited_.boundaryField()[patchI].patchNeighbourField();
            const vectorField TGrad_neigh   = TGradLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2RhoLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Ux_neigh = d2UxLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uy_neigh = d2UyLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uz_neigh = d2UzLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2T_neigh = d2TLimited_.boundaryField()[patchI].patchNeighbourField();
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
                    scalar rhoFluxTemp;
                    vector rhoUFluxTemp = vector::zero;
                    scalar rhoEFluxTemp;
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];
                    const scalar rho_i = rho_[i]      + scalarPoly(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar rho_j = rho_neigh[j] + scalarPoly(delta_j, rhoGrad_neigh[j], d2Rho_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    const vector U_i = U_[i]      + vectorPoly(delta_i, UGradLimited_[i], d2UxLimited_[i], d2UyLimited_[i], d2UzLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const vector U_j = U_neigh[j] + vectorPoly(delta_j, UGrad_neigh[j], d2Ux_neigh[j], d2Uy_neigh[j], d2Uz_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    const scalar T_i = T_[i]      + scalarPoly(delta_i, TGradLimited_[i], d2TLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar T_j = T_neigh[j] + scalarPoly(delta_j, TGrad_neigh[j], d2T_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                 rho_i, rho_j, U_i, U_j, T_i, T_j,
                                 normal[j]);
                    rhoFlux  += rhoFluxTemp  * Quad4::w[gaussI];
                    rhoUFlux += rhoUFluxTemp * Quad4::w[gaussI];
                    rhoEFlux += rhoEFluxTemp * Quad4::w[gaussI];
                }
                const scalar temp = 0.25 * magSf[j];
                resRho[i]  -= rhoFlux  * temp;
                resRhoU[i] -= rhoUFlux * temp;
                resRhoE[i] -= rhoEFlux * temp;
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
                    const scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar T_i = T_[i] + scalarPoly(delta_i, TGradLimited_[i], d2TLimited_[i], rLengthScale_[i], basisConst_[i]);
                    rhoUFlux += (T_i*rho_i/fluidProps_.gamma) * normal[faceI] * Quad4::w[gaussI];
                }
                resRhoU[i] -= 0.25 * rhoUFlux * magSf[faceI];
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
                    scalar rhoFluxTemp;
                    vector rhoUFluxTemp = vector::zero;
                    scalar rhoEFluxTemp;
                    const vector delta_i = quadPoints[gaussI] - Cn[faceI];
                    const scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const vector U_i = U_[i] + vectorPoly(delta_i, UGradLimited_[i], d2UxLimited_[i], d2UyLimited_[i], d2UzLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar T_i = T_[i] + scalarPoly(delta_i, TGradLimited_[i], d2TLimited_[i], rLengthScale_[i], basisConst_[i]);
                    evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                 rho_i, rhoBound[faceI], U_i, UBound[faceI], T_i, TBound[faceI],
                                 normal[faceI]);
                    rhoFlux  += rhoFluxTemp  * Quad4::w[gaussI];
                    rhoUFlux += rhoUFluxTemp * Quad4::w[gaussI];
                    rhoEFlux += rhoEFluxTemp * Quad4::w[gaussI];
                }
                const scalar temp = 0.25 * magSf[faceI];
                resRho[i]  -= rhoFlux  * temp;
                resRhoU[i] -= rhoUFlux * temp;
                resRhoE[i] -= rhoEFlux * temp;
            }
        }
    }
}

void Foam::euler3rdSolver::evaluateFlux
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
    const vector& normal
) const
{
    const scalar p_L = T_L*rho_L/fluidProps_.gamma;
    const scalar p_R = T_R*rho_R/fluidProps_.gamma;
    riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
}
