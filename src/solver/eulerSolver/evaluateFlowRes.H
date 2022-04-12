void Foam::eulerSolver::evaluateFlowRes
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
    pGrad_   = fvc::grad(p_);
    limitGrad();
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const scalar magSf = mesh_.magSf()[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i)*rhoLimit_[i];
        const scalar rho_j = rho_[j] + (rhoGrad_[j]&delta_j)*rhoLimit_[j];
        const vector U_i = U_[i] + cmptMultiply(UGrad_[i]&delta_i, ULimit_[i]);
        const vector U_j = U_[j] + cmptMultiply(UGrad_[j]&delta_j, ULimit_[j]);
        const scalar p_i = p_[i] + (pGrad_[i]&delta_i)*pLimit_[i];
        const scalar p_j = p_[j] + (pGrad_[j]&delta_j)*pLimit_[j];
        riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                               rho_i, rho_j, U_i, U_j, p_i, p_j,
                               normal_[faceI], fluidProps_.gamma);
        rhoFlux  *= magSf;
        rhoUFlux *= magSf;
        rhoEFlux *= magSf;
        resRho[i]  -= rhoFlux;
        resRhoU[i] -= rhoUFlux;
        resRhoE[i] -= rhoEFlux;
        resRho[j]  += rhoFlux;
        resRhoU[j] += rhoUFlux;
        resRhoE[j] += rhoEFlux;
    }
    rhoLimit_.correctBoundaryConditions();
    ULimit_.correctBoundaryConditions();
    pLimit_.correctBoundaryConditions();
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        const word name = mesh_.boundary()[patchI].name();
        const word type = mesh_.boundary()[patchI].type();
        if (isA<cyclicFvPatch>(mesh_.boundary()[patchI]))
        {
            const cyclicFvPatch& ownerPatch = refCast<const cyclicFvPatch>(mesh_.boundary()[patchI]);
            const cyclicFvPatch& neighPatch = ownerPatch.neighbFvPatch();
            const vectorField ownerCn = ownerPatch.Cn();
            const vectorField neighCn = neighPatch.Cn();
            const vectorField& ownerCf = ownerPatch.Cf();
            const vectorField& neighCf = neighPatch.Cf();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoLimit_neigh = rhoLimit_.boundaryField()[patchI].patchNeighbourField();
            const vectorField ULimit_neigh   = ULimit_.boundaryField()[patchI].patchNeighbourField();
            const scalarField pLimit_neigh   = pLimit_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, j)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[j];
                const vector delta_i = ownerCf[j] - ownerCn[j];
                const vector delta_j = neighCf[j] - neighCn[j];
                const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i)*rhoLimit_[i];
                const scalar rho_j = rho_neigh[j] + (rhoGrad_neigh[j]&delta_j)*rhoLimit_neigh[j];
                const vector U_i = U_[i] + cmptMultiply(UGrad_[i]&delta_i, ULimit_[i]);
                const vector U_j = U_neigh[j] + cmptMultiply(UGrad_neigh[j]&delta_j, ULimit_neigh[j]);
                const scalar p_i = p_[i] + (pGrad_[i]&delta_i)*pLimit_[i];
                const scalar p_j = p_neigh[j] + (pGrad_neigh[j]&delta_j)*pLimit_neigh[j];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, p_i, p_j,
                                       normal[j], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[j];
                resRhoU[i] -= rhoUFlux * magSf[j];
                resRhoE[i] -= rhoEFlux * magSf[j];
            }
        }
        if (isA<processorFvPatch>(mesh_.boundary()[patchI]))
        {
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoLimit_neigh = rhoLimit_.boundaryField()[patchI].patchNeighbourField();
            const vectorField ULimit_neigh   = ULimit_.boundaryField()[patchI].patchNeighbourField();
            const scalarField pLimit_neigh   = pLimit_.boundaryField()[patchI].patchNeighbourField();
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const vectorField neighCn = mesh_.C().boundaryField()[patchI].patchNeighbourField();
            const vectorField& ownerCf = mesh_.boundary()[patchI].Cf();
            forAll(bfaceCells, j)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[j];
                const vector delta_i = ownerCf[j] - ownerCn[j];
                const vector delta_j = ownerCf[j] - neighCn[j];
                const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i)*rhoLimit_[i];
                const scalar rho_j = rho_neigh[j] + (rhoGrad_neigh[j]&delta_j)*rhoLimit_neigh[j];
                const vector U_i = U_[i] + cmptMultiply((UGrad_[i]&delta_i), ULimit_[i]);
                const vector U_j = U_neigh[j] + cmptMultiply(UGrad_neigh[j]&delta_j, ULimit_neigh[j]);
                const scalar p_i = p_[i] + (pGrad_[i]&delta_i)*pLimit_[i];
                const scalar p_j = p_neigh[j] + (pGrad_neigh[j]&delta_j)*pLimit_neigh[j];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, p_i, p_j,
                                       normal[j], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[j];
                resRhoU[i] -= rhoUFlux * magSf[j];
                resRhoE[i] -= rhoEFlux * magSf[j];
            }
        }
        if (type == "symmetryPlane" || type == "wall")
        {
            const scalarField& p_patch = p_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux = p_patch[faceI]*normal[faceI];
                const label i = bfaceCells[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
            }
        }
        if (name == "subInlet")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const scalarField rhoInter = rhoBound.patchInternalField();
            const vectorField UInter = UBound.patchInternalField();
            const scalarField pInter = pBound.patchInternalField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rhoInter[faceI], rhoBound[faceI],
                                       UInter[faceI], UBound[faceI],
                                       pInter[faceI], pBound[faceI],
                                       normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
        if (name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const scalarField rhoInter = rhoBound.patchInternalField();
            const vectorField UInter = UBound.patchInternalField();
            const scalarField pInter = pBound.patchInternalField();
            const scalarField& Vn = normal&UBound;
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const scalar machInf2 = magSqr(UBound[faceI]);
                const label i = bfaceCells[faceI];
                if (machInf2 < 1.0)
                {
                    if (Vn[faceI] < 0)
                        riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                               rhoInter[faceI], rhoBound[faceI],
                                               UInter[faceI], UBound[faceI],
                                               pInter[faceI], pBound[faceI],
                                               normal[faceI], fluidProps_.gamma);
                    else
                        Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                           rhoInter[faceI], UInter[faceI], pBound[faceI],
                                           normal[faceI], fluidProps_.gamma);
                }
                else
                {
                    if (Vn[faceI] < 0)
                        Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                           rhoBound[faceI], UBound[faceI], pBound[faceI],
                                           normal[faceI], fluidProps_.gamma);
                    else
                        Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                           rhoInter[faceI], UInter[faceI], pInter[faceI],
                                           normal[faceI], fluidProps_.gamma);
                }
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
        if (name == "supInlet" || name == "supOutlet" || name == "subOutlet")
        {
            const scalarField& rho_patch = rho_.boundaryField()[patchI];
            const vectorField& U_patch = U_.boundaryField()[patchI];
            const scalarField& p_patch = p_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                   rho_patch[faceI], U_patch[faceI], p_patch[faceI],
                                   normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
    }
}
