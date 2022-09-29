void Foam::euler2ndSolver::evaluateFlowRes
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
        const label start = mesh_.boundary()[patchI].start();
        if (isA<cyclicFvPatch>(mesh_.boundary()[patchI]))
        {
            const cyclicFvPatch& ownerPatch = refCast<const cyclicFvPatch>(mesh_.boundary()[patchI]);
            const cyclicFvPatch& neighPatch = ownerPatch.neighbFvPatch();
            const vectorField delta_i = ownerPatch.Cf() - ownerPatch.Cn();
            const vectorField delta_j = neighPatch.Cf() - neighPatch.Cn();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoLimit_neigh = rhoLimit_.boundaryField()[patchI].patchNeighbourField();
            const vectorField ULimit_neigh   = ULimit_.boundaryField()[patchI].patchNeighbourField();
            const scalarField pLimit_neigh   = pLimit_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i[faceI])*rhoLimit_[i];
                const scalar rho_j = rho_neigh[faceI] + (rhoGrad_neigh[faceI]&delta_j[faceI])*rhoLimit_neigh[faceI];
                const vector U_i = U_[i] + cmptMultiply(UGrad_[i]&delta_i[faceI], ULimit_[i]);
                const vector U_j = U_neigh[faceI] + cmptMultiply(UGrad_neigh[faceI]&delta_j[faceI], ULimit_neigh[faceI]);
                const scalar p_i = p_[i] + (pGrad_[i]&delta_i[faceI])*pLimit_[i];
                const scalar p_j = p_neigh[faceI] + (pGrad_neigh[faceI]&delta_j[faceI])*pLimit_neigh[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, p_i, p_j,
                                       normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
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
            const vectorField delta_i = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            const vectorField delta_j = mesh_.boundary()[patchI].Cf() - mesh_.C().boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i[faceI])*rhoLimit_[i];
                const scalar rho_j = rho_neigh[faceI] + (rhoGrad_neigh[faceI]&delta_j[faceI])*rhoLimit_neigh[faceI];
                const vector U_i = U_[i] + cmptMultiply((UGrad_[i]&delta_i[faceI]), ULimit_[i]);
                const vector U_j = U_neigh[faceI] + cmptMultiply(UGrad_neigh[faceI]&delta_j[faceI], ULimit_neigh[faceI]);
                const scalar p_i = p_[i] + (pGrad_[i]&delta_i[faceI])*pLimit_[i];
                const scalar p_j = p_neigh[faceI] + (pGrad_neigh[faceI]&delta_j[faceI])*pLimit_neigh[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, p_i, p_j,
                                       normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const vectorField delta = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const scalar p_i = p_[i] + (pGrad_[i]&delta[faceI])*pLimit_[i];
                vector rhoUFlux = p_i*normal[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
            }
            // std::vector<vector> quadPoints(4, vector::zero);
            // forAll(bfaceCells, faceI)
            // {
            //     const label i = bfaceCells[faceI];
            //     scalar rhoFlux  = 0.0;
            //     vector rhoUFlux = vector::zero;
            //     scalar rhoEFlux = 0.0;
            //     gaussQuad4(start+faceI, quadPoints);
            //     for (label gaussI = 0; gaussI != 4; ++gaussI)
            //     {
            //         scalar rhoFluxTemp;
            //         vector rhoUFluxTemp = vector::zero;
            //         scalar rhoEFluxTemp;
            //         const vector delta = quadPoints[gaussI] - mesh_.C()[i];
            //         const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta)*rhoLimit_[i];
            //         const scalar p_i = p_[i] + (pGrad_[i]&delta)*pLimit_[i];
            //         vector Normal = -vector(quadPoints[gaussI].x(), quadPoints[gaussI].y(), 0.0);
            //         Normal /= mag(Normal);
            //         const vector U_i = U_[i];
            //         const vector U_correct = vector(+(U_i.x()*Normal.y()-U_i.y()*Normal.x())*Normal.y(),
            //                                         -(U_i.x()*Normal.y()-U_i.y()*Normal.x())*Normal.x(), 0.0);
            //         Foam::evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
            //                            rho_i, U_correct, p_i,
            //                            normal[faceI], fluidProps_.gamma);
            //         rhoFlux  += rhoFluxTemp  * Quad4::w[gaussI];
            //         rhoUFlux += rhoUFluxTemp * Quad4::w[gaussI];
            //         rhoEFlux += rhoEFluxTemp * Quad4::w[gaussI];
            //     }
            //     const scalar temp = 0.25 * magSf[faceI];
            //     resRho[i]  -= rhoFlux  * temp;
            //     resRhoU[i] -= rhoUFlux * temp;
            //     resRhoE[i] -= rhoEFlux * temp;
            // }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const vectorField delta = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta[faceI])*rhoLimit_[i];
                const vector U_i = U_[i] + cmptMultiply((UGrad_[i]&delta[faceI]), ULimit_[i]);
                const scalar p_i = p_[i] + (pGrad_[i]&delta[faceI])*pLimit_[i];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rhoBound[faceI], U_i, UBound[faceI], p_i, pBound[faceI],
                                       normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
    }
}
