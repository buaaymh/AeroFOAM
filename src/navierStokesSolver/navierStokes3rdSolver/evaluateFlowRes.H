void Foam::navierStokes3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    scalar rho_i, rho_j, rhoE_i, rhoE_j;
    vector rhoU_i, rhoU_j, rhoGrad_i, rhoGrad_j, rhoEGrad_i, rhoEGrad_j;
    tensor rhoUGrad_i, rhoUGrad_j;
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    iterationStep();
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux = 0.0;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux = 0.0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const scalar deltaTilde = min(mesh_.V()[i], mesh_.V()[j])/mesh_.magSf()[faceI];
        for (label gaussI = 0; gaussI != faceQuad_[faceI]->size(); ++gaussI)
        {
            const vector delta_i = faceQuad_[faceI]->at(gaussI) - mesh_.C()[i];
            const vector delta_j = faceQuad_[faceI]->at(gaussI) - mesh_.C()[j];
            evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i,
                                 rho_[i], rhoU_[i], rhoE_[i], coefs_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
            evaluateVarsAndGrads(rho_j, rhoU_j, rhoE_j, rhoGrad_j, rhoUGrad_j, rhoEGrad_j,
                                 rho_[j], rhoU_[j], rhoE_[j], coefs_[j], basisMean_[j], rDeltaXYZ_[j], delta_j);
            evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j,
                         rhoGrad_i, rhoGrad_j, rhoUGrad_i, rhoUGrad_j, rhoEGrad_i, rhoEGrad_j,
                         deltaTilde, faceQuad_[faceI]->normal(gaussI), faceQuad_[faceI]->weight(gaussI));
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
            const scalarField rho_neigh  = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoU_neigh = rhoU_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoE_neigh = rhoE_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1Rho_neigh  = d1Rho_.boundaryField()[patchI].patchNeighbourField();
            const tensorField d1RhoU_neigh = d1RhoU_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1RhoE_neigh = d1RhoE_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUx_neigh = d2RhoUx_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUy_neigh = d2RhoUy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUz_neigh = d2RhoUz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoE_neigh = d2RhoE_.boundaryField()[patchI].patchNeighbourField();
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch2ndQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[faceI];
                const scalar deltaTilde = mesh_.V()[i]/magSf[faceI];
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - ownerCn[faceI];
                    const vector delta_j = boundaryQuad[faceI]->at(gaussI) - neighCn[faceI];
                    evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i,
                                        rho_[i], rhoU_[i], rhoE_[i], d1Rho_[i], d1RhoU_[i], d1RhoE_[i],
                                        d2Rho_[i], d2RhoUx_[i], d2RhoUy_[i], d2RhoUz_[i], d2RhoE_[i],
                                        basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateVarsAndGrads(rho_j, rhoU_j, rhoE_j, rhoGrad_j, rhoUGrad_j, rhoEGrad_j,
                                        rho_neigh[faceI], rhoU_neigh[faceI], rhoE_neigh[faceI], d1Rho_neigh[faceI], d1RhoU_neigh[faceI], d1RhoE_neigh[faceI],
                                        d2Rho_neigh[faceI], d2RhoUx_neigh[faceI], d2RhoUy_neigh[faceI], d2RhoUz_neigh[faceI], d2RhoE_neigh[faceI],
                                        basisMean_neigh[faceI], rDeltaXYZ_neigh[faceI], delta_j);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j,
                                 rhoGrad_i, rhoGrad_j, rhoUGrad_i, rhoUGrad_j, rhoEGrad_i, rhoEGrad_j,
                                 deltaTilde, boundaryQuad[faceI]->normal(gaussI), boundaryQuad[faceI]->weight(gaussI));
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
        if (type == "symmetryPlane" || type == "symmetry")
        {
            const vectorField  Cn = mesh_.boundary()[patchI].Cn();
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch2ndQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux = vector::zero;
                const label i = bfaceCells[faceI];
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - Cn[faceI];
                    evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], coefs_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    const scalar p = (rhoE_i-0.5*magSqr(rhoU_i)/rho_i)*(Gamma-1.0);
                    rhoUFlux += p * boundaryQuad[faceI]->normal(gaussI) * boundaryQuad[faceI]->weight(gaussI);
                }
                resRhoU[i] -= rhoUFlux;
            }
        }
        if (type == "wall")
        {
            const vectorField  Cn = mesh_.boundary()[patchI].Cn();
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch2ndQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux = vector::zero;
                const label i = bfaceCells[faceI];
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - Cn[faceI];
                    evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i,
                                         rho_[i], rhoU_[i], rhoE_[i], coefs_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateWallFlux(rhoUFlux, rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i,
                                     boundaryQuad[faceI]->normal(gaussI), boundaryQuad[faceI]->weight(gaussI));
                }
                resRhoU[i] -= rhoUFlux;
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const vectorField  Cn = mesh_.boundary()[patchI].Cn();
            std::unique_ptr<Face> faceQuad;
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[faceI];
                build2ndFace(mesh_, start+faceI, faceQuad);
                for (label gaussI = 0; gaussI != faceQuad->size(); ++gaussI)
                {
                    const vector delta_i = faceQuad->at(gaussI) - Cn[faceI];
                    evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i,
                                         rho_[i], rhoU_[i], rhoE_[i], coefs_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateInletOutletFlux(rhoFlux, rhoUFlux, rhoEFlux, rhoBound[faceI], UBound[faceI], pBound[faceI],
                                            rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i,
                                            faceQuad->normal(gaussI), faceQuad->weight(gaussI));
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
    }
}