void Foam::euler3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    scalar rho_i, rho_j, rhoE_i, rhoE_j;
    vector rhoU_i, rhoU_j;
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
        for (label gaussI = 0; gaussI != faceQuad_[faceI]->size(); ++gaussI)
        {
            const vector delta_i = faceQuad_[faceI]->at(gaussI) - mesh_.C()[i];
            const vector delta_j = faceQuad_[faceI]->at(gaussI) - mesh_.C()[j];
            evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], coefs_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
            evaluateVars(rho_j, rhoU_j, rhoE_j, rho_[j], rhoU_[j], rhoE_[j], coefs_[j], basisMean_[j], rDeltaXYZ_[j], delta_j);
            evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j,
                            faceQuad_[faceI]->normal(gaussI), faceQuad_[faceI]->weight(gaussI));
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
            forAll(bfaceCells, j)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[j];
                for (label gaussI = 0; gaussI != boundaryQuad[j]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[j]->at(gaussI) - ownerCn[j];
                    const vector delta_j = boundaryQuad[j]->at(gaussI) - neighCn[j];
                    evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], d1Rho_[i], d1RhoU_[i], d1RhoE_[i],
                                 d2Rho_[i], d2RhoUx_[i], d2RhoUy_[i], d2RhoUz_[i], d2RhoE_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateVars(rho_j, rhoU_j, rhoE_j, rho_neigh[j], rhoU_neigh[j], rhoE_neigh[j], d1Rho_neigh[j], d1RhoU_neigh[j], d1RhoE_neigh[j],
                                 d2Rho_neigh[j], d2RhoUx_neigh[j], d2RhoUy_neigh[j], d2RhoUz_neigh[j], d2RhoE_neigh[j], basisMean_neigh[j], rDeltaXYZ_neigh[j], delta_j);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j,
                                 boundaryQuad[j]->normal(gaussI), boundaryQuad[j]->weight(gaussI));
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const vectorField  Cn = mesh_.boundary()[patchI].Cn();
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch2ndQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[faceI];
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - Cn[faceI];
                    evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], coefs_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    rhoU_j = rhoU_i - 2*(rhoU_i&boundaryQuad[faceI]->normal(gaussI))*boundaryQuad[faceI]->normal(gaussI);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, rho_i, rho_i, rhoU_i, rhoU_j, rhoE_i, rhoE_i,
                                 boundaryQuad[faceI]->normal(gaussI), boundaryQuad[faceI]->weight(gaussI));
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
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
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                Foam::evaluateConvectiveFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                             rhoBound[faceI], UBound[faceI], pBound[faceI],
                                             normal[faceI]*magSf[faceI]);
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
    }
}
