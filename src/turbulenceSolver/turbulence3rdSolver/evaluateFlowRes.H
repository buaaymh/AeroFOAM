void Foam::turbulence3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE,
    scalarField& resNuTilda
)
{
    scalar rho_i, rho_j, rhoE_i, rhoE_j, nuTilda_i, nuTilda_j;
    vector rhoU_i, rhoU_j, rhoGrad_i, rhoGrad_j, rhoEGrad_i, rhoEGrad_j, nuTildaGrad_i, nuTildaGrad_j;
    tensor rhoUGrad_i, rhoUGrad_j;
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    resNuTilda = scalarField(mesh_.nCells(), 0.0);
    iterationStep();
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux = 0.0;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux = 0.0;
        scalar nuTildaFlux = 0.0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const scalar deltaTilde = min(mesh_.V()[i], mesh_.V()[j])/mesh_.magSf()[faceI];
        for (label gaussI = 0; gaussI != faceQuad_[faceI]->size(); ++gaussI)
        {
            const vector delta_i = faceQuad_[faceI]->at(gaussI) - mesh_.C()[i];
            const vector delta_j = faceQuad_[faceI]->at(gaussI) - mesh_.C()[j];
            evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, nuTilda_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i, nuTildaGrad_i,
                                 rho_[i], rhoU_[i], rhoE_[i], nuTilda_[i], coefs_[i], coefsNuTilda_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
            evaluateVarsAndGrads(rho_j, rhoU_j, rhoE_j, nuTilda_j, rhoGrad_j, rhoUGrad_j, rhoEGrad_j, nuTildaGrad_j,
                                 rho_[j], rhoU_[j], rhoE_[j], nuTilda_[j], coefs_[j], coefsNuTilda_[j], basisMean_[j], rDeltaXYZ_[j], delta_j);
            evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, nuTildaFlux, 
                         rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j, nuTilda_i, nuTilda_j,
                         rhoGrad_i, rhoGrad_j, rhoUGrad_i, rhoUGrad_j, rhoEGrad_i, rhoEGrad_j, nuTildaGrad_i, nuTildaGrad_j,
                         deltaTilde, faceQuad_[faceI]->normal(gaussI), faceQuad_[faceI]->weight(gaussI));
        }
        resRho[i]  -= rhoFlux;
        resRhoU[i] -= rhoUFlux;
        resRhoE[i] -= rhoEFlux;
        resNuTilda[i] -= nuTildaFlux;
        resRho[j]  += rhoFlux;
        resRhoU[j] += rhoUFlux;
        resRhoE[j] += rhoEFlux;
        resNuTilda[j] += nuTildaFlux;
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
            const scalarField nuTilda_neigh= nuTilda_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1Rho_neigh  = d1Rho_.boundaryField()[patchI].patchNeighbourField();
            const tensorField d1RhoU_neigh = d1RhoU_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1RhoE_neigh = d1RhoE_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1NuTilda_neigh = d1NuTilda_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUx_neigh = d2RhoUx_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUy_neigh = d2RhoUy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUz_neigh = d2RhoUz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoE_neigh = d2RhoE_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2NuTilda_neigh = d2NuTilda_.boundaryField()[patchI].patchNeighbourField();
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch2ndQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                scalar nuTildaFlux = 0.0;
                const label i = bfaceCells[faceI];
                const scalar deltaTilde = mesh_.V()[i]/magSf[faceI];
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - ownerCn[faceI];
                    const vector delta_j = boundaryQuad[faceI]->at(gaussI) - neighCn[faceI];
                    evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, nuTilda_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i, nuTildaGrad_i,
                                         rho_[i], rhoU_[i], rhoE_[i], nuTilda_[i], d1Rho_[i], d1RhoU_[i], d1RhoE_[i], d1NuTilda_[i],
                                         d2Rho_[i], d2RhoUx_[i], d2RhoUy_[i], d2RhoUz_[i], d2RhoE_[i], d2NuTilda_[i],
                                         basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateVarsAndGrads(rho_j, rhoU_j, rhoE_j, nuTilda_j, rhoGrad_j, rhoUGrad_j, rhoEGrad_j, nuTildaGrad_j,
                                        rho_neigh[faceI], rhoU_neigh[faceI], rhoE_neigh[faceI],  nuTilda_neigh[faceI],
                                        d1Rho_neigh[faceI], d1RhoU_neigh[faceI], d1RhoE_neigh[faceI], d1NuTilda_neigh[faceI],
                                        d2Rho_neigh[faceI], d2RhoUx_neigh[faceI], d2RhoUy_neigh[faceI], d2RhoUz_neigh[faceI], d2RhoE_neigh[faceI], d2NuTilda_neigh[faceI],
                                        basisMean_neigh[faceI], rDeltaXYZ_neigh[faceI], delta_j);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux, nuTildaFlux,
                                 rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j, nuTilda_i, nuTilda_j,
                                 rhoGrad_i, rhoGrad_j, rhoUGrad_i, rhoUGrad_j, rhoEGrad_i, rhoEGrad_j, nuTildaGrad_i, nuTildaGrad_j,
                                 deltaTilde, boundaryQuad[faceI]->normal(gaussI), boundaryQuad[faceI]->weight(gaussI));
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
                resNuTilda[i] -= nuTildaFlux;
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
                scalar nuTildaFlux = 0.0;
                const label i = bfaceCells[faceI];
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - Cn[faceI];
                    evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, nuTilda_i,
                                         rhoGrad_i, rhoUGrad_i, rhoEGrad_i, nuTildaGrad_i,
                                         rho_[i], rhoU_[i], rhoE_[i], nuTilda_[i], coefs_[i], coefsNuTilda_[i],
                                         basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateWallFlux(rhoUFlux, nuTildaFlux,
                                     rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i, nuTildaGrad_i,
                                     boundaryQuad[faceI]->normal(gaussI), boundaryQuad[faceI]->weight(gaussI));
                }
                resRhoU[i] -= rhoUFlux;
                resNuTilda[i] -= nuTildaFlux;
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const fvPatchScalarField& nuTildaBound = nuTilda_.boundaryField()[patchI];
            const fvPatchScalarField& muLamBound = muLam_.boundaryField()[patchI];
            const fvPatchScalarField& muTurbBound = muTurb_.boundaryField()[patchI];
            const vectorField  Cn = mesh_.boundary()[patchI].Cn();
            std::unique_ptr<Face> faceQuad;
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                scalar nuTildaFlux = 0.0;
                const label i = bfaceCells[faceI];
                build2ndFace(mesh_, start+faceI, faceQuad);
                for (label gaussI = 0; gaussI != faceQuad->size(); ++gaussI)
                {
                    const vector delta_i = faceQuad->at(gaussI) - Cn[faceI];
                    evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, nuTilda_i,
                                         rhoGrad_i, rhoUGrad_i, rhoEGrad_i, nuTildaGrad_i,
                                         rho_[i], rhoU_[i], rhoE_[i], nuTilda_[i], coefs_[i], coefsNuTilda_[i],
                                         basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateInletOutletFlux(rhoFlux, rhoUFlux, rhoEFlux, nuTildaFlux,
                                            rhoBound[faceI], UBound[faceI], pBound[faceI],
                                            nuTildaBound[faceI], muLamBound[faceI], muTurbBound[faceI],
                                            rho_i, rhoU_i, rhoE_i, rhoGrad_i, rhoUGrad_i, rhoEGrad_i, nuTildaGrad_i, 
                                            faceQuad->normal(gaussI), faceQuad->weight(gaussI));
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
                resNuTilda[i] -= nuTildaFlux;
            }
        }
    }
    forAll(mesh_.C(), cellI)
    {
        const Cell& quadCell = cellQuad_[cellI];
        for (label gaussI = 0; gaussI != quadCell.size(); ++gaussI)
        {
            const vector delta_i = quadCell.at(gaussI) - mesh_.C()[cellI];
            evaluateVarsAndGrads(rho_i, rhoU_i, rhoE_i, nuTilda_i,
                                 rhoGrad_i, rhoUGrad_i, rhoEGrad_i, nuTildaGrad_i,
                                 rho_[cellI], rhoU_[cellI], rhoE_[cellI], nuTilda_[cellI],
                                 coefs_[cellI], coefsNuTilda_[cellI], basisMean_[cellI], rDeltaXYZ_[cellI], delta_i);
            const tensor UGrad = (rho_i*rhoUGrad_i-rhoU_i*rhoGrad_i)/sqr(rho_i);
            const scalar Omega = mag(vorticity(UGrad));
            const scalar p = (rhoE_i-0.5*magSqr(rhoU_i)/rho_i)*(Gamma-1.0);
            const scalar T = Gamma*p/rho_i;
            const scalar nuLam = muLam(T, S_T_)/rho_i;
            resNuTilda[cellI] += evaluateSource(quadCell.distToWall(gaussI), Omega, nuLam, nuTilda_i,
                                                nuTildaGrad_i, quadCell.weight(gaussI));
        }
    }
}