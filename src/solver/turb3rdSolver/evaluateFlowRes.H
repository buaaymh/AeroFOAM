void Foam::turb3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    // reconstruction
    rhoGrad_ = fvc::grad(rho_);
    UGrad_   = fvc::grad(U_);
    TGrad_   = fvc::grad(T_);
    vr_.updateCoefficients(rho_, U_, T_, rhoGrad_, UGrad_, TGrad_,
                           d2Rho_, d2Ux_, d2Uy_, d2Uz_, d2T_, fluidProps_.gamma);
    // evaluate flux
    const volVectorField& rLengthScale = vr_.rLengthScale();
    const volSymmTensorField& basisConst = vr_.basisConst();
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux = 0.0;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux = 0.0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const scalar dxTilde = min(mesh_.V()[i], mesh_.V()[j])/mesh_.magSf()[faceI];
        for (label gaussI = 0; gaussI != 4; ++gaussI)
        {
            scalar rhoFluxTemp;
            vector rhoUFluxTemp = vector::zero;
            scalar rhoEFluxTemp;
            const vector delta_i = vr_.quadrature(faceI)[gaussI] - mesh_.C()[i];
            const vector delta_j = vr_.quadrature(faceI)[gaussI] - mesh_.C()[j];
            const scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGrad_[i], d2Rho_[i], rLengthScale[i], basisConst[i]);
            const scalar rho_j = rho_[j] + scalarPoly(delta_j, rhoGrad_[j], d2Rho_[j], rLengthScale[j], basisConst[j]);
            const vector U_i = U_[i] + vectorPoly(delta_i, UGrad_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], rLengthScale[i], basisConst[i]);
            const vector U_j = U_[j] + vectorPoly(delta_j, UGrad_[j], d2Ux_[j], d2Uy_[j], d2Uz_[j], rLengthScale[j], basisConst[j]);
            const scalar T_i = T_[i] + scalarPoly(delta_i, TGrad_[i], d2T_[i], rLengthScale[i], basisConst[i]);
            const scalar T_j = T_[j] + scalarPoly(delta_j, TGrad_[j], d2T_[j], rLengthScale[j], basisConst[j]);
            const tensor UGrad_i = vectorPolyGrad(delta_i, UGrad_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], rLengthScale[i]);
            const tensor UGrad_j = vectorPolyGrad(delta_j, UGrad_[j], d2Ux_[j], d2Uy_[j], d2Uz_[j], rLengthScale[j]);
            const vector TGrad_i = scalarPolyGrad(delta_i, TGrad_[i], d2T_[i], rLengthScale[i]);
            const vector TGrad_j = scalarPolyGrad(delta_j, TGrad_[j], d2T_[j], rLengthScale[j]);
            evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                         rho_i, rho_j, U_i, U_j, T_i, T_j, UGrad_i, UGrad_j, TGrad_i, TGrad_j, nuTilde_[i], nuTilde_[j],
                         dxTilde , normal_[faceI]);
            rhoFlux  += rhoFluxTemp  * Quad4::w[gaussI];
            rhoUFlux += rhoUFluxTemp * Quad4::w[gaussI];
            rhoEFlux += rhoEFluxTemp * Quad4::w[gaussI];
        }
        const scalar temp = mesh_.magSf()[faceI] * 0.25;
        rhoFlux  *= temp;
        rhoUFlux *= temp;
        rhoEFlux *= temp;
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
        const label faceI = mesh_.boundary()[patchI].start();
        if (isA<processorFvPatch>(mesh_.boundary()[patchI]))
        {
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const vectorField neighCn = mesh_.C().boundaryField()[patchI].patchNeighbourField();
            const vectorField rLengthScale_neigh = rLengthScale.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisConst_neigh = basisConst.boundaryField()[patchI].patchNeighbourField();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh   = T_.boundaryField()[patchI].patchNeighbourField();
            const scalarField nuTilde_neigh = nuTilde_.boundaryField()[patchI].patchNeighbourField();
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
                vr_.gaussQuad4(faceI+j, quadPoints);
                const scalar dxTilde = mesh_.V()[i]/magSf[j];
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    scalar rhoFluxTemp;
                    vector rhoUFluxTemp = vector::zero;
                    scalar rhoEFluxTemp;
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];

                    const scalar rho_i = rho_[i]      + scalarPoly(delta_i, rhoGrad_[i],      d2Rho_[i],      rLengthScale[i],       basisConst[i]);
                    const scalar rho_j = rho_neigh[j] + scalarPoly(delta_j, rhoGrad_neigh[j], d2Rho_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    const vector U_i = U_[i]      + vectorPoly(delta_i, UGrad_[i],      d2Ux_[i],      d2Uy_[i],      d2Uz_[i],      rLengthScale[i],       basisConst[i]);
                    const vector U_j = U_neigh[j] + vectorPoly(delta_j, UGrad_neigh[j], d2Ux_neigh[j], d2Uy_neigh[j], d2Uz_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    const scalar T_i = T_[i]      + scalarPoly(delta_i, TGrad_[i],      d2T_[i],      rLengthScale[i],       basisConst[i]);
                    const scalar T_j = T_neigh[j] + scalarPoly(delta_j, TGrad_neigh[j], d2T_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    const tensor UGrad_i = vectorPolyGrad(delta_i, UGrad_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], rLengthScale[i]);
                    const tensor UGrad_j = vectorPolyGrad(delta_j, UGrad_neigh[j], d2Ux_neigh[j], d2Uy_neigh[j], d2Uz_neigh[j], rLengthScale_neigh[j]);
                    const vector TGrad_i = scalarPolyGrad(delta_i, TGrad_[i], d2T_[i], rLengthScale[i]);
                    const vector TGrad_j = scalarPolyGrad(delta_j, TGrad_neigh[j], d2T_neigh[j], rLengthScale_neigh[j]);
            
                    evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                 rho_i, rho_j, U_i, U_j, T_i, T_j, UGrad_i, UGrad_j, TGrad_i, TGrad_j, nuTilde_[i], nuTilde_neigh[j],
                                 dxTilde , normal[j]);
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
        if (type == "symmetryPlane")
        {
            const scalarField& p_patch = p_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux = p_patch[faceI]*normal[faceI];
                const label i = bfaceCells[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
            }
        }
        if (type == "wall")
        {
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const fvPatchScalarField& TBound = T_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                vector rhoUFlux = vector::zero;
                const label i = bfaceCells[faceI];
                evaluateWallFlux(rhoUFlux, pBound[faceI], TBound[faceI], UGrad_[i], normal[faceI]);
                resRhoU[i] -= rhoUFlux * magSf[faceI];              
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const fvPatchScalarField& nuTildeBound = nuTilde_.boundaryField()[patchI];
            const scalarField nuTildeInter = nuTildeBound.patchInternalField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar nuTilde = 0.5*(nuTildeBound[faceI]+nuTildeInter[faceI]);
                evaluateInletOutletFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                        rhoBound[faceI], UBound[faceI], pBound[faceI], nuTilde,
                                        UGrad_[i], TGrad_[i], normal[faceI]);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
    }
}
