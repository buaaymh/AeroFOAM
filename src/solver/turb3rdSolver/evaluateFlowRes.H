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
    rhoGrad_ = fvc::grad(rho_);
    UGrad_   = fvc::grad(U_);
    pGrad_   = fvc::grad(p_);
    updateCoefficients();
    limitCoefficients();
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
            const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
            const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
            scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i], basisConst_[i]);
            scalar rho_j = rho_[j] + scalarPoly(delta_j, rhoGradLimited_[j], d2RhoLimited_[j], rLengthScale_[j], basisConst_[j]);
            vector U_i = U_[i] + vectorPoly(delta_i, UGradLimited_[i], d2UxLimited_[i], d2UyLimited_[i], d2UzLimited_[i], rLengthScale_[i], basisConst_[i]);
            vector U_j = U_[j] + vectorPoly(delta_j, UGradLimited_[j], d2UxLimited_[j], d2UyLimited_[j], d2UzLimited_[j], rLengthScale_[j], basisConst_[j]);
            scalar p_i = p_[i] + scalarPoly(delta_i, pGradLimited_[i], d2PLimited_[i], rLengthScale_[i], basisConst_[i]);
            scalar p_j = p_[j] + scalarPoly(delta_j, pGradLimited_[j], d2PLimited_[j], rLengthScale_[j], basisConst_[j]);
            const vector rhoGrad_i = scalarPolyGrad(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i]);
            const vector rhoGrad_j = scalarPolyGrad(delta_j, rhoGradLimited_[j], d2RhoLimited_[j], rLengthScale_[j]);
            const tensor UGrad_i = vectorPolyGrad(delta_i, UGradLimited_[i], d2UxLimited_[i], d2UyLimited_[i], d2UzLimited_[i], rLengthScale_[i]);
            const tensor UGrad_j = vectorPolyGrad(delta_j, UGradLimited_[j], d2UxLimited_[j], d2UyLimited_[j], d2UzLimited_[j], rLengthScale_[j]);
            const vector pGrad_i = scalarPolyGrad(delta_i, pGradLimited_[i], d2PLimited_[i], rLengthScale_[i]);
            const vector pGrad_j = scalarPolyGrad(delta_j, pGradLimited_[j], d2PLimited_[j], rLengthScale_[j]);
            evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                         rho_i, rho_j, U_i, U_j, p_i, p_j,
                         rhoGrad_i, rhoGrad_j, UGrad_i, UGrad_j, pGrad_i, pGrad_j, nuTilde_[i], nuTilde_[j],
                         dxTilde, normal_[faceI]);
            rhoFlux  += rhoFluxTemp  * Quad4::w[gaussI];
            rhoUFlux += rhoUFluxTemp * Quad4::w[gaussI];
            rhoEFlux += rhoEFluxTemp * Quad4::w[gaussI];
        }
        const scalar temp = 0.25 * mesh_.magSf()[faceI];
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
        const label start = mesh_.boundary()[patchI].start();
        if (isA<processorFvPatch>(mesh_.boundary()[patchI]))
        {
            
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const vectorField neighCn = mesh_.C().boundaryField()[patchI].patchNeighbourField();
            const vectorField rLengthScale_neigh = rLengthScale_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisConst_neigh = basisConst_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            const scalarField nuTilde_neigh = nuTilde_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGradLimited_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGradLimited_.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGradLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2RhoLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Ux_neigh = d2UxLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uy_neigh = d2UyLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uz_neigh = d2UzLimited_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2P_neigh = d2PLimited_.boundaryField()[patchI].patchNeighbourField();
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, j)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[j];
                gaussQuad4(start+j, quadPoints);
                const scalar dxTilde = mesh_.V()[i]/magSf[j];
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
                    const scalar p_i = p_[i]      + scalarPoly(delta_i, pGradLimited_[i], d2PLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar p_j = p_neigh[j] + scalarPoly(delta_j, pGrad_neigh[j], d2P_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    
                    const vector rhoGrad_i = scalarPolyGrad(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i]);
                    const vector rhoGrad_j = scalarPolyGrad(delta_j, rhoGrad_neigh[j], d2Rho_neigh[j], rLengthScale_neigh[j]);
                    const tensor UGrad_i = vectorPolyGrad(delta_i, UGradLimited_[i], d2UxLimited_[i], d2UyLimited_[i], d2UzLimited_[i], rLengthScale_[i]);
                    const tensor UGrad_j = vectorPolyGrad(delta_j, UGrad_neigh[j], d2Ux_neigh[j], d2Uy_neigh[j], d2Uz_neigh[j], rLengthScale_neigh[j]);
                    const vector pGrad_i = scalarPolyGrad(delta_i, pGradLimited_[i], d2PLimited_[i], rLengthScale_[i]);
                    const vector pGrad_j = scalarPolyGrad(delta_j, pGrad_neigh[j], d2P_neigh[j], rLengthScale_neigh[j]);
                    evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                 rho_i, rho_j, U_i, U_j, p_i, p_j, rhoGrad_i, rhoGrad_j, UGrad_i, UGrad_j, pGrad_i, pGrad_j, nuTilde_[i], nuTilde_neigh[j],
                                 dxTilde, normal[j]);
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
        if (type == "symmetryPlane" || type == "symmetry")
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
                    const scalar p_i = p_[i] + scalarPoly(delta_i, pGradLimited_[i], d2PLimited_[i], rLengthScale_[i], basisConst_[i]);
                    rhoUFlux += p_i * normal[faceI] * Quad4::w[gaussI];
                }
                resRhoU[i] -= 0.25 * rhoUFlux * magSf[faceI];
            }
        }
        if (type == "wall")
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
                    vector rhoUFluxTemp = vector::zero;
                    const vector delta_i = quadPoints[gaussI] - Cn[faceI];
                    const scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar p_i = p_[i] + scalarPoly(delta_i, pGradLimited_[i], d2PLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const tensor UGrad_i = vectorPolyGrad(delta_i, UGrad_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], rLengthScale_[i]);
                    evaluateWallFlux(rhoUFluxTemp, rho_i, p_i, UGrad_i, normal[faceI]);
                    rhoUFlux += rhoUFluxTemp * Quad4::w[gaussI];
                }
                resRhoU[i] -= 0.25 * rhoUFlux * magSf[faceI];
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const fvPatchScalarField& nuTildeBound = nuTilde_.boundaryField()[patchI];
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                gaussQuad4(start+faceI, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    scalar rhoFluxTemp;
                    vector rhoUFluxTemp = vector::zero;
                    scalar rhoEFluxTemp;
                    const vector delta_i = quadPoints[gaussI] - Cn[faceI];
                    const scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGradLimited_[i], d2RhoLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const vector U_i = U_[i] + vectorPoly(delta_i, UGradLimited_[i], d2UxLimited_[i], d2UyLimited_[i], d2UzLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar p_i = p_[i] + scalarPoly(delta_i, pGradLimited_[i], d2PLimited_[i], rLengthScale_[i], basisConst_[i]);
                    const scalar nuTilde = 0.5*(nuTildeBound[faceI]+nuTilde_[i]);
                    const vector rhoGrad_i = scalarPolyGrad(delta_i, rhoGrad_[i], d2Rho_[i], rLengthScale_[i]);
                    const tensor UGrad_i = vectorPolyGrad(delta_i, UGrad_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], rLengthScale_[i]);
                    const vector pGrad_i = scalarPolyGrad(delta_i, pGrad_[i], d2P_[i], rLengthScale_[i]);
                    evaluateInletOutletFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                            rho_i, rhoBound[faceI], U_i, UBound[faceI], p_i, pBound[faceI], nuTilde,
                                            rhoGrad_i, UGrad_i, pGrad_i, normal[faceI]);
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
