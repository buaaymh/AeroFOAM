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
    vr_.updateCoefficients(rho_, U_, T_, rhoGrad_, UGrad_, TGrad_,
                           d2Rho_, d2Ux_, d2Uy_, d2Uz_, d2T_);
    const volVectorField& rLengthScale = vr_.rLengthScale();
    const volSymmTensorField& basisConst = vr_.basisConst();
    const std::vector<bool>& p0 = vr_.p0();
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux = 0.0;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux = 0.0;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        if (p0[i] && p0[j])
        {
            const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
            const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
            const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i);
            const scalar rho_j = rho_[j] + (rhoGrad_[j]&delta_j);
            const vector U_i = U_[i] + (UGrad_[i]&delta_i);
            const vector U_j = U_[j] + (UGrad_[j]&delta_j);
            const scalar T_i = T_[i] + (TGrad_[i]&delta_i);
            const scalar T_j = T_[j] + (TGrad_[j]&delta_j);
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
                const vector delta_i = vr_.quadrature(faceI)[gaussI] - mesh_.C()[i];
                const vector delta_j = vr_.quadrature(faceI)[gaussI] - mesh_.C()[j];
                const scalar rho_i = rho_[i] + scalarPoly(delta_i, rhoGrad_[i], d2Rho_[i], rLengthScale[i], basisConst[i]);
                const scalar rho_j = rho_[j] + scalarPoly(delta_j, rhoGrad_[j], d2Rho_[j], rLengthScale[j], basisConst[j]);
                const vector U_i = U_[i] + vectorPoly(delta_i, UGrad_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], rLengthScale[i], basisConst[i]);
                const vector U_j = U_[j] + vectorPoly(delta_j, UGrad_[j], d2Ux_[j], d2Uy_[j], d2Uz_[j], rLengthScale[j], basisConst[j]);
                const scalar T_i = T_[i] + scalarPoly(delta_i, TGrad_[i], d2T_[i], rLengthScale[i], basisConst[i]);
                const scalar T_j = T_[j] + scalarPoly(delta_j, TGrad_[j], d2T_[j], rLengthScale[j], basisConst[j]);
                riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                        rho_i, rho_j, U_i, U_j, T_i, T_j,
                                        normal_[faceI], fluidProps_.gamma);
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
        const label faceI = mesh_.boundary()[patchI].start();
        if (isA<cyclicFvPatch>(mesh_.boundary()[patchI]))
        {
            const cyclicFvPatch& ownerPatch = refCast<const cyclicFvPatch>(mesh_.boundary()[patchI]);
            const cyclicFvPatch& neighPatch = ownerPatch.neighbFvPatch();
            const vectorField ownerCn = ownerPatch.Cn();
            const vectorField neighCn = neighPatch.Cn();
            const vectorField ownToNei = neighPatch.Cf() - ownerPatch.Cf();
            const vectorField rLengthScale_neigh = rLengthScale.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisConst_neigh = basisConst.boundaryField()[patchI].patchNeighbourField();
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
                vr_.gaussQuad4(faceI+j, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    scalar rhoFluxTemp;
                    vector rhoUFluxTemp = vector::zero;
                    scalar rhoEFluxTemp;
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - (neighCn[j]-ownToNei[j]);
                    const scalar rho_i = rho_[i]      + scalarPoly(delta_i, rhoGrad_[i],      d2Rho_[i],      rLengthScale[i],       basisConst[i]);
                    const scalar rho_j = rho_neigh[j] + scalarPoly(delta_j, rhoGrad_neigh[j], d2Rho_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    const vector U_i = U_[i]      + vectorPoly(delta_i, UGrad_[i],      d2Ux_[i],      d2Uy_[i],      d2Uz_[i],      rLengthScale[i],       basisConst[i]);
                    const vector U_j = U_neigh[j] + vectorPoly(delta_j, UGrad_neigh[j], d2Ux_neigh[j], d2Uy_neigh[j], d2Uz_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    const scalar T_i = T_[i]      + scalarPoly(delta_i, TGrad_[i],      d2T_[i],      rLengthScale[i],       basisConst[i]);
                    const scalar T_j = T_neigh[j] + scalarPoly(delta_j, TGrad_neigh[j], d2T_neigh[j], rLengthScale_neigh[j], basisConst_neigh[j]);
                    riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                           rho_i, rho_j, U_i, U_j, T_i, T_j,
                                           normal[j], fluidProps_.gamma);
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
        if (isA<processorFvPatch>(mesh_.boundary()[patchI]))
        {
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const vectorField neighCn = mesh_.C().boundaryField()[patchI].patchNeighbourField();
            const vectorField rLengthScale_neigh = rLengthScale.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisConst_neigh = basisConst.boundaryField()[patchI].patchNeighbourField();
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
                vr_.gaussQuad4(faceI+j, quadPoints);
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
                    riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                                           rho_i, rho_j, U_i, U_j, T_i, T_j,
                                           normal[j], fluidProps_.gamma);
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
        if (name == "inlet")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const fvPatchScalarField& Ma = Ma_.boundaryField()[patchI];
            const scalarField rhoInter = rhoBound.patchInternalField();
            const vectorField UInter = UBound.patchInternalField();
            const scalarField pInter = pBound.patchInternalField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                scalar rho = rhoBound[faceI];
                vector U = UBound[faceI];
                scalar p = pBound[faceI];
                const label i = bfaceCells[faceI]; 
                if (Ma[faceI] < 1)
                {
                    const scalar rhoc = rhoInter[faceI]*c_[i];
                    p = 0.5*(p+pInter[faceI]-rhoc*(normal[faceI]&(U-UInter[faceI])));
                    rho += (p-pBound[faceI])/sqr(c_[i]);
                    U += normal[faceI]*(p-pBound[faceI])/rhoc;
                }
                Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                   rho, U, p,
                                   normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
        if (name == "outlet")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const scalarField pInter = pBound.patchInternalField();
            const scalarField Ma = Ma_.boundaryField()[patchI].patchInternalField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                scalar rho = rhoBound[faceI];
                vector U = UBound[faceI];
                scalar p = pBound[faceI];
                const label i = bfaceCells[faceI];
                if (Ma[faceI] < 1.0)
                {
                    rho -= (pInter[faceI]-p)/sqr(c_[i]);
                    U += normal[faceI]*(pInter[faceI]-p)/(rho*c_[i]);
                }
                Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                    rho, U, p,
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
            const scalarField& Ma = Ma_.boundaryField()[patchI];
            const scalarField& Vn = normal&UBound;
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                if (Ma[faceI] < 1.0)
                {
                    const scalar c2 = pInter[faceI]*fluidProps_.gamma/rhoInter[faceI];
                    const scalar rhoc = rhoInter[faceI]*sqrt(c2);
                    scalar rho, p; vector U;
                    if (Vn[faceI] < 0)
                    {
                        p = 0.5*(pBound[faceI]+pInter[faceI]-rhoc*(normal[faceI]&(UBound[faceI]-UInter[faceI])));
                        rho = rhoBound[faceI] + (p-pBound[faceI])/c2;
                        U = UBound[faceI] + normal[faceI]*(p-pBound[faceI])/rhoc;
                    }
                    else
                    {
                        p = pBound[faceI];
                        rho = rhoInter[faceI] + (p-pInter[faceI])/c2;
                        U = UInter[faceI] - normal[faceI]*(p-pInter[faceI])/rhoc;
                    }
                    Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho, U, p,
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
    }
}
