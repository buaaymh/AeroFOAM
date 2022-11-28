void Foam::eulerPrimVar3rdSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    scalar rho_i, rho_j, p_i, p_j;
    vector U_i, U_j;
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
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
            rho_i = rho_[i] + (rhoGrad_[i]&delta_i);
            rho_j = rho_[j] + (rhoGrad_[j]&delta_j);
            U_i = U_[i] + (UGrad_[i]&delta_i);
            U_j = U_[j] + (UGrad_[j]&delta_j);
            p_i = p_[i] + (pGrad_[i]&delta_i);
            p_j = p_[j] + (pGrad_[j]&delta_j);
            evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                         rho_i, rho_j, U_i, U_j, p_i, p_j,
                         normal_[faceI], mesh_.magSf()[faceI]);
        }
        else
        {
            for (label gaussI = 0; gaussI != 4; ++gaussI)
            {
                const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
                const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
                evaluateVars(rho_i, U_i, p_i, rho_[i], U_[i], p_[i], rhoGrad_[i], UGrad_[i], pGrad_[i],
                             d2Rho_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], d2P_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                evaluateVars(rho_j, U_j, p_j, rho_[j], U_[j], p_[j], rhoGrad_[j], UGrad_[j], pGrad_[j],
                             d2Rho_[j], d2Ux_[j], d2Uy_[j], d2Uz_[j], d2P_[j], basisMean_[j], rDeltaXYZ_[j], delta_j);
                evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_i, rho_j, U_i, U_j, p_i, p_j,
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
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Ux_neigh = d2Ux_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uy_neigh = d2Uy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Uz_neigh = d2Uz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2P_neigh = d2P_.boundaryField()[patchI].patchNeighbourField();
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, j)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[j];
                gaussQuad4(mesh_, start+j, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - ownerCn[j];
                    const vector delta_j = quadPoints[gaussI] - neighCn[j];
                    evaluateVars(rho_i, U_i, p_i, rho_[i], U_[i], p_[i], rhoGrad_[i], UGrad_[i], pGrad_[i],
                                 d2Rho_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], d2P_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateVars(rho_j, U_j, p_j, rho_neigh[j], U_neigh[j], p_neigh[j], rhoGrad_neigh[j], UGrad_neigh[j], pGrad_neigh[j],
                                 d2Rho_neigh[j], d2Ux_neigh[j], d2Uy_neigh[j], d2Uz_neigh[j], d2P_neigh[j], basisMean_neigh[j], rDeltaXYZ_neigh[j], delta_j);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                 rho_i, rho_j, U_i, U_j, p_i, p_j,
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
                gaussQuad4(mesh_, start+faceI, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - Cn[faceI];
                    p_i   = p_[i]   + scalarPoly(delta_i, pGrad_[i], d2P_[i], rDeltaXYZ_[i], basisMean_[i]);
                    rhoUFlux += p_i*normal[faceI]*Quad4::w[gaussI];
                }
                resRhoU[i] -= rhoUFlux*magSf[faceI];
            }
            // forAll(bfaceCells, faceI)
            // {
            //     const label i = bfaceCells[faceI];
            //     scalar rhoFlux  = 0.0;
            //     vector rhoUFlux = vector::zero;
            //     scalar rhoEFlux = 0.0;
            //     gaussQuad9(mesh_, start+faceI, quadPoints);
            //     for (label gaussI = 0; gaussI != 9; ++gaussI)
            //     {
            //         const vector delta_i = quadPoints[gaussI] - Cn[faceI];
            //         evaluateVars(rho_i, U_i, p_i, rho_[i], U_[i], p_[i], rhoGrad_[i], UGrad_[i], pGrad_[i],
            //                      d2Rho_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], d2P_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
            //         vector N = vector(-quadPoints[gaussI].x(), -quadPoints[gaussI].y(), 0); N /= mag(N);
            //         U_j = vector(U_i[0]*(sqr(N[1])-sqr(N[0]))-2*N[0]*N[1]*U_i[1],
            //                      U_i[1]*(sqr(N[0])-sqr(N[1]))-2*N[0]*N[1]*U_i[0], 0);
            //         evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
            //                      rho_i, rho_i, U_i, U_j, p_i, p_i,
            //                      normal[faceI], magSf[faceI]*Quad9::w[gaussI]);
            //     }
            //     resRho[i]  -= rhoFlux;
            //     resRhoU[i] -= rhoUFlux;
            //     resRhoE[i] -= rhoEFlux;
            // }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            std::vector<vector> quadPoints(4, vector::zero);
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux  = 0.0;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux = 0.0;
                const label i = bfaceCells[faceI]; 
                gaussQuad4(mesh_, start+faceI, quadPoints);
                for (label gaussI = 0; gaussI != 4; ++gaussI)
                {
                    const vector delta_i = quadPoints[gaussI] - Cn[faceI];
                    evaluateVars(rho_i, U_i, p_i, rho_[i], U_[i], p_[i], rhoGrad_[i], UGrad_[i], pGrad_[i],
                                 d2Rho_[i], d2Ux_[i], d2Uy_[i], d2Uz_[i], d2P_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                 rho_i, rhoBound[faceI], U_i, UBound[faceI], p_i, pBound[faceI],
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
    const scalar& p_L,
    const scalar& p_R,
    const vector& normal,
    const scalar& magSf
) const
{
    scalar rhoFluxTemp;
    vector rhoUFluxTemp;
    scalar rhoEFluxTemp;
    riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
    rhoFlux  += magSf * rhoFluxTemp;
    rhoUFlux += magSf * rhoUFluxTemp;
    rhoEFlux += magSf * rhoEFluxTemp;
}
