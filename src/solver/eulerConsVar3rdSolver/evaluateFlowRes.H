void Foam::eulerConsVar3rdSolver::evaluateFlowRes
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
            rhoU_i = rhoU_[i] + (rhoUGrad_[i]&delta_i);
            rhoU_j = rhoU_[j] + (rhoUGrad_[j]&delta_j);
            rhoE_i = rhoE_[i] + (rhoEGrad_[i]&delta_i);
            rhoE_j = rhoE_[j] + (rhoEGrad_[j]&delta_j);
            evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                         rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j,
                         normal_[faceI], mesh_.magSf()[faceI]);
        }
        else
        {
            for (label gaussI = 0; gaussI != 4; ++gaussI)
            {
                const vector delta_i = quad_[faceI][gaussI] - mesh_.C()[i];
                const vector delta_j = quad_[faceI][gaussI] - mesh_.C()[j];
                evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], rhoGrad_[i], rhoUGrad_[i], rhoEGrad_[i],
                             d2Rho_[i], d2RhoUx_[i], d2RhoUy_[i], d2RhoUz_[i], d2RhoE_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                evaluateVars(rho_j, rhoU_j, rhoE_j, rho_[j], rhoU_[j], rhoE_[j], rhoGrad_[j], rhoUGrad_[j], rhoEGrad_[j],
                             d2Rho_[j], d2RhoUx_[j], d2RhoUy_[j], d2RhoUz_[j], d2RhoE_[j], basisMean_[j], rDeltaXYZ_[j], delta_j);
                evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                             rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j,
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
            const scalarField rho_neigh  = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoU_neigh = rhoU_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoE_neigh = rhoE_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh  = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField rhoUGrad_neigh = rhoUGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoEGrad_neigh = rhoEGrad_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUx_neigh = d2RhoUx_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUy_neigh = d2RhoUy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUz_neigh = d2RhoUz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoE_neigh = d2RhoE_.boundaryField()[patchI].patchNeighbourField();
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
                    evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], rhoGrad_[i], rhoUGrad_[i], rhoEGrad_[i],
                                 d2Rho_[i], d2RhoUx_[i], d2RhoUy_[i], d2RhoUz_[i], d2RhoE_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateVars(rho_j, rhoU_j, rhoE_j, rho_neigh[j], rhoU_neigh[j], rhoE_neigh[j], rhoGrad_neigh[j], rhoUGrad_neigh[j], rhoEGrad_neigh[j],
                                 d2Rho_neigh[j], d2RhoUx_neigh[j], d2RhoUy_neigh[j], d2RhoUz_neigh[j], d2RhoE_neigh[j], basisMean_neigh[j], rDeltaXYZ_neigh[j], delta_j);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                 rho_i, rho_j, rhoU_i, rhoU_j, rhoE_i, rhoE_j,
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
                    evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], rhoGrad_[i], rhoUGrad_[i], rhoEGrad_[i],
                                 d2Rho_[i], d2RhoUx_[i], d2RhoUy_[i], d2RhoUz_[i], d2RhoE_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    rhoUFlux += (rhoE_i-0.5*magSqr(rhoU_i)/rho_i)*(fluidProps_.gamma-1.0)*normal[faceI]*Quad4::w[gaussI];
                }
                resRhoU[i] -= rhoUFlux*magSf[faceI];
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            const fvPatchScalarField& rhoBound  = rho_.boundaryField()[patchI];
            const fvPatchVectorField& rhoUBound = rhoU_.boundaryField()[patchI];
            const fvPatchScalarField& rhoEBound = rhoE_.boundaryField()[patchI];
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
                    evaluateVars(rho_i, rhoU_i, rhoE_i, rho_[i], rhoU_[i], rhoE_[i], rhoGrad_[i], rhoUGrad_[i], rhoEGrad_[i],
                                 d2Rho_[i], d2RhoUx_[i], d2RhoUy_[i], d2RhoUz_[i], d2RhoE_[i], basisMean_[i], rDeltaXYZ_[i], delta_i);
                    evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                 rho_i, rhoBound[faceI], rhoU_i, rhoUBound[faceI], rhoE_i, rhoEBound[faceI],
                                 normal[faceI], magSf[faceI]*Quad4::w[gaussI]);
                }
                resRho[i]  -= rhoFlux;
                resRhoU[i] -= rhoUFlux;
                resRhoE[i] -= rhoEFlux;
            }
        }
    }
}

void Foam::eulerConsVar3rdSolver::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar& rho_L,
    const scalar& rho_R,
    const vector& rhoU_L,
    const vector& rhoU_R,
    const scalar& rhoE_L,
    const scalar& rhoE_R,
    const vector& normal,
    const scalar& magSf
) const
{
    scalar rhoFluxTemp;
    vector rhoUFluxTemp;
    scalar rhoEFluxTemp;

    scalar p_L, p_R, T_L, T_R;
    vector U_L, U_R;
    consToPrim(rho_L, rhoU_L, rhoE_L, fluidProps_.gamma, U_L, p_L, T_L);
    consToPrim(rho_R, rhoU_R, rhoE_R, fluidProps_.gamma, U_R, p_R, T_R);

    riemann_->evaluateFlux(rhoFluxTemp, rhoUFluxTemp, rhoEFluxTemp,
                           rho_L, rho_R, U_L, U_R, p_L, p_R,
                           normal, fluidProps_.gamma);
    rhoFlux  += magSf * rhoFluxTemp;
    rhoUFlux += magSf * rhoUFluxTemp;
    rhoEFlux += magSf * rhoEFluxTemp;
}
