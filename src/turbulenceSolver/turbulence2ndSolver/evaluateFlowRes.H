void Foam::turbulence2ndSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    rhoGrad_ = fvc::grad(rho_);
    pGrad_   = fvc::grad(p_);
    TGrad_   = fvc::grad(T_);
    const word limiter = mesh_.schemesDict().subDict("gradSchemes").lookupOrDefault<word>("Limiter", "Venkat");
    if (limiter == "Venkat") limiterVenkatakrishnan(rho_, U_, p_, rhoGrad_, UGrad_, pGrad_);
    else if (limiter == "Barth") limiterBarthJespersen(rho_, U_, p_, rhoGrad_, UGrad_, pGrad_);
    else Info << "Error in limiter type!" << endl;
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux;
        vector rhoUFlux;
        scalar rhoEFlux;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i)*rhoLimit_[i];
        const scalar rho_j = rho_[j] + (rhoGrad_[j]&delta_j)*rhoLimit_[j];
        const vector U_i = U_[i] + cmptMultiply(UGrad_[i]&delta_i, ULimit_[i]);
        const vector U_j = U_[j] + cmptMultiply(UGrad_[j]&delta_j, ULimit_[j]);
        const scalar p_i = p_[i] + (pGrad_[i]&delta_i)*pLimit_[i];
        const scalar p_j = p_[j] + (pGrad_[j]&delta_j)*pLimit_[j];
        riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                               rho_i, rho_j, U_i, U_j, p_i, p_j, normal_[faceI]);
        resRho[i]  -= rhoFlux * mesh_.magSf()[faceI];
        resRhoU[i] -= rhoUFlux* mesh_.magSf()[faceI];
        resRhoE[i] -= rhoEFlux* mesh_.magSf()[faceI];
        resRho[j]  += rhoFlux * mesh_.magSf()[faceI];
        resRhoU[j] += rhoUFlux* mesh_.magSf()[faceI];
        resRhoE[j] += rhoEFlux* mesh_.magSf()[faceI];
        scalar dx_L  = mag(delta_i&normal_[faceI]) + SMALL;
        scalar dx_R  = mag(delta_j&normal_[faceI]) + SMALL;
        scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx;
        evaluateViscousFlux(rhoUFlux, rhoEFlux,
                            dx_L*U_[i]     +dx_R*U_[j],
                            dx_L*TGrad_[i] +dx_R*TGrad_[j],
                            dx_L*UGrad_[i] +dx_R*UGrad_[j],
                            dx_L*muLam_[i] +dx_R*muLam_[j],
                            dx_L*muTurb_[i]+dx_R*muTurb_[j], normal_[faceI]);
        resRhoU[i] += rhoUFlux*mesh_.magSf()[faceI];
        resRhoE[i] += rhoEFlux*mesh_.magSf()[faceI];
        resRhoU[j] -= rhoUFlux*mesh_.magSf()[faceI];
        resRhoE[j] -= rhoEFlux*mesh_.magSf()[faceI];
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
            const vectorField delta = mesh_.boundary()[patchI].delta();
            const vectorField delta_i = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            const vectorField delta_j = delta_i - mesh_.boundary()[patchI].delta();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh = T_.boundaryField()[patchI].patchNeighbourField();
            const scalarField muLam_neigh = muLam_.boundaryField()[patchI].patchNeighbourField();
            const scalarField muTurb_neigh= muTurb_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad_.boundaryField()[patchI].patchNeighbourField();
            const vectorField TGrad_neigh   = TGrad_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoLimit_neigh = rhoLimit_.boundaryField()[patchI].patchNeighbourField();
            const vectorField ULimit_neigh   = ULimit_.boundaryField()[patchI].patchNeighbourField();
            const scalarField pLimit_neigh   = pLimit_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar rho_i = rho_[i] + (rhoGrad_[i]&delta_i[faceI])*rhoLimit_[i];
                const scalar rho_j = rho_neigh[faceI] + (rhoGrad_neigh[faceI]&delta_j[faceI])*rhoLimit_neigh[faceI];
                const vector U_i = U_[i] + cmptMultiply(UGrad_[i]&delta_i[faceI], ULimit_[i]);
                const vector U_j = U_neigh[faceI] + cmptMultiply(UGrad_neigh[faceI]&delta_j[faceI], ULimit_neigh[faceI]);
                const scalar p_i = p_[i] + (pGrad_[i]&delta_i[faceI])*pLimit_[i];
                const scalar p_j = p_neigh[faceI] + (pGrad_neigh[faceI]&delta_j[faceI])*pLimit_neigh[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, p_i, p_j, normal[faceI]);
                resRho[i]  -= rhoFlux *magSf[faceI];
                resRhoU[i] -= rhoUFlux*magSf[faceI];
                resRhoE[i] -= rhoEFlux*magSf[faceI];
                scalar dx_L  = mag(delta_i[faceI]&normal[faceI]) + SMALL;
                scalar dx_R  = mag(delta_j[faceI]&normal[faceI]) + SMALL;
                scalar dx = dx_L + dx_R; dx_L /= dx; dx_R /= dx; 
                evaluateViscousFlux(rhoUFlux, rhoEFlux,
                                    dx_L*U_[i]     +dx_R*U_neigh[faceI],
                                    dx_L*TGrad_[i] +dx_R*TGrad_neigh[faceI],
                                    dx_L*UGrad_[i] +dx_R*UGrad_neigh[faceI],
                                    dx_L*muLam_[i] +dx_R*muLam_neigh[faceI],
                                    dx_L*muTurb_[i]+dx_R*muTurb_neigh[faceI], normal[faceI]);
                resRhoU[i] += rhoUFlux*magSf[faceI];
                resRhoE[i] += rhoEFlux*magSf[faceI];
            }
        }
        if (type == "symmetryPlane" || type == "symmetry")
        {
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                vector rhoUFlux = pBound[faceI]*normal[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
            }
        }
        if (type == "wall")
        {
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const fvPatchScalarField& muLamWall = muLam_.boundaryField()[patchI];
            const tensorField UGradWall = normal*U_.boundaryField()[patchI].snGrad();
            forAll(bfaceCells, faceI)
            {
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                vector rhoUFlux = pBound[faceI]*normal[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                evaluateViscousFlux(rhoUFlux, rhoEFlux, vector::zero, vector::zero, UGradWall[faceI],
                                    muLamWall[faceI], 0, normal[faceI]);
                resRhoU[i] += rhoUFlux * magSf[faceI];
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho_.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            const fvPatchScalarField& muLamBound = muLam_.boundaryField()[patchI];
            const fvPatchScalarField& muTurbBound = muTurb_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                Foam::evaluateConvectiveFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                             rhoBound[faceI], UBound[faceI], pBound[faceI],
                                             normal[faceI]);
                resRho[i]  -= rhoFlux  * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
                evaluateViscousFlux(rhoUFlux, rhoEFlux, UBound[faceI], TGrad_[i], UGrad_[i],
                                    muLamBound[faceI], muTurbBound[faceI], normal[faceI]);
                resRhoU[i] += rhoUFlux * magSf[faceI];
                resRhoE[i] += rhoEFlux * magSf[faceI];
            }
        }
    }
}
