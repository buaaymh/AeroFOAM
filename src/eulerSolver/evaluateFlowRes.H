void Foam::eulerSolver::evaluateFlowResForGMRES
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
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
        riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                               rho_[i], rho_[j], U_[i], U_[j], p_[i], p_[j],
                               normal_[faceI]);
        resRho[i]  -= rhoFlux*mesh_.magSf()[faceI];
        resRhoU[i] -= rhoUFlux*mesh_.magSf()[faceI];
        resRhoE[i] -= rhoEFlux*mesh_.magSf()[faceI];
        resRho[j]  += rhoFlux*mesh_.magSf()[faceI];
        resRhoU[j] += rhoUFlux*mesh_.magSf()[faceI];
        resRhoE[j] += rhoEFlux*mesh_.magSf()[faceI];
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
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                        rho_[i], rho_neigh[faceI], U_[i], U_neigh[faceI], p_[i], p_neigh[faceI],
                                        normal[faceI]);
                resRho[i]  -= rhoFlux*magSf[faceI];
                resRhoU[i] -= rhoUFlux*magSf[faceI];
                resRhoE[i] -= rhoEFlux*magSf[faceI];
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const fvPatchScalarField& pBound = p_.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                vector rhoUFlux = pBound[faceI]*normal[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
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
                                            normal[faceI]);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
    }
}
