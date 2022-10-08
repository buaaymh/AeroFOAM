void Foam::euler2ndSolver::evaluateFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{   
    rhoGrad_ = fvc::grad(rho_);
    UGrad_   = fvc::grad(U_);
    pGrad_   = fvc::grad(p_);
    evaluateFlowRes(rho_, U_, p_, rhoGrad_, UGrad_, pGrad_, resRho, resRhoU, resRhoE);
}

void Foam::euler2ndSolver::evaluateFlowRes
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p,
    const volVectorField& rhoGrad,
    const volTensorField& UGrad,
    const volVectorField& pGrad,
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    volScalarField rhoLimit
    (
        IOobject
        (
            "rhoLimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    );
    volVectorField ULimit
    (
        IOobject
        (
            "ULimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless, vector::one)
    );
    volScalarField pLimit
    (
        IOobject
        (
            "pLimit",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 1)
    );
    limitGrad(rho, U, p, rhoGrad, UGrad, pGrad, rhoLimit, ULimit, pLimit);
    resRho  = scalarField(mesh_.nCells(), 0.0);
    resRhoU = vectorField(mesh_.nCells(), vector::zero);
    resRhoE = scalarField(mesh_.nCells(), 0.0);
    forAll(mesh_.owner(), faceI)
    {
        scalar rhoFlux;
        vector rhoUFlux = vector::zero;
        scalar rhoEFlux;
        const label i = mesh_.owner()[faceI];
        const label j = mesh_.neighbour()[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        const scalar rho_i = rho[i] + (rhoGrad[i]&delta_i)*rhoLimit[i];
        const scalar rho_j = rho[j] + (rhoGrad[j]&delta_j)*rhoLimit[j];
        const vector U_i = U[i] + cmptMultiply(UGrad[i]&delta_i, ULimit[i]);
        const vector U_j = U[j] + cmptMultiply(UGrad[j]&delta_j, ULimit[j]);
        const scalar p_i = p[i] + (pGrad_[i]&delta_i)*pLimit[i];
        const scalar p_j = p[j] + (pGrad_[j]&delta_j)*pLimit[j];
        riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                               rho_i, rho_j, U_i, U_j, p_i, p_j,
                               normal_[faceI], fluidProps_.gamma);
        resRho[i]  -= rhoFlux  * mesh_.magSf()[faceI];
        resRhoU[i] -= rhoUFlux * mesh_.magSf()[faceI];
        resRhoE[i] -= rhoEFlux * mesh_.magSf()[faceI];
        resRho[j]  += rhoFlux  * mesh_.magSf()[faceI];
        resRhoU[j] += rhoUFlux * mesh_.magSf()[faceI];
        resRhoE[j] += rhoEFlux * mesh_.magSf()[faceI];
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& normal = normal_.boundaryField()[patchI];
        const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
        const word name = mesh_.boundary()[patchI].name();
        const word type = mesh_.boundary()[patchI].type();
        if (isA<cyclicFvPatch>(mesh_.boundary()[patchI]))
        {
            const cyclicFvPatch& ownerPatch = refCast<const cyclicFvPatch>(mesh_.boundary()[patchI]);
            const cyclicFvPatch& neighPatch = ownerPatch.neighbFvPatch();
            const vectorField delta_i = ownerPatch.Cf() - ownerPatch.Cn();
            const vectorField delta_j = neighPatch.Cf() - neighPatch.Cn();
            const scalarField rho_neigh = rho.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoLimit_neigh = rhoLimit.boundaryField()[patchI].patchNeighbourField();
            const vectorField ULimit_neigh   = ULimit.boundaryField()[patchI].patchNeighbourField();
            const scalarField pLimit_neigh   = pLimit.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar rho_i = rho[i] + (rhoGrad[i]&delta_i[faceI])*rhoLimit[i];
                const scalar rho_j = rho_neigh[faceI] + (rhoGrad_neigh[faceI]&delta_j[faceI])*rhoLimit_neigh[faceI];
                const vector U_i = U[i] + cmptMultiply(UGrad[i]&delta_i[faceI], ULimit[i]);
                const vector U_j = U_neigh[faceI] + cmptMultiply(UGrad_neigh[faceI]&delta_j[faceI], ULimit_neigh[faceI]);
                const scalar p_i = p[i] + (pGrad[i]&delta_i[faceI])*pLimit[i];
                const scalar p_j = p_neigh[faceI] + (pGrad_neigh[faceI]&delta_j[faceI])*pLimit_neigh[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, p_i, p_j,
                                       normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
        if (isA<processorFvPatch>(mesh_.boundary()[patchI]))
        {
            const scalarField rho_neigh = rho.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh   = U.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh   = p.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad.boundaryField()[patchI].patchNeighbourField();
            const tensorField UGrad_neigh   = UGrad.boundaryField()[patchI].patchNeighbourField();
            const vectorField pGrad_neigh   = pGrad.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoLimit_neigh = rhoLimit.boundaryField()[patchI].patchNeighbourField();
            const vectorField ULimit_neigh   = ULimit.boundaryField()[patchI].patchNeighbourField();
            const scalarField pLimit_neigh   = pLimit.boundaryField()[patchI].patchNeighbourField();
            const vectorField delta_i = mesh_.boundary()[patchI].Cf() - mesh_.boundary()[patchI].Cn();
            const vectorField delta_j = mesh_.boundary()[patchI].Cf() - mesh_.C().boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                const scalar rho_i = rho[i] + (rhoGrad[i]&delta_i[faceI])*rhoLimit[i];
                const scalar rho_j = rho_neigh[faceI] + (rhoGrad_neigh[faceI]&delta_j[faceI])*rhoLimit_neigh[faceI];
                const vector U_i = U[i] + cmptMultiply((UGrad[i]&delta_i[faceI]), ULimit[i]);
                const vector U_j = U_neigh[faceI] + cmptMultiply(UGrad_neigh[faceI]&delta_j[faceI], ULimit_neigh[faceI]);
                const scalar p_i = p[i] + (pGrad[i]&delta_i[faceI])*pLimit[i];
                const scalar p_j = p_neigh[faceI] + (pGrad_neigh[faceI]&delta_j[faceI])*pLimit_neigh[faceI];
                riemann_->evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                       rho_i, rho_j, U_i, U_j, p_i, p_j,
                                       normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
        if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
        {
            const fvPatchScalarField& pBound = p.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                vector rhoUFlux = pBound[faceI]*normal[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
            }
        }
        if (name == "inlet" || name == "outlet" || name == "farField")
        {
            const fvPatchScalarField& rhoBound = rho.boundaryField()[patchI];
            const fvPatchVectorField& UBound = U.boundaryField()[patchI];
            const fvPatchScalarField& pBound = p.boundaryField()[patchI];
            forAll(bfaceCells, faceI)
            {
                scalar rhoFlux;
                vector rhoUFlux = vector::zero;
                scalar rhoEFlux;
                const label i = bfaceCells[faceI];
                Foam::evaluateFlux(rhoFlux, rhoUFlux, rhoEFlux,
                                   rhoBound[faceI], UBound[faceI], pBound[faceI],
                                   normal[faceI], fluidProps_.gamma);
                resRho[i]  -= rhoFlux * magSf[faceI];
                resRhoU[i] -= rhoUFlux * magSf[faceI];
                resRhoE[i] -= rhoEFlux * magSf[faceI];
            }
        }
    }
}
