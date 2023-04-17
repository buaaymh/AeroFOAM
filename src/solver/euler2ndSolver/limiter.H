void Foam::euler2ndSolver::deltaMaxMin
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p,
    scalarField& rhoMin,
    scalarField& rhoMax,
    vectorField& UMin,
    vectorField& UMax,
    scalarField& pMin,
    scalarField& pMax
) const
{
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        // neighbors of node i
        rhoMin[i] = min(rhoMin[i], rho[j]);
        rhoMax[i] = max(rhoMax[i], rho[j]);
        UMin[i] = min(UMin[i], U[j]);
        UMax[i] = max(UMax[i], U[j]);
        pMin[i] = min(pMin[i], p[j]);
        pMax[i] = max(pMax[i], p[j]);

        // neighbors of node j
        rhoMin[j] = min(rhoMin[j], rho[i]);
        rhoMax[j] = max(rhoMax[j], rho[i]);
        UMin[j] = min(UMin[j], U[i]);
        UMax[j] = max(UMax[j], U[i]);
        pMin[j] = min(pMin[j], p[i]);
        pMax[j] = max(pMax[j], p[i]);
    }
    forAll(mesh_.boundary(), patchI)
    {
        if(mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const scalarField rho_neigh  = rho.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh = U.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh = p.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                // neighbors of node i
                rhoMin[i] = min(rhoMin[i], rho_neigh[faceI]);
                rhoMax[i] = max(rhoMax[i], rho_neigh[faceI]);
                UMin[i] = min(UMin[i], U_neigh[faceI]);
                UMax[i] = max(UMax[i], U_neigh[faceI]);
                pMin[i] = min(pMin[i], p_neigh[faceI]);
                pMax[i] = max(pMax[i], p_neigh[faceI]);
            }
        }
    }
    rhoMin -= rho.field();
    UMin -= U.field();
    pMin -= p.field();
    rhoMax -= rho.field();
    UMax -= U.field();
    pMax -= p.field();
}

void Foam::euler2ndSolver::limiterVenkatakrishnan()
{
    
    const scalar K = mesh_.schemesDict().subDict("gradSchemes").lookupOrDefault<scalar>("Venkat", 1);
    scalarField rhoMin = rho_.field();
    scalarField rhoMax = rho_.field();
    vectorField UMin = U_.field();
    vectorField UMax = U_.field();
    scalarField pMin = p_.field();
    scalarField pMax = p_.field();
    deltaMaxMin(rho_, U_, p_, rhoMin, rhoMax, UMin, UMax, pMin, pMax);
    rhoLimit_.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    ULimit_.primitiveFieldRef() = vectorField(mesh_.nCells(), vector::one);
    pLimit_.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    const scalarField delta = mag(mesh_.delta());
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        const scalar eps2 = Foam::pow3(K*delta[faceI]);
        rhoLimit_[i] = min(rhoLimit_[i], Venkat(rhoGrad_[i]&delta_i, rhoMin[i], rhoMax[i], eps2));
        rhoLimit_[j] = min(rhoLimit_[j], Venkat(rhoGrad_[j]&delta_j, rhoMin[j], rhoMax[j], eps2));
        vector UProj = UGrad_[i]&delta_i;
        ULimit_[i][0] = min(ULimit_[i][0], Venkat(UProj[0], UMin[i][0], UMax[i][0], eps2));
        ULimit_[i][1] = min(ULimit_[i][1], Venkat(UProj[1], UMin[i][1], UMax[i][1], eps2));
        ULimit_[i][2] = min(ULimit_[i][2], Venkat(UProj[2], UMin[i][2], UMax[i][2], eps2));
        UProj = UGrad_[j]&delta_j;
        ULimit_[j][0] = min(ULimit_[j][0], Venkat(UProj[0], UMin[j][0], UMax[j][0], eps2));
        ULimit_[j][1] = min(ULimit_[j][1], Venkat(UProj[1], UMin[j][1], UMax[j][1], eps2));
        ULimit_[j][2] = min(ULimit_[j][2], Venkat(UProj[2], UMin[j][2], UMax[j][2], eps2));
        pLimit_[i] = min(pLimit_[i], Venkat(pGrad_[i]&delta_i, pMin[i], pMax[i], eps2));
        pLimit_[j] = min(pLimit_[j], Venkat(pGrad_[j]&delta_j, pMin[j], pMax[j], eps2));
    }
    forAll(mesh_.boundary(), patchI)
    {
        if(mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const vectorField& Cf = mesh_.boundary()[patchI].Cf();
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            const scalarField delta = mag(mesh_.boundary()[patchI].delta());
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const vector delta_i = Cf[faceI] - Cn[faceI];
                const scalar eps2 = Foam::pow3(K*delta[faceI]);
                rhoLimit_[i] = min(rhoLimit_[i], Venkat(rhoGrad_[i]&delta_i, rhoMin[i], rhoMax[i], eps2));
                const vector UProj = UGrad_[i]&delta_i;
                ULimit_[i][0] = min(ULimit_[i][0], Venkat(UProj[0], UMin[i][0], UMax[i][0], eps2));
                ULimit_[i][1] = min(ULimit_[i][1], Venkat(UProj[1], UMin[i][1], UMax[i][1], eps2));
                ULimit_[i][2] = min(ULimit_[i][2], Venkat(UProj[2], UMin[i][2], UMax[i][2], eps2));
                pLimit_[i] = min(pLimit_[i], Venkat(pGrad_[i]&delta_i, pMin[i], pMax[i], eps2));
            }
        }
    }
    rhoLimit_.correctBoundaryConditions();
    ULimit_.correctBoundaryConditions();
    pLimit_.correctBoundaryConditions();
}

void Foam::euler2ndSolver::limiterBarthJespersen()
{
    scalarField rhoMin = rho_.field();
    scalarField rhoMax = rho_.field();
    vectorField UMin = U_.field();
    vectorField UMax = U_.field();
    scalarField pMin = p_.field();
    scalarField pMax = p_.field();
    deltaMaxMin(rho_, U_, p_, rhoMin, rhoMax, UMin, UMax, pMin, pMax);
    rhoLimit_.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    ULimit_.primitiveFieldRef() = vectorField(mesh_.nCells(), vector::one);
    pLimit_.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    const scalarField delta = mag(mesh_.delta());
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        rhoLimit_[i] = min(rhoLimit_[i], Barth(rhoGrad_[i]&delta_i, rhoMin[i], rhoMax[i]));
        rhoLimit_[j] = min(rhoLimit_[j], Barth(rhoGrad_[j]&delta_j, rhoMin[j], rhoMax[j]));
        vector UProj = UGrad_[i]&delta_i;
        ULimit_[i][0] = min(ULimit_[i][0], Barth(UProj[0], UMin[i][0], UMax[i][0]));
        ULimit_[i][1] = min(ULimit_[i][1], Barth(UProj[1], UMin[i][1], UMax[i][1]));
        ULimit_[i][2] = min(ULimit_[i][2], Barth(UProj[2], UMin[i][2], UMax[i][2]));
        UProj = UGrad_[j]&delta_j;
        ULimit_[j][0] = min(ULimit_[j][0], Barth(UProj[0], UMin[j][0], UMax[j][0]));
        ULimit_[j][1] = min(ULimit_[j][1], Barth(UProj[1], UMin[j][1], UMax[j][1]));
        ULimit_[j][2] = min(ULimit_[j][2], Barth(UProj[2], UMin[j][2], UMax[j][2]));
        pLimit_[i] = min(pLimit_[i], Barth(pGrad_[i]&delta_i, pMin[i], pMax[i]));
        pLimit_[j] = min(pLimit_[j], Barth(pGrad_[j]&delta_j, pMin[j], pMax[j]));
    }
    forAll(mesh_.boundary(), patchI)
    {
        if(mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const vectorField& Cf = mesh_.boundary()[patchI].Cf();
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            const scalarField delta = mag(mesh_.boundary()[patchI].delta());
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const vector delta_i = Cf[faceI] - Cn[faceI];
                rhoLimit_[i] = min(rhoLimit_[i], Barth(rhoGrad_[i]&delta_i, rhoMin[i], rhoMax[i]));
                const vector UProj = UGrad_[i]&delta_i;
                ULimit_[i][0] = min(ULimit_[i][0], Barth(UProj[0], UMin[i][0], UMax[i][0]));
                ULimit_[i][1] = min(ULimit_[i][1], Barth(UProj[1], UMin[i][1], UMax[i][1]));
                ULimit_[i][2] = min(ULimit_[i][2], Barth(UProj[2], UMin[i][2], UMax[i][2]));
                pLimit_[i] = min(pLimit_[i], Barth(pGrad_[i]&delta_i, pMin[i], pMax[i]));
            }
        }
    }
    rhoLimit_.correctBoundaryConditions();
    ULimit_.correctBoundaryConditions();
    pLimit_.correctBoundaryConditions();
}
