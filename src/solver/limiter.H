void Foam::solver::deltaMaxMin
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

void Foam::solver::limiterVenkatakrishnan
(
    const scalar K,
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p,
    volVectorField& rhoGrad,
    volTensorField& UGrad,
    volVectorField& pGrad
) const
{
    scalarField rhoMin = rho.field();
    scalarField rhoMax = rho.field();
    vectorField UMin = U.field();
    vectorField UMax = U.field();
    scalarField pMin = p.field();
    scalarField pMax = p.field();
    deltaMaxMin(rho, U, p, rhoMin, rhoMax, UMin, UMax, pMin, pMax);
    scalarField rhoLimit(mesh_.nCells(), 1.0);
    vectorField ULimit(mesh_.nCells(), vector::one);
    scalarField pLimit(mesh_.nCells(), 1.0);
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
        rhoLimit[i] = min(rhoLimit[i], Venkat(rhoGrad[i]&delta_i, rhoMin[i], rhoMax[i], eps2));
        rhoLimit[j] = min(rhoLimit[j], Venkat(rhoGrad[j]&delta_j, rhoMin[j], rhoMax[j], eps2));
        vector UProj = UGrad[i]&delta_i;
        ULimit[i][0] = min(ULimit[i][0], Venkat(UProj[0], UMin[i][0], UMax[i][0], eps2));
        ULimit[i][1] = min(ULimit[i][1], Venkat(UProj[1], UMin[i][1], UMax[i][1], eps2));
        ULimit[i][2] = min(ULimit[i][2], Venkat(UProj[2], UMin[i][2], UMax[i][2], eps2));
        UProj = UGrad[j]&delta_j;
        ULimit[j][0] = min(ULimit[j][0], Venkat(UProj[0], UMin[j][0], UMax[j][0], eps2));
        ULimit[j][1] = min(ULimit[j][1], Venkat(UProj[1], UMin[j][1], UMax[j][1], eps2));
        ULimit[j][2] = min(ULimit[j][2], Venkat(UProj[2], UMin[j][2], UMax[j][2], eps2));
        pLimit[i] = min(pLimit[i], Venkat(pGrad[i]&delta_i, pMin[i], pMax[i], eps2));
        pLimit[j] = min(pLimit[j], Venkat(pGrad[j]&delta_j, pMin[j], pMax[j], eps2));
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
                rhoLimit[i] = min(rhoLimit[i], Venkat(rhoGrad[i]&delta_i, rhoMin[i], rhoMax[i], eps2));
                const vector UProj = UGrad[i]&delta_i;
                ULimit[i][0] = min(ULimit[i][0], Venkat(UProj[0], UMin[i][0], UMax[i][0], eps2));
                ULimit[i][1] = min(ULimit[i][1], Venkat(UProj[1], UMin[i][1], UMax[i][1], eps2));
                ULimit[i][2] = min(ULimit[i][2], Venkat(UProj[2], UMin[i][2], UMax[i][2], eps2));
                pLimit[i] = min(pLimit[i], Venkat(pGrad[i]&delta_i, pMin[i], pMax[i], eps2));
            }
        }
    }
    rhoGrad.primitiveFieldRef() *= rhoLimit;
    forAll(mesh_.C(), cellI)
    {
        UGrad[cellI] = tensor(UGrad[cellI].x()*ULimit[cellI].x(),
                              UGrad[cellI].y()*ULimit[cellI].y(),
                              UGrad[cellI].z()*ULimit[cellI].z());
    }
    pGrad.primitiveFieldRef() *= pLimit;
    rhoGrad.correctBoundaryConditions();
    UGrad.correctBoundaryConditions();
    pGrad.correctBoundaryConditions();
}

void Foam::solver::limiterBarthJespersen
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& p,
    volVectorField& rhoGrad,
    volTensorField& UGrad,
    volVectorField& pGrad
) const
{
    scalarField rhoMin = rho.field();
    scalarField rhoMax = rho.field();
    vectorField UMin = U.field();
    vectorField UMax = U.field();
    scalarField pMin = p.field();
    scalarField pMax = p.field();
    deltaMaxMin(rho, U, p, rhoMin, rhoMax, UMin, UMax, pMin, pMax);
    scalarField rhoLimit(mesh_.nCells(), 1.0);
    vectorField ULimit(mesh_.nCells(), vector::one);
    scalarField pLimit(mesh_.nCells(), 1.0);
    const scalarField delta = mag(mesh_.delta());
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        rhoLimit[i] = min(rhoLimit[i], Barth(rhoGrad[i]&delta_i, rhoMin[i], rhoMax[i]));
        rhoLimit[j] = min(rhoLimit[j], Barth(rhoGrad[j]&delta_j, rhoMin[j], rhoMax[j]));
        vector UProj = UGrad[i]&delta_i;
        ULimit[i][0] = min(ULimit[i][0], Barth(UProj[0], UMin[i][0], UMax[i][0]));
        ULimit[i][1] = min(ULimit[i][1], Barth(UProj[1], UMin[i][1], UMax[i][1]));
        ULimit[i][2] = min(ULimit[i][2], Barth(UProj[2], UMin[i][2], UMax[i][2]));
        UProj = UGrad[j]&delta_j;
        ULimit[j][0] = min(ULimit[j][0], Barth(UProj[0], UMin[j][0], UMax[j][0]));
        ULimit[j][1] = min(ULimit[j][1], Barth(UProj[1], UMin[j][1], UMax[j][1]));
        ULimit[j][2] = min(ULimit[j][2], Barth(UProj[2], UMin[j][2], UMax[j][2]));
        pLimit[i] = min(pLimit[i], Barth(pGrad[i]&delta_i, pMin[i], pMax[i]));
        pLimit[j] = min(pLimit[j], Barth(pGrad[j]&delta_j, pMin[j], pMax[j]));
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
                rhoLimit[i] = min(rhoLimit[i], Barth(rhoGrad[i]&delta_i, rhoMin[i], rhoMax[i]));
                const vector UProj = UGrad[i]&delta_i;
                ULimit[i][0] = min(ULimit[i][0], Barth(UProj[0], UMin[i][0], UMax[i][0]));
                ULimit[i][1] = min(ULimit[i][1], Barth(UProj[1], UMin[i][1], UMax[i][1]));
                ULimit[i][2] = min(ULimit[i][2], Barth(UProj[2], UMin[i][2], UMax[i][2]));
                pLimit[i] = min(pLimit[i], Barth(pGrad[i]&delta_i, pMin[i], pMax[i]));
            }
        }
    }
    rhoGrad.primitiveFieldRef() *= rhoLimit;
    forAll(mesh_.C(), cellI)
    {
        UGrad[cellI] = tensor(UGrad[cellI].x()*ULimit[cellI].x(),
                              UGrad[cellI].y()*ULimit[cellI].y(),
                              UGrad[cellI].z()*ULimit[cellI].z());
    }
    pGrad.primitiveFieldRef() *= pLimit;
    rhoGrad.correctBoundaryConditions();
    UGrad.correctBoundaryConditions();
    pGrad.correctBoundaryConditions();
}


scalar Foam::Venkat
(
    const scalar& project,
    const scalar& deltaMin,
    const scalar& deltaMax,
    const scalar& eps2
)
{
    if (project > SMALL)
    {
        scalar y = deltaMax*(deltaMax+project) + eps2;
        return (y + deltaMax*project) / (y + 2*project*project);
    }
    else if (project < -SMALL)
    {
        scalar y = deltaMin*(deltaMin+project) + eps2;
        return (y + deltaMin*project) / (y + 2*project*project);
    }
    return 1.0;
}

scalar Foam::Barth
(
    const scalar& project,
    const scalar& deltaMin,
    const scalar& deltaMax
)
{
    if (project > SMALL)
    {
        return deltaMax / project;
    }
    else if (project < -SMALL)
    {
        return deltaMin / project;
    }
    return 1.0;
}