void Foam::euler2ndSolver::limitGrad
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    const volVectorField& rhoGrad,
    const volTensorField& UGrad,
    const volVectorField& TGrad,
    volScalarField& rhoLimit,
    volVectorField& ULimit,
    volScalarField& TLimit
) const
{
    scalarField rhoMin = rho.field();
    scalarField rhoMax = rho.field();
    vectorField UMin = U.field();
    vectorField UMax = U.field();
    scalarField TMin = T.field();
    scalarField TMax = T.field();
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        // neighbors of node i
        rhoMin[i] = min(rhoMin[i], rho[j]);
        UMin[i] = min(UMin[i], U[j]);
        TMin[i] = min(TMin[i], T[j]);
        
        rhoMax[i] = max(rhoMax[i], rho[j]);
        UMax[i] = max(UMax[i], U[j]);
        TMax[i] = max(TMax[i], T[j]);
        // neighbors of node j
        rhoMin[j] = min(rhoMin[j], rho[i]);
        UMin[j] = min(UMin[j], U[i]);
        TMin[j] = min(TMin[j], T[i]);
        
        rhoMax[j] = max(rhoMax[j], rho[i]);
        UMax[j] = max(UMax[j], U[i]);
        TMax[j] = max(TMax[j], T[i]);
    }
    forAll(mesh_.boundary(), patchI)
    {
        if(mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const scalarField rho_neigh  = rho.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh = U.boundaryField()[patchI].patchNeighbourField();
            const scalarField T_neigh = T.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                // neighbors of node i
                rhoMin[i] = min(rhoMin[i], rho_neigh[faceI]);
                UMin[i] = min(UMin[i], U_neigh[faceI]);
                TMin[i] = min(TMin[i], T_neigh[faceI]);
                
                rhoMax[i] = max(rhoMax[i], rho_neigh[faceI]);
                UMax[i] = max(UMax[i], U_neigh[faceI]);
                TMax[i] = max(TMax[i], T_neigh[faceI]);
            }
        }
    }
    rhoMin -= rho_.field();
    UMin -= U.field();
    TMin -= T.field();
    rhoMax -= rho.field();
    UMax -= U.field();
    TMax -= T.field();
    rhoLimit.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    ULimit.primitiveFieldRef() = vectorField(mesh_.nCells(), vector::one);
    TLimit.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        const scalar eps2_i = 125*pow(magSqr(delta_i), 1.5);
        const scalar eps2_j = 125*pow(magSqr(delta_j), 1.5);
        rhoLimit[i] = min(rhoLimit[i], Venkat(rhoGrad[i]&delta_i, rhoMin[i], rhoMax[i], eps2_i));
        rhoLimit[j] = min(rhoLimit[j], Venkat(rhoGrad[j]&delta_j, rhoMin[j], rhoMax[j], eps2_j));
        vector UProj = UGrad[i]&delta_i;
        ULimit[i][0] = min(ULimit[i][0], Venkat(UProj[0], UMin[i][0], UMax[i][0], eps2_i));
        ULimit[i][1] = min(ULimit[i][1], Venkat(UProj[1], UMin[i][1], UMax[i][1], eps2_i));
        ULimit[i][2] = min(ULimit[i][2], Venkat(UProj[2], UMin[i][2], UMax[i][2], eps2_i));
        UProj = UGrad[j]&delta_j;
        ULimit[j][0] = min(ULimit[j][0], Venkat(UProj[0], UMin[j][0], UMax[j][0], eps2_j));
        ULimit[j][1] = min(ULimit[j][1], Venkat(UProj[1], UMin[j][1], UMax[j][1], eps2_j));
        ULimit[j][2] = min(ULimit[j][2], Venkat(UProj[2], UMin[j][2], UMax[j][2], eps2_j));
        TLimit[i] = min(TLimit[i], Venkat(TGrad[i]&delta_i, TMin[i], TMax[i], eps2_i));
        TLimit[j] = min(TLimit[j], Venkat(TGrad[j]&delta_j, TMin[j], TMax[j], eps2_j));
    }
    forAll(mesh_.boundary(), patchI)
    {
        if(mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const vectorField& Cf = mesh_.boundary()[patchI].Cf();
            const vectorField Cn = mesh_.boundary()[patchI].Cn();
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const vector delta = Cf[faceI] - Cn[faceI];
                const scalar eps2 = 125*pow(magSqr(delta), 1.5);
                rhoLimit[i] = min(rhoLimit[i], Venkat(rhoGrad[i]&delta, rhoMin[i], rhoMax[i], eps2));
                const vector UProj = UGrad[i]&delta;
                ULimit[i][0] = min(ULimit[i][0], Venkat(UProj[0], UMin[i][0], UMax[i][0], eps2));
                ULimit[i][1] = min(ULimit[i][1], Venkat(UProj[1], UMin[i][1], UMax[i][1], eps2));
                ULimit[i][2] = min(ULimit[i][2], Venkat(UProj[2], UMin[i][2], UMax[i][2], eps2));
                TLimit[i] = min(TLimit[i], Venkat(TGrad[i]&delta, TMin[i], TMax[i], eps2));
            }
        }
    }
    rhoLimit.correctBoundaryConditions();
    ULimit.correctBoundaryConditions();
    TLimit.correctBoundaryConditions();
}

scalar Foam::euler2ndSolver::Venkat
(
    const scalar& project,
    const scalar& deltaMin,
    const scalar& deltaMax,
    const scalar& eps2
) const
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