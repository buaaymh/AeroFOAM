void Foam::reconstructionVenkatakrishnan
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    volVectorField& rhoGrad,
    volTensorField& UGrad,
    volVectorField& TGrad
)
{
    const fvMesh& mesh = rho.mesh();
    rhoGrad = fvc::grad(rho);
    UGrad   = fvc::grad(U);
    TGrad   = fvc::grad(T);

    scalarField rhoMin = rho.field();
    scalarField rhoMax = rho.field();
    vectorField UMin = U.field();
    vectorField UMax = U.field();
    scalarField TMin = T.field();
    scalarField TMax = T.field();
    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();
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
    forAll(mesh.boundary(), patchI)
    {
        if(mesh.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh.boundary()[patchI].faceCells();
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
    rhoMin -= rho.field();
    UMin -= U.field();
    TMin -= T.field();
    rhoMax -= rho.field();
    UMax -= U.field();
    TMax -= T.field();

    scalarField rhoLimit(mesh.nCells(), 1.0);
    vectorField ULimit(mesh.nCells(), vector::one);
    scalarField TLimit(mesh.nCells(), 1.0);
    const scalarField delta = mag(mesh.delta());
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const vector delta_i = mesh.Cf()[faceI] - mesh.C()[i];
        const vector delta_j = mesh.Cf()[faceI] - mesh.C()[j];
        const scalar eps2 = 125*Foam::pow3(delta[faceI]);
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
        TLimit[i] = min(TLimit[i], Venkat(TGrad[i]&delta_i, TMin[i], TMax[i], eps2));
        TLimit[j] = min(TLimit[j], Venkat(TGrad[j]&delta_j, TMin[j], TMax[j], eps2));
    }
    forAll(mesh.boundary(), patchI)
    {
        if(mesh.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh.boundary()[patchI].faceCells();
            const vectorField& Cf = mesh.boundary()[patchI].Cf();
            const vectorField Cn = mesh.boundary()[patchI].Cn();
            const scalarField delta = mag(mesh.boundary()[patchI].delta());
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const vector delta_i = Cf[faceI] - Cn[faceI];
                const scalar eps2 = 125*Foam::pow3(delta[faceI]);
                rhoLimit[i] = min(rhoLimit[i], Venkat(rhoGrad[i]&delta_i, rhoMin[i], rhoMax[i], eps2));
                const vector UProj = UGrad[i]&delta_i;
                ULimit[i][0] = min(ULimit[i][0], Venkat(UProj[0], UMin[i][0], UMax[i][0], eps2));
                ULimit[i][1] = min(ULimit[i][1], Venkat(UProj[1], UMin[i][1], UMax[i][1], eps2));
                ULimit[i][2] = min(ULimit[i][2], Venkat(UProj[2], UMin[i][2], UMax[i][2], eps2));
                TLimit[i] = min(TLimit[i], Venkat(TGrad[i]&delta_i, TMin[i], TMax[i], eps2));
            }
        }
    }
    rhoGrad.primitiveFieldRef() *= rhoLimit;
    forAll(mesh.C(), cellI)
    {
        UGrad[cellI].x() *= ULimit[cellI].x();
        UGrad[cellI].y() *= ULimit[cellI].y();
        UGrad[cellI].z() *= ULimit[cellI].z();
    }
    TGrad.primitiveFieldRef() *= TLimit;
    rhoGrad.correctBoundaryConditions();
    UGrad.correctBoundaryConditions();
    TGrad.correctBoundaryConditions();
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
