void Foam::euler2ndSolver::limitGrad()
{
    limiterInit();
    rhoLimit_.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    ULimit_.primitiveFieldRef() = vectorField(mesh_.nCells(), vector::one);
    pLimit_.primitiveFieldRef() = scalarField(mesh_.nCells(), 1.0);
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const vector delta_i = mesh_.Cf()[faceI] - mesh_.C()[i];
        const vector delta_j = mesh_.Cf()[faceI] - mesh_.C()[j];
        const scalar eps2_i = 16*pow(magSqr(delta_i), 1.5);
        const scalar eps2_j = 16*pow(magSqr(delta_j), 1.5);

        rhoLimit_[i] = min(rhoLimit_[i], Venkat(rhoGrad_[i]&delta_i, rhoMin_[i], rhoMax_[i], eps2_i));
        rhoLimit_[j] = min(rhoLimit_[j], Venkat(rhoGrad_[j]&delta_j, rhoMin_[j], rhoMax_[j], eps2_j));
        vector UProj = UGrad_[i]&delta_i;
        vector UMax = UMax_[i];
        vector UMin = UMin_[i];
        ULimit_[i][0] = min(ULimit_[i][0], Venkat(UProj[0], UMin[0], UMax[0], eps2_i));
        ULimit_[i][1] = min(ULimit_[i][1], Venkat(UProj[1], UMin[1], UMax[1], eps2_i));
        ULimit_[i][2] = min(ULimit_[i][2], Venkat(UProj[2], UMin[2], UMax[2], eps2_i));
        UProj = UGrad_[j]&delta_j;
        UMax = UMax_[j];
        UMin = UMin_[j];
        ULimit_[j][0] = min(ULimit_[j][0], Venkat(UProj[0], UMin[0], UMax[0], eps2_j));
        ULimit_[j][1] = min(ULimit_[j][1], Venkat(UProj[1], UMin[1], UMax[1], eps2_j));
        ULimit_[j][2] = min(ULimit_[j][2], Venkat(UProj[2], UMin[2], UMax[2], eps2_j));
        pLimit_[i] = min(pLimit_[i], Venkat(pGrad_[i]&delta_i, pMin_[i], pMax_[i], eps2_i));
        pLimit_[j] = min(pLimit_[j], Venkat(pGrad_[j]&delta_j, pMin_[j], pMax_[j], eps2_j));
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
                const scalar eps2 = 16*pow(magSqr(delta), 1.5);
                rhoLimit_[i] = min(rhoLimit_[i], Venkat(rhoGrad_[i]&delta, rhoMin_[i], rhoMax_[i], eps2));
                const vector UProj = UGrad_[i]&delta;
                const vector& UMax = UMax_[i];
                const vector& UMin = UMin_[i];
                ULimit_[i][0] = min(ULimit_[i][0], Venkat(UProj[0], UMin[0], UMax[0], eps2));
                ULimit_[i][1] = min(ULimit_[i][1], Venkat(UProj[1], UMin[1], UMax[1], eps2));
                ULimit_[i][2] = min(ULimit_[i][2], Venkat(UProj[2], UMin[2], UMax[2], eps2));
                pLimit_[i] = min(pLimit_[i], Venkat(pGrad_[i]&delta, pMin_[i], pMax_[i], eps2));
            }
        }
    }
}

void Foam::euler2ndSolver::limiterInit()
{
    label i, j;
    scalar rho_i, rho_j, p_i, p_j;
    vector U_i, U_j;

    rhoMin_ = rho_.field();
    rhoMax_ = rho_.field();
    UMin_ = U_.field();
    UMax_ = U_.field();
    pMin_ = p_.field();
    pMax_ = p_.field();
    // loop over interior faces
    forAll(mesh_.owner(), faceI)
    {
        i = mesh_.owner()[faceI];
        j = mesh_.neighbour()[faceI];
        // left state
        rho_i = rho_[i];
        U_i = U_[i];
        p_i = p_[i];
        // right state
        rho_j = rho_[j];
        U_j = U_[j];
        p_j = p_[j];
        // neighbors of node i
        rhoMin_[i] = min(rhoMin_[i], rho_j);
        UMin_[i] = min(UMin_[i], U_j);
        pMin_[i] = min(pMin_[i], p_j);
        
        rhoMax_[i] = max(rhoMax_[i], rho_j);
        UMax_[i] = max(UMax_[i], U_j);
        pMax_[i] = max(pMax_[i], p_j);
        // neighbors of node j
        rhoMin_[j] = min(rhoMin_[j], rho_i);
        UMin_[j] = min(UMin_[j], U_i);
        pMin_[j] = min(pMin_[j], p_i);
        
        rhoMax_[j] = max(rhoMax_[j], rho_i);
        UMax_[j] = max(UMax_[j], U_i);
        pMax_[j] = max(pMax_[j], p_i);
    }
    forAll(mesh_.boundary(), patchI)
    {
        if(mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const scalarField rho_neigh  = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField U_neigh = U_.boundaryField()[patchI].patchNeighbourField();
            const scalarField p_neigh = p_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, faceI)
            {
                i = bfaceCells[faceI];
                // right state
                rho_j = rho_neigh[faceI];
                U_j = U_neigh[faceI];
                p_j = p_neigh[faceI];
                // neighbors of node i
                rhoMin_[i] = min(rhoMin_[i], rho_j);
                UMin_[i] = min(UMin_[i], U_j);
                pMin_[i] = min(pMin_[i], p_j);
                
                rhoMax_[i] = max(rhoMax_[i], rho_j);
                UMax_[i] = max(UMax_[i], U_j);
                pMax_[i] = max(pMax_[i], p_j);
            }
        }
    }
    rhoMin_ -= rho_.field();
    UMin_ -= U_.field();
    pMin_ -= p_.field();
    rhoMax_ -= rho_.field();
    UMax_ -= U_.field();
    pMax_ -= p_.field();
}

scalar Foam::euler2ndSolver::Venkat
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