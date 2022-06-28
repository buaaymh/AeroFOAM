void Foam::solver::evaluateMatrixLDU
(
    scalarField& L,
    scalarField& U,
    scalarField& D
)
{
    const scalar omega = mesh_.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("relaxOmega", 1.2);
    // Calculate L, U and D
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        L[faceI] = 0.5*omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        U[faceI] = 0.5*omega*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI]);
        D[i] += L[faceI];
        D[j] += U[faceI];
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& Sf = mesh_.boundary()[patchI].Sf();
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        forAll(bfaceCells, faceI)
        {
            const label i = bfaceCells[faceI];
            D[i] += 0.5*omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        }
    }
}

void Foam::solver::matrixVectorProduct
(
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
    const scalarField& v_Rho,
    const vectorField& v_RhoU,
    const scalarField& v_RhoE,
    scalarField& y_Rho,
    vectorField& y_RhoU,
    scalarField& y_RhoE
)
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    y_Rho  = D*v_Rho;
    y_RhoU = D*v_RhoU;
    y_RhoE = D*v_RhoE;

    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];

        scalar rho1  = rho_[i]+v_Rho[i];
        vector rhoU1 = rhoU_[i]+v_RhoU[i];
        scalar rhoE1 = rhoE_[i]+v_RhoE[i];
        scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
        scalar phi0 = U_[i]&Sf[faceI];
        scalar phi1 = (rhoU1/rho1)&Sf[faceI];
        y_Rho[j]  -= L[faceI]*v_Rho[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0);
        y_RhoU[j] -= L[faceI]*v_RhoU[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI]);
        y_RhoE[j] -= L[faceI]*v_RhoE[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0);

        rho1  = rho_[j]+v_Rho[j];
        rhoU1 = rhoU_[j]+v_RhoU[j];
        rhoE1 = rhoE_[j]+v_RhoE[j];
        p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
        phi0 = U_[j]&Sf[faceI];
        phi1 = (rhoU1/rho1)&Sf[faceI];
        y_Rho[i]  -= U[faceI]*v_Rho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0);
        y_RhoU[i] -= U[faceI]*v_RhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI]);
        y_RhoE[i] -= U[faceI]*v_RhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0);
    }
}

void Foam::solver::precondition
(
    const scalarField& L,
    const scalarField& U,
    const scalarField& D,
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    // Forward Sweep
    scalarField dRhoStar(resRho*localDtDv_);
    vectorField dRhoUStar(resRhoU*localDtDv_);
    scalarField dRhoEStar(resRhoE*localDtDv_);
    forAll(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rRho = resRho[cellI];
        vector rRhoU = resRhoU[cellI];
        scalar rRhoE = resRhoE[cellI];
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label i = owner[faceI];
                if (i != cellI) /* L(cellI) < cellI */
                {
                    const scalar rho1  = rho_[i]+dRhoStar[i];
                    const vector rhoU1 = rhoU_[i]+dRhoUStar[i];
                    const scalar rhoE1 = rhoE_[i]+dRhoEStar[i];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[i]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += L[faceI]*dRhoStar[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0);
                    rRhoU += L[faceI]*dRhoUStar[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI]);
                    rRhoE += L[faceI]*dRhoEStar[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0);
                }
            }
        }
        dRhoStar[cellI]  = rRho  / D[cellI];
        dRhoUStar[cellI] = rRhoU / D[cellI];
        dRhoEStar[cellI] = rRhoE / D[cellI];
    }
    // Backward Sweep
    scalarField& dRho  = resRho;
    vectorField& dRhoU = resRhoU;
    scalarField& dRhoE = resRhoE;
    forAllReverse(mesh_.cells(), cellI)
    {
        const labelList& cellFaces = mesh_.cells()[cellI];
        scalar rRho  = 0;
        vector rRhoU = vector::zero;
        scalar rRhoE = 0;
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                const label j = neighbour[faceI];
                if (j != cellI) /* U(cellI) > cellI */
                {
                    const scalar rho1  = rho_[j]+dRho[j];
                    const vector rhoU1 = rhoU_[j]+dRhoU[j];
                    const scalar rhoE1 = rhoE_[j]+dRhoE[j];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[j]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += U[faceI]*dRho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0);
                    rRhoU += U[faceI]*dRhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI]);
                    rRhoE += U[faceI]*dRhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0);
                }
            }
        }
        dRho[cellI]  = dRhoStar[cellI]  + rRho  / D[cellI];
        dRhoU[cellI] = dRhoUStar[cellI] + rRhoU / D[cellI];
        dRhoE[cellI] = dRhoEStar[cellI] + rRhoE / D[cellI];
    }
}

void Foam::ApplyGivens
(
    const scalar s,
    const scalar c,
    scalar& h1,
    scalar& h2
)
{
    scalar temp = c*h1 + s*h2;
    h2 = c*h2 - s*h1;
    h1 = temp;
}

scalar Foam::Sign(scalar x, scalar y)
{
    if (y == 0.0) return 0.0;
    return mag(x) * (y < 0.0 ? -1.0 : 1.0);
}

void Foam::GenerateGivens
(
    scalar& dx,
    scalar& dy,
    scalar& s,
    scalar& c
)
{
    if ( (dx == 0.0) && (dy == 0.0) ) {
        c = 1.0;
        s = 0.0;
    }
    else if ( mag(dy) > mag(dx) ) {
        scalar tmp = dx/dy;
        dx = sqrt(1.0 + tmp*tmp);
        s = Sign(1.0/dx, dy);
        c = tmp*s;
    }
    else if ( mag(dy) <= mag(dx) ) {
        scalar tmp = dy/dx;
        dy = sqrt(1.0 + tmp*tmp);
        c = Sign(1.0/dy, dx);
        s = tmp*c;
    }
    else {
        // dx and/or dy must be invalid
        dx = 0.0;
        dy = 0.0;
        c = 1.0;
        s = 0.0;
    }
    dx = mag(dx*dy);
    dy = 0.0;
}

