void Foam::solver::solveFlowPseudoTimeSystem
(
    const scalar& dt,
    const scalar& alpha,
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE,
    scalar& resL1
)
{
    updateLTS();
    const scalar omega = mesh_.solutionDict().subDict("LUSGS").lookupOrDefault<scalar>("relaxOmega", 1.2);
    // Calculate D
    scalarField D(1.0/localDtDv_ + mesh_.V().field()/dt);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        D[i] += 0.5*omega*alpha*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        D[j] += 0.5*omega*alpha*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI]);
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const vectorField& Sf = mesh_.boundary()[patchI].Sf();
        const scalarField& magSf = mesh_.boundary()[patchI].magSf();
        forAll(bfaceCells, faceI)
        {
            const label i = bfaceCells[faceI];
            D[i] += 0.5*omega*alpha*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
        }
    }
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
                    const scalar ac = 0.5*omega*(mag(U_[i]&Sf[faceI])+c_[i]*magSf[faceI]);
                    const scalar rho1  = rho_[i]+dRhoStar[i];
                    const vector rhoU1 = rhoU_[i]+dRhoUStar[i];
                    const scalar rhoE1 = rhoE_[i]+dRhoEStar[i];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[i]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += alpha*(ac*dRhoStar[i]  + 0.5*(rho1*phi1       - rho_[i]*phi0));
                    rRhoU += alpha*(ac*dRhoUStar[i] + 0.5*(rhoU1*phi1      - rhoU_[i]*phi0 + (p1-p_[i])*Sf[faceI]));
                    rRhoE += alpha*(ac*dRhoEStar[i] + 0.5*((rhoE1+p1)*phi1 - (rhoE_[i]+p_[i])*phi0));
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
                    const scalar ac = 0.5*omega*(mag(U_[j]&Sf[faceI])+c_[j]*magSf[faceI]);
                    const scalar rho1  = rho_[j]+dRho[j];
                    const vector rhoU1 = rhoU_[j]+dRhoU[j];
                    const scalar rhoE1 = rhoE_[j]+dRhoE[j];
                    const scalar p1    = (rhoE1-0.5*magSqr(rhoU1)/rho1)*(fluidProps_.gamma-1.0);
                    const scalar phi0 = U_[j]&Sf[faceI];
                    const scalar phi1 = (rhoU1/rho1)&Sf[faceI];
                    rRho  += ac*dRho[j]  - 0.5*(rho1*phi1       - rho_[j]*phi0);
                    rRhoU += ac*dRhoU[j] - 0.5*(rhoU1*phi1      - rhoU_[j]*phi0 + (p1-p_[j])*Sf[faceI]);
                    rRhoE += ac*dRhoE[j] - 0.5*((rhoE1+p1)*phi1 - (rhoE_[j]+p_[j])*phi0);
                }
            }
        }
        dRho[cellI]  = dRhoStar[cellI]  + alpha*rRho  / D[cellI];
        dRhoU[cellI] = dRhoUStar[cellI] + alpha*rRhoU / D[cellI];
        dRhoE[cellI] = dRhoEStar[cellI] + alpha*rRhoE / D[cellI];
    }
    resL1 = gSum(mag(dRho));
    rho_.primitiveFieldRef()  += dRho;
    rhoU_.primitiveFieldRef() += dRhoU;
    rhoE_.primitiveFieldRef() += dRhoE;
}