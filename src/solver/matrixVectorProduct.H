void Foam::solver::matrixVectorProduct
(
    const scalarField& resRho,
    const vectorField& resRhoU,
    const scalarField& resRhoE,
    const scalarField& v_Rho,
    const vectorField& v_RhoU,
    const scalarField& v_RhoE,
    scalarField& y_Rho,
    vectorField& y_RhoU,
    scalarField& y_RhoE
)
{
    const scalar eps = 1e-8/Foam::sqrt(gSumSqr(v_Rho)+gSum(magSqr(v_RhoU))+gSumSqr(v_RhoE)+1e-16);
    volScalarField rho  = rho_;
    volVectorField rhoU = rhoU_;
    volScalarField rhoE = rhoE_;
    rho_.primitiveFieldRef()  += v_Rho*eps;
    rhoU_.primitiveFieldRef() += v_RhoU*eps;
    rhoE_.primitiveFieldRef() += v_RhoE*eps;
    correctFields();
    solver::evaluateFlowRes(y_Rho, y_RhoU, y_RhoE);
    y_Rho  = v_Rho /localDtDv_ - (y_Rho  - resRho) /eps;
    y_RhoU = v_RhoU/localDtDv_ - (y_RhoU - resRhoU)/eps;
    y_RhoE = v_RhoE/localDtDv_ - (y_RhoE - resRhoE)/eps;
    rho_  = rho;
    rhoU_ = rhoU;
    rhoE_ = rhoE;
}

void Foam::solver::evaluateFlowRes
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
    if (fluidProps_.simulationType != "Euler")
    {
        UGrad_ = fvc::grad(U_);
        TGrad_ = fvc::grad(T_);
        solver::evaluateViscousFlowRes(resRho, resRhoU, resRhoE);
    }
}

void Foam::solver::evaluateViscousFlowRes
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{

}