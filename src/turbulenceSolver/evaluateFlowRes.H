void Foam::turbulenceSolver::evaluateFlowResForGMRES
(
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    eulerSolver::evaluateFlowResForGMRES(resRho, resRhoU, resRhoE);
    UGrad_   = fvc::grad(U_);
    TGrad_   = fvc::grad(T_);
    // forAll(mesh_.owner(), faceI)
    // {
    //     vector rhoUFlux = vector::zero;
    //     scalar rhoEFlux = 0;
    //     const label i = mesh_.owner()[faceI];
    //     const label j = mesh_.neighbour()[faceI];
    //     vector TGrad = AverageGrad(T_[i], T_[j], TGrad_[i], TGrad_[j], delta_[faceI]);
    //     tensor UGrad = AverageGrad(U_[i], U_[j], UGrad_[i], UGrad_[j], delta_[faceI]);
    //     evaluateViscousFlux(rhoUFlux, rhoEFlux, 0.5*(U_[i]+U_[j]), TGrad, UGrad,
    //                         0.5*(muLam_[i]+muLam_[j]), normal_[faceI]);
    //     resRhoU[i] += rhoUFlux*mesh_.magSf()[faceI];
    //     resRhoE[i] += rhoEFlux*mesh_.magSf()[faceI];
    //     resRhoU[j] -= rhoUFlux*mesh_.magSf()[faceI];
    //     resRhoE[j] -= rhoEFlux*mesh_.magSf()[faceI];
    // }
    // forAll(mesh_.boundary(), patchI)
    // {
    //     const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
    //     const vectorField& normal = normal_.boundaryField()[patchI];
    //     const scalarField& magSf = mesh_.magSf().boundaryField()[patchI];
    //     const word name = mesh_.boundary()[patchI].name();
    //     const word type = mesh_.boundary()[patchI].type();
    //     if (mesh_.boundary()[patchI].coupled())
    //     {
    //         const vectorField delta = mesh_.boundary()[patchI].delta();
    //         const vectorField U_neigh = U_.boundaryField()[patchI].patchNeighbourField();
    //         const scalarField T_neigh = T_.boundaryField()[patchI].patchNeighbourField();
    //         const scalarField muLam_neigh = muLam_.boundaryField()[patchI].patchNeighbourField();
    //         const tensorField UGrad_neigh = UGrad_.boundaryField()[patchI].patchNeighbourField();
    //         const vectorField TGrad_neigh = TGrad_.boundaryField()[patchI].patchNeighbourField();
    //         forAll(bfaceCells, faceI)
    //         {
    //             vector rhoUFlux;
    //             scalar rhoEFlux;
    //             const label i = bfaceCells[faceI];
    //             vector TGrad = AverageGrad(T_[i], T_neigh[faceI], TGrad_[i], TGrad_neigh[faceI], delta[faceI]);
    //             tensor UGrad = AverageGrad(U_[i], U_neigh[faceI], UGrad_[i], UGrad_neigh[faceI], delta[faceI]);
    //             evaluateViscousFlux(rhoUFlux, rhoEFlux, 0.5*(U_[i]+U_neigh[faceI]), TGrad, UGrad,
    //                                 0.5*(muLam_[i]+muLam_neigh[faceI]), normal[faceI]);
    //             resRhoU[i] += rhoUFlux*magSf[faceI];
    //             resRhoE[i] += rhoEFlux*magSf[faceI];
    //         }
    //     }
    //     if (type == "wall")
    //     {
    //         forAll(bfaceCells, faceI)
    //         {
    //             vector rhoUFlux;
    //             scalar rhoEFlux;
    //             const label i = bfaceCells[faceI];
    //             evaluateViscousFlux(rhoUFlux, rhoEFlux, vector::zero, vector::zero, UGrad_[i],
    //                                 muLam_[i], normal[faceI]);
    //             resRhoU[i] += rhoUFlux * magSf[faceI];
    //         }
    //     }
    //     if (name == "inlet" || name == "outlet" || name == "farField")
    //     {
    //         const fvPatchVectorField& UBound = U_.boundaryField()[patchI];
    //         const fvPatchScalarField& muLamBound = muLam_.boundaryField()[patchI];
    //         forAll(bfaceCells, faceI)
    //         {
    //             vector rhoUFlux;
    //             scalar rhoEFlux;
    //             const label i = bfaceCells[faceI];
    //             evaluateViscousFlux(rhoUFlux, rhoEFlux, UBound[faceI], TGrad_[i], UGrad_[i],
    //                                 muLamBound[faceI], normal[faceI]);
    //             resRhoU[i] += rhoUFlux * magSf[faceI];
    //             resRhoE[i] += rhoEFlux * magSf[faceI];
    //         }
    //     }
    // }
}

// void Foam::navierStokesSolver::evaluateViscousFlux
// (
//     vector& rhoUFlux,
//     scalar& rhoEFlux,
//     const vector& U,
//     const vector& TGrad,
//     const tensor& UGrad,
//     const scalar& muLam,
//     const vector& normal
// ) const
// {
//     const scalar mu = muLam;
//     const scalar lambda = -2.0/3.0*mu;
//     const scalar kappa  = (muLam/Pr_Lam)/(Gamma-1.0);
//     const vector StrainRateU = (UGrad + UGrad.T())&normal;
//     const scalar divU = UGrad.xx() + UGrad.yy() + UGrad.zz();
//     const vector rhoUViscousFlux = mu*StrainRateU + lambda*divU*normal;
//     rhoUFlux = Ma_Re_*rhoUViscousFlux;
//     rhoEFlux = Ma_Re_*((rhoUViscousFlux&U)+kappa*(TGrad&normal));
// }
