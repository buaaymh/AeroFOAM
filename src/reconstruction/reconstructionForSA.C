/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "reconstructionForSA.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::ReconstructionForSA::ReconstructionForSA
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& rhoU,
    volScalarField& rhoE,
    volScalarField& nuTilda
)
:
    Reconstruction(fluidProps, rho, rhoU, rhoE),
    nuTilda_(nuTilda),
    coefsNuTilda_(mesh_.nCells(), Col9X1::Zero()),
    d1NuTilda_
    (
        IOobject
        (
            "d1NuTilda",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    d2NuTilda_
    (
        IOobject
        (
            "d2NuTilda",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    )
{
    cellQuad_.reserve(mesh_.nCells());
    // Wall-reflection vectors
    const volVectorField& n = wallDist::New(mesh_).n();
    // Wall distance
    const volScalarField& y = wallDist::New(mesh_).y();
    forAll(mesh_.C(), cellI)
    {
        cellQuad_.emplace_back(mesh_, cellI, y[cellI], n[cellI]);
    }
}

void Foam::ReconstructionForSA::iterationStep()
{
    /* Calculate b of vr linear system */
    std::vector<Mat9X5> bW(mesh_.nCells(), Mat9X5::Zero());
    std::vector<Col9X1> bN(mesh_.nCells(), Col9X1::Zero());
    const auto& owner = mesh_.owner();
    const auto& neighbour = mesh_.neighbour();
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        const scalar dRho0  = rho_[j]  - rho_[i];
        const vector dRhoU0 = rhoU_[j] - rhoU_[i];
        const scalar dRhoE0 = rhoE_[j] - rhoE_[i];
        const Col5X1 dVars(dRho0, dRhoU0[0], dRhoU0[1], dRhoU0[2], dRhoE0);
        const scalar dNuTilda0 = nuTilda_[j] - nuTilda_[i];
        bW[i] += lowerb_[faceI]*dVars.transpose();
        bN[i] += lowerb_[faceI]*dNuTilda0;
        bW[j] -= upperb_[faceI]*dVars.transpose();
        bN[j] -= upperb_[faceI]*dNuTilda0;
    }
    forAll(mesh_.boundary(), patchI)
    {
        const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
        const word type = mesh_.boundary()[patchI].type();
        if (mesh_.boundary()[patchI].coupled())
        {
            const scalarField rho_neigh  = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoU_neigh = rhoU_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rhoE_neigh = rhoE_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1Rho_neigh  = d1Rho_.boundaryField()[patchI].patchNeighbourField();
            const tensorField d1RhoU_neigh = d1RhoU_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1RhoE_neigh = d1RhoE_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUx_neigh = d2RhoUx_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUy_neigh = d2RhoUy_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoUz_neigh = d2RhoUz_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2RhoE_neigh  = d2RhoE_.boundaryField()[patchI].patchNeighbourField();
            const scalarField nuTilda_neigh  = nuTilda_.boundaryField()[patchI].patchNeighbourField();
            const vectorField d1NuTilda_neigh  = d1NuTilda_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2NuTilda_neigh = d2NuTilda_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                const scalar dRho0  = rho_neigh[j]  - rho_[i];
                const vector dRhoU0 = rhoU_neigh[j] - rhoU_[i];
                const scalar dRhoE0 = rhoE_neigh[j] - rhoE_[i];
                const Col5X1 dVars(dRho0, dRhoU0[0], dRhoU0[1], dRhoU0[2], dRhoE0);
                const scalar dNuTilda0 = nuTilda_neigh[j] - nuTilda_[i];
                bW[i] += coupledb_[patchI][j]*dVars.transpose();
                bW[i] += (Mat5X9
                {
                    {d1Rho_neigh[j][0], d1Rho_neigh[j][1], d1Rho_neigh[j][2],
                     d2Rho_neigh[j][0], d2Rho_neigh[j][1], d2Rho_neigh[j][2],
                     d2Rho_neigh[j][3], d2Rho_neigh[j][4], d2Rho_neigh[j][5]},
                    {d1RhoU_neigh[j][0],  d1RhoU_neigh[j][1],  d1RhoU_neigh[j][2],
                     d2RhoUx_neigh[j][0], d2RhoUx_neigh[j][1], d2RhoUx_neigh[j][2],
                     d2RhoUx_neigh[j][3], d2RhoUx_neigh[j][4], d2RhoUx_neigh[j][5]},
                    {d1RhoU_neigh[j][3],  d1RhoU_neigh[j][4],  d1RhoU_neigh[j][5],
                     d2RhoUy_neigh[j][0], d2RhoUy_neigh[j][1], d2RhoUy_neigh[j][2],
                     d2RhoUy_neigh[j][3], d2RhoUy_neigh[j][4], d2RhoUy_neigh[j][5]},
                    {d1RhoU_neigh[j][6],  d1RhoU_neigh[j][7],  d1RhoU_neigh[j][8],
                     d2RhoUz_neigh[j][0], d2RhoUz_neigh[j][1], d2RhoUz_neigh[j][2],
                     d2RhoUz_neigh[j][3], d2RhoUz_neigh[j][4], d2RhoUz_neigh[j][5]},
                    {d1RhoE_neigh[j][0], d1RhoE_neigh[j][1], d1RhoE_neigh[j][2],
                     d2RhoE_neigh[j][0], d2RhoE_neigh[j][1], d2RhoE_neigh[j][2],
                     d2RhoE_neigh[j][3], d2RhoE_neigh[j][4], d2RhoE_neigh[j][5]}
                } * coupledB_[patchI][j]).transpose();
                bN[i] += coupledb_[patchI][j]*dNuTilda0;
                bN[i] += coupledB_[patchI][j].transpose() *
                         Col9X1(d1NuTilda_neigh[j][0], d1NuTilda_neigh[j][1], d1NuTilda_neigh[j][2],
                                d2NuTilda_neigh[j][0], d2NuTilda_neigh[j][1], d2NuTilda_neigh[j][2],
                                d2NuTilda_neigh[j][3], d2NuTilda_neigh[j][4], d2NuTilda_neigh[j][5]);
            }
        }
        if (type == "symmetryPlane" || type == "symmetry" || (type == "wall" && fluidProps_.simulationType == "Euler"))
        {
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const scalarField d_ij = mag(mesh_.boundary()[patchI].delta());
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch4stQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const vector dp = (0.5/d_ij[faceI])*vector(parameter_.vrWeightSqr[0], 0, 0);
                Mat9X5 b = Mat9X5::Zero();
                Col9X1 b_scalar = Col9X1::Zero();
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - ownerCn[faceI];
                    const vector dpwTemp = boundaryQuad[faceI]->weight(gaussI)*dp;
                    Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    Col5X1 varDn0 = coefs_[i].transpose()*Dn0_i;
                    // Dn0
                    const vector rhoUDn0_i = rhoU_[i] + vector(varDn0(1), varDn0(2), varDn0(3));
                    const vector rhoUDn0_j = rhoUDn0_i - 2*(rhoUDn0_i&boundaryQuad[faceI]->normal(gaussI))
                                           * boundaryQuad[faceI]->normal(gaussI);
                    const scalar dRhoDn0  = dpwTemp[0]*varDn0(0);
                    const vector dRhoUDn0 = dpwTemp[0]*(rhoUDn0_j - rhoU_[i]);
                    const scalar dRhoEDn0 = dpwTemp[0]*varDn0(4);
                    const scalar dNuTilda0 = dpwTemp[0]*(Dn0_i.dot(coefsNuTilda_[i]));
                    b += Dn0_i*Col5X1(dRhoDn0, dRhoUDn0[0], dRhoUDn0[1], dRhoUDn0[2], dRhoEDn0).transpose();
                    b_scalar += Dn0_i*dNuTilda0;
                }
                bW[i] += b;
                bN[i] += b_scalar;
            }
        }
        if (type == "wall" && fluidProps_.simulationType != "Euler")
        {
            const vectorField ownerCn = mesh_.boundary()[patchI].Cn();
            const scalarField d_ij = mag(mesh_.boundary()[patchI].delta());
            const std::vector<std::unique_ptr<Face>>& boundaryQuad = patch4stQuad_[patchI];
            forAll(bfaceCells, faceI)
            {
                const label i = bfaceCells[faceI];
                const vector dp = (0.5/d_ij[faceI])*vector(parameter_.vrWeightSqr[0], 0, 0);
                Mat9X5 b = Mat9X5::Zero();
                Col9X1 b_scalar = Col9X1::Zero();
                for (label gaussI = 0; gaussI != boundaryQuad[faceI]->size(); ++gaussI)
                {
                    const vector delta_i = boundaryQuad[faceI]->at(gaussI) - ownerCn[faceI];
                    const vector dpwTemp = boundaryQuad[faceI]->weight(gaussI)*dp;
                    Col9X1 Dn0_i = polynomialDn0(delta_i, rDeltaXYZ_[i], basisMean_[i]);
                    Col5X1 varDn0 = coefs_[i].transpose()*Dn0_i;
                    // Dn0
                    const vector rhoUDn0_i = rhoU_[i] + vector(varDn0(1), varDn0(2), varDn0(3));
                    const vector rhoUDn0_j = -rhoUDn0_i;
                    const scalar dRhoDn0  = dpwTemp[0]*varDn0(0);
                    const vector dRhoUDn0 = dpwTemp[0]*(rhoUDn0_j - rhoU_[i]);
                    const scalar dRhoEDn0 = dpwTemp[0]*varDn0(4);
                    const scalar nuTilda_i = nuTilda_[i] + Dn0_i.dot(coefsNuTilda_[i]);
                    const scalar dNuTilda0 = dpwTemp[0]*(-nuTilda_i - nuTilda_[i]);
                    b += Dn0_i*Col5X1(dRhoDn0, dRhoUDn0[0], dRhoUDn0[1], dRhoUDn0[2], dRhoEDn0).transpose();
                    b_scalar += Dn0_i*dNuTilda0;
                }
                bW[i] += b;
                bN[i] += b_scalar;
            }
        }
    }
    /* Block Gauss–Seidel iteration */
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        bW[i] += B_[faceI].transpose() * coefs_[j];
        bW[j] += B_[faceI] * coefs_[i];
        bN[i] += B_[faceI].transpose() * coefsNuTilda_[j];
        bN[j] += B_[faceI] * coefsNuTilda_[i];

    }
    forAll(mesh_.cells(), cellI)
    {
        coefs_[cellI] = rA_[cellI]*bW[cellI];
        coefsNuTilda_[cellI] = rA_[cellI]*bN[cellI];
    }
    /* Limiting coefficients of polynomials */
    limitCoefficients();
    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                d1Rho_[i]  = vector(coefs_[i](0,0), coefs_[i](1,0), coefs_[i](2,0));
                d1RhoU_[i] = tensor(coefs_[i](0,1), coefs_[i](1,1), coefs_[i](2,1),
                                    coefs_[i](0,2), coefs_[i](1,2), coefs_[i](2,2),
                                    coefs_[i](0,3), coefs_[i](1,3), coefs_[i](2,3));
                d1RhoE_[i] = vector(coefs_[i](0,4), coefs_[i](1,4), coefs_[i](2,4));
                d2Rho_[i]   = symmTensor(coefs_[i](3,0), coefs_[i](4,0), coefs_[i](5,0), coefs_[i](6,0), coefs_[i](7,0), coefs_[i](8,0));
                d2RhoUx_[i] = symmTensor(coefs_[i](3,1), coefs_[i](4,1), coefs_[i](5,1), coefs_[i](6,1), coefs_[i](7,1), coefs_[i](8,1));
                d2RhoUy_[i] = symmTensor(coefs_[i](3,2), coefs_[i](4,2), coefs_[i](5,2), coefs_[i](6,2), coefs_[i](7,2), coefs_[i](8,2));
                d2RhoUz_[i] = symmTensor(coefs_[i](3,3), coefs_[i](4,3), coefs_[i](5,3), coefs_[i](6,3), coefs_[i](7,3), coefs_[i](8,3));
                d2RhoE_[i]  = symmTensor(coefs_[i](3,4), coefs_[i](4,4), coefs_[i](5,4), coefs_[i](6,4), coefs_[i](7,4), coefs_[i](8,4));
                d1NuTilda_[i] = vector(coefsNuTilda_[i](0), coefsNuTilda_[i](1), coefsNuTilda_[i](2));
                d2NuTilda_[i] = symmTensor(coefsNuTilda_[i](3), coefsNuTilda_[i](4), coefsNuTilda_[i](5),
                                           coefsNuTilda_[i](6), coefsNuTilda_[i](7), coefsNuTilda_[i](8));
            }
        }
    }
    d1Rho_.correctBoundaryConditions();
    d1RhoU_.correctBoundaryConditions();
    d1RhoE_.correctBoundaryConditions();
    d2Rho_.correctBoundaryConditions();
    d2RhoUx_.correctBoundaryConditions();
    d2RhoUy_.correctBoundaryConditions();
    d2RhoUz_.correctBoundaryConditions();
    d2RhoE_.correctBoundaryConditions();
    d1NuTilda_.correctBoundaryConditions();
    d2NuTilda_.correctBoundaryConditions();
}

void Foam::evaluateVarsAndGrads
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    scalar& nuTilda,
    vector& rhoGrad,
    tensor& rhoUGrad,
    vector& rhoEGrad,
    vector& nuTildaGrad,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const scalar& nuTilda_0,
    const Mat9X5& coefs,
    const Col9X1& coefsNuTilda_,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
)
{
    Col9X1 Dn0 = polynomialDn0(delta, rDeltaXYZ, basisMean);
    Col5X1 varProj = coefs.transpose()*Dn0;
    rho  = rho_0  + varProj(0);
    rhoU = rhoU_0 + vector(varProj(1), varProj(2), varProj(3));
    rhoE = rhoE_0 + varProj(4);
    nuTilda = nuTilda_0 + Dn0.dot(coefsNuTilda_);
    Mat9X3 polyGrad = polynomialGrad(delta, rDeltaXYZ);
    Mat5X3 consGrad = coefs.transpose()*polyGrad;
    Col3X1 turbGrad = polyGrad.transpose()*coefsNuTilda_;
    rhoGrad  = vector(consGrad(0,0), consGrad(0,1), consGrad(0,2));
    rhoUGrad = tensor(consGrad(1,0), consGrad(1,1), consGrad(1,2),
                      consGrad(2,0), consGrad(2,1), consGrad(2,2),
                      consGrad(3,0), consGrad(3,1), consGrad(3,2));
    rhoEGrad = vector(consGrad(4,0), consGrad(4,1), consGrad(4,2));
    nuTildaGrad = vector(turbGrad(0), turbGrad(1), turbGrad(2));
}

void Foam::evaluateVarsAndGrads
(
    scalar& rho,
    vector& rhoU,
    scalar& rhoE,
    scalar& nuTilda,
    vector& rhoGrad,
    tensor& rhoUGrad,
    vector& rhoEGrad,
    vector& nuTildaGrad,
    const scalar& rho_0,
    const vector& rhoU_0,
    const scalar& rhoE_0,
    const scalar& nuTilda_0,
    const vector& d1Rho,
    const tensor& d1RhoU,
    const vector& d1RhoE,
    const vector& d1NuTilda,
    const symmTensor& d2Rho,
    const symmTensor& d2RhoUx,
    const symmTensor& d2RhoUy,
    const symmTensor& d2RhoUz,
    const symmTensor& d2RhoE,
    const symmTensor& d2NuTilda,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
)
{
    Col9X1 Dn0 = polynomialDn0(delta, rDeltaXYZ, basisMean);
    Col5X1 varProj = Mat5X9
    {
        {d1Rho[0],  d1Rho[1],  d1Rho[2],  d2Rho[0],   d2Rho[1],   d2Rho[2],   d2Rho[3],   d2Rho[4],   d2Rho[5]},
        {d1RhoU[0], d1RhoU[1], d1RhoU[2], d2RhoUx[0], d2RhoUx[1], d2RhoUx[2], d2RhoUx[3], d2RhoUx[4], d2RhoUx[5]},
        {d1RhoU[3], d1RhoU[4], d1RhoU[5], d2RhoUy[0], d2RhoUy[1], d2RhoUy[2], d2RhoUy[3], d2RhoUy[4], d2RhoUy[5]},
        {d1RhoU[6], d1RhoU[7], d1RhoU[8], d2RhoUz[0], d2RhoUz[1], d2RhoUz[2], d2RhoUz[3], d2RhoUz[4], d2RhoUz[5]},
        {d1RhoE[0], d1RhoE[1], d1RhoE[2], d2RhoE[0],  d2RhoE[1],  d2RhoE[2],  d2RhoE[3],  d2RhoE[4],  d2RhoE[5]}
    } * Dn0;
    rho  = rho_0  + varProj(0);
    rhoU = rhoU_0 + vector(varProj(1), varProj(2), varProj(3));
    rhoE = rhoE_0 + varProj(4);
    nuTilda = nuTilda_0 + Dn0.dot(Col9X1(d1NuTilda[0], d1NuTilda[1], d1NuTilda[2],
                                         d2NuTilda[0], d2NuTilda[1], d2NuTilda[2],
                                         d2NuTilda[3], d2NuTilda[4], d2NuTilda[5]));
    
    Mat9X3 polyGrad = polynomialGrad(delta, rDeltaXYZ);
    Mat5X3 consGrad = Mat5X9
    {
        {d1Rho[0],  d1Rho[1],  d1Rho[2],  d2Rho[0],   d2Rho[1],   d2Rho[2],   d2Rho[3],   d2Rho[4],   d2Rho[5]},
        {d1RhoU[0], d1RhoU[1], d1RhoU[2], d2RhoUx[0], d2RhoUx[1], d2RhoUx[2], d2RhoUx[3], d2RhoUx[4], d2RhoUx[5]},
        {d1RhoU[3], d1RhoU[4], d1RhoU[5], d2RhoUy[0], d2RhoUy[1], d2RhoUy[2], d2RhoUy[3], d2RhoUy[4], d2RhoUy[5]},
        {d1RhoU[6], d1RhoU[7], d1RhoU[8], d2RhoUz[0], d2RhoUz[1], d2RhoUz[2], d2RhoUz[3], d2RhoUz[4], d2RhoUz[5]},
        {d1RhoE[0], d1RhoE[1], d1RhoE[2], d2RhoE[0],  d2RhoE[1],  d2RhoE[2],  d2RhoE[3],  d2RhoE[4],  d2RhoE[5]}
    } * polyGrad;
    Col3X1 turbGrad = polyGrad.transpose()*Col9X1(d1NuTilda[0], d1NuTilda[1], d1NuTilda[2],
                                                  d2NuTilda[0], d2NuTilda[1], d2NuTilda[2],
                                                  d2NuTilda[3], d2NuTilda[4], d2NuTilda[5]);
    rhoGrad  = vector(consGrad(0,0), consGrad(0,1), consGrad(0,2));
    rhoUGrad = tensor(consGrad(1,0), consGrad(1,1), consGrad(1,2),
                      consGrad(2,0), consGrad(2,1), consGrad(2,2),
                      consGrad(3,0), consGrad(3,1), consGrad(3,2));
    rhoEGrad = vector(consGrad(4,0), consGrad(4,1), consGrad(4,2));
    nuTildaGrad = vector(turbGrad(0), turbGrad(1), turbGrad(2));
}