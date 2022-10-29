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

#include "eulerPrimVar3rdSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::eulerPrimVar3rdSolver::eulerPrimVar3rdSolver
(
    const fluidProperties& fluidProps,
    volScalarField& rho,
    volVectorField& U,
    volScalarField& p
)
:
    vr3rdSolver(fluidProps, rho, U, p),
    rhoGrad_
    (
        IOobject
        (
            "rhoGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    UGrad_
    (
        IOobject
        (
            "UGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor(dimless/dimLength, tensor::zero)
    ),
    TGrad_
    (
        IOobject
        (
            "TGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    d2Rho_
    (
        IOobject
        (
            "d2Rho",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2Ux_
    (
        IOobject
        (
            "d2Ux",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2Uy_
    (
        IOobject
        (
            "d2Uy",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2Uz_
    (
        IOobject
        (
            "d2Uz",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d2T_
    (
        IOobject
        (
            "d2T",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    d1Var_(mesh_.nCells(), Mat5X3::Zero()),
    d2Var_(mesh_.nCells(), Mat5X6::Zero())
{
    Info << "Ths solver is 3rd order for Euler flow." << nl
         << "Ths scheme is VR with AV." << nl
         << "====================================================" << endl << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "reconstructionIter.H"

#include "limitCoefficients.H"

#include "evaluateFlowRes.H"

void Foam::eulerPrimVar3rdSolver::evaluateVars
(
    const vector& delta,
    const scalar& rho_A,
    const vector& U_A,
    const scalar& T_A,
    const vector& rhoGrad,
    const tensor& UGrad,
    const vector& TGrad,
    const symmTensor& d2Rho,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const symmTensor& d2T,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean,
    scalar& rho,
    vector& U,
    scalar& T
)
{
    symmTensor basisPoly = Foam::basisPoly(delta, rDeltaXYZ, basisMean);
    rho = rho_A + (rhoGrad&delta) + cmptSum(cmptMultiply(d2Rho, basisPoly));
    U   = U_A   + (UGrad&delta)   + vector(cmptSum(cmptMultiply(d2Ux, basisPoly)),
                                         cmptSum(cmptMultiply(d2Uy, basisPoly)),
                                         cmptSum(cmptMultiply(d2Uz, basisPoly)));
    T = T_A + (TGrad&delta) + cmptSum(cmptMultiply(d2T, basisPoly));
}

// ************************************************************************* //
