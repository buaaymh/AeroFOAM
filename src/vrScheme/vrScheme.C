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

#include "vrScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::vrScheme::vrScheme
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    rLengthScale_
    (
        IOobject
        (
            "rLengthScale",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless, vector::zero)
    ),
    basisConst_
    (
        IOobject
        (
            "basisConst",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedSymmTensor(dimless, symmTensor::zero)
    ),
    p0_(mesh_.nCells(), false),
    trouble_(mesh_.nCells(), false),
    rA_(mesh_.nCells(), vrScheme::Matrix::Zero()),
    B_(mesh_.nInternalFaces(), vrScheme::Matrix::Zero()),
    bRho_(mesh_.nCells(), vrScheme::Column::Zero()),
    bUx_(mesh_.nCells(), vrScheme::Column::Zero()),
    bUy_(mesh_.nCells(), vrScheme::Column::Zero()),
    bUz_(mesh_.nCells(), vrScheme::Column::Zero()),
    bP_(mesh_.nCells(), vrScheme::Column::Zero()),
    quad_(mesh_.nInternalFaces(), std::vector<vector>(4, vector::zero)),
    delta_(mag(mesh_.delta())),
    limit_(mesh_.nCells(), 0)
{
    vrWeight_ = mesh_.schemesDict().subDict("vrSchemes").lookup<vector>("weightList");
    adaptive_ = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("adaptive", false);
    positive_ = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    Info << "The vrWeight is " << vrWeight_ << nl
         << "The adaptive is " << adaptive_ << nl
         << "The positive is " << positive_ << endl;
    rLengthScaleInit();
    basisConstInit();
    matInit();
    forAll(mesh_.owner(), faceI) { gaussQuad4(faceI, quad_[faceI]); }
}

#include "initFunctions.H"

#include "quadPoints.H"

#include "basisFunc.H"

#include "WAPFunc.H"

#include "limitFunc.H"

#include "updateCoefficients.H"

// ************************************************************************* //