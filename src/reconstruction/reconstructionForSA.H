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

Class
    Reconstruction

Description
    Variation reconstruction approach.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef reconstructionForSA_H
#define reconstructionForSA_H

#include "reconstruction.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ReconstructionForSA: public Reconstruction
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    ReconstructionForSA(const ReconstructionForSA&);
    //- Disallow default bitwise assignment
    void operator=(const ReconstructionForSA&);

protected:

    volScalarField& nuTilda_;

    std::vector<Col9X1> coefsNuTilda_;
    //- Quadrature points
    std::vector<Cell> cellQuad_;
    //- First polynomial coefficients
    volVectorField d1NuTilda_;
    //- Quadratic polynomial coefficients
    volSymmTensorField d2NuTilda_;

public:
    ReconstructionForSA
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& rhoU,
        volScalarField& rhoE,
        volScalarField& nuTilda
    );

    void iterationStep();

    // Destructor - default
    virtual ~ReconstructionForSA() {};

};

void evaluateVarsAndGrads
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
);

void evaluateVarsAndGrads
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
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
