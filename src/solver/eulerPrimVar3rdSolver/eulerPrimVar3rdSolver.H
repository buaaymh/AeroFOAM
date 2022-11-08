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
    eulerPrimVar3rdSolver

Description
    3rd order solver for Euler flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef eulerPrimVar3rdSolver_H
#define eulerPrimVar3rdSolver_H

#include "vr3rdSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class eulerPrimVar3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class eulerPrimVar3rdSolver: public vr3rdSolver
{

private:
    // Private Member Functions
    //- Disallow default bitwise copy construct
    eulerPrimVar3rdSolver(const eulerPrimVar3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const eulerPrimVar3rdSolver&);

    void limitCoefficients();

    void evaluateFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho_L,
        const scalar& rho_R,
        const vector& U_L,
        const vector& U_R,
        const scalar& p_L,
        const scalar& p_R,
        const vector& normal,
        const scalar& magSf
    ) const;

protected:
    //- Density grad
    volVectorField rhoGrad_;
    //- Velocity grad
    volTensorField UGrad_;
    //- pressure grad
    volVectorField pGrad_;
    //- Quadratic polynomial coefficients
    volSymmTensorField d2Rho_;
    volSymmTensorField d2Ux_;
    volSymmTensorField d2Uy_;
    volSymmTensorField d2Uz_;
    volSymmTensorField d2P_;
    
    void reconstructionIter();

public:

    // Constructors - default
    eulerPrimVar3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~eulerPrimVar3rdSolver() {};

    // Member Functions
    void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
