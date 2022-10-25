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
    eulerConsVar3rdSolver

Description
    3rd order solver for Euler flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef eulerConsVar3rdSolver_H
#define eulerConsVar3rdSolver_H

#include "vr3rdSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class eulerConsVar3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class eulerConsVar3rdSolver: public vr3rdSolver
{

private:
    // Private Member Functions
    //- Disallow default bitwise copy construct
    eulerConsVar3rdSolver(const eulerConsVar3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const eulerConsVar3rdSolver&);

    void reconstructionIter();

    void limitCoefficients();

    void evaluateFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho_L,
        const scalar& rho_R,
        const vector& rhoU_L,
        const vector& rhoU_R,
        const scalar& rhoE_L,
        const scalar& rhoE_R,
        const vector& normal
    ) const;

protected:
    //- Density grad
    volVectorField rhoGrad_;
    volVectorField rhoGradLimited_;
    //- Velocity grad
    volTensorField rhoUGrad_;
    volTensorField rhoUGradLimited_;
    //- pressure grad
    volVectorField rhoEGrad_;
    volVectorField rhoEGradLimited_;
    //- Quadratic polynomial coefficients
    volSymmTensorField d2Rho_;
    volSymmTensorField d2RhoUx_;
    volSymmTensorField d2RhoUy_;
    volSymmTensorField d2RhoUz_;
    volSymmTensorField d2RhoE_;
    volSymmTensorField d2RhoLimited_;
    volSymmTensorField d2RhoUxLimited_;
    volSymmTensorField d2RhoUyLimited_;
    volSymmTensorField d2RhoUzLimited_;
    volSymmTensorField d2RhoELimited_;

public:

    // Constructors - default
    eulerConsVar3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~eulerConsVar3rdSolver() {};

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
