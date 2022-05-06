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
    euler2ndSolver

Description
    2nd order solver for Euler flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef euler2ndSolver_H
#define euler2ndSolver_H

#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class euler2ndSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class euler2ndSolver: public solver
{
    // Private Data
    //- Density limit
    volScalarField rhoLimit_;
    //- Velocity limit
    volVectorField ULimit_;
    //- pressure limit
    volScalarField pLimit_;
    //- Density min and max
    scalarField rhoMin_;
    scalarField rhoMax_;
    //- Velocity min and max
    vectorField UMin_;
    vectorField UMax_;
    //- pressure min and max
    scalarField pMin_;
    scalarField pMax_;

    // Private Member Functions
    //- Disallow default bitwise copy construct
    euler2ndSolver(const euler2ndSolver&);
    //- Disallow default bitwise assignment
    void operator=(const euler2ndSolver&);

public:

    // Constructors - default
    euler2ndSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~euler2ndSolver() {};

    // Member Functions

    void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

private:

    void limiterInit();

    void limitGrad();

    scalar Venkat
    (
        const scalar& project,
        const scalar& deltaMin,
        const scalar& deltaMax,
        const scalar& cellVolume
    );

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
