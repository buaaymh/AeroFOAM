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

#ifndef turbulence2ndSolver_H
#define turbulence2ndSolver_H

#include "turbulenceSolver.H"
#include "limiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class euler2ndSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class turbulence2ndSolver: public turbulenceSolver, public Limiter
{
    using turbulenceSolver::mesh_;
    // Private Data
    volVectorField rhoGrad_;
    //- temperature grad
    volVectorField pGrad_;

    // Private Member Functions
    //- Disallow default bitwise copy construct
    turbulence2ndSolver(const turbulence2ndSolver&);
    //- Disallow default bitwise assignment
    void operator=(const turbulence2ndSolver&);

public:

    // Constructors - default
    turbulence2ndSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p,
        volScalarField& nuTilda
    );

    // Destructor - default
    virtual ~turbulence2ndSolver() {};

    // Member Functions
    virtual void evaluateFlowRes
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
