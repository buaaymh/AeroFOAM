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

#ifndef navierStokes2ndSolver_H
#define navierStokes2ndSolver_H

#include "navierStokesSolver.H"
#include "limiter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class navierStokesSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class navierStokes2ndSolver: public navierStokesSolver, public Limiter
{
    using navierStokesSolver::mesh_;
    // Private Data
    volVectorField rhoGrad_;
    //- temperature grad
    volVectorField pGrad_;

    // Private Member Functions
    //- Disallow default bitwise copy construct
    navierStokes2ndSolver(const navierStokes2ndSolver&);
    //- Disallow default bitwise assignment
    void operator=(const navierStokes2ndSolver&);

public:

    // Constructors - default
    navierStokes2ndSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~navierStokes2ndSolver() {};

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
