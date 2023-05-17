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
    euler3rdSolver

Description
    3rd order solver for Euler flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef navierStokes3rdSolver_H
#define navierStokes3rdSolver_H

#include "navierStokesSolver.H"
#include "reconstruction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class euler3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class navierStokes3rdSolver: public navierStokesSolver, public Reconstruction
{
    using navierStokesSolver::mesh_;
    using navierStokesSolver::rho_;
    using navierStokesSolver::rhoU_;
    using navierStokesSolver::rhoE_;

private:
    // Private Member Functions
    //- Disallow default bitwise copy construct
    navierStokes3rdSolver(const navierStokes3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const navierStokes3rdSolver&);
    
    // Member Functions

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
        const vector& normal,
        const scalar& magSf
    ) const;

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
        const vector& rhoGrad_L,
        const vector& rhoGrad_R,
        const tensor& rhoUGrad_L,
        const tensor& rhoUGrad_R,
        const vector& rhoEGrad_L,
        const vector& rhoEGrad_R,
        const scalar& deltaTilde,
        const vector& normal,
        const scalar& magSf
    ) const;

    void evaluateWallFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho,
        const vector& rhoU,
        const scalar& rhoE,
        const vector& rhoGrad,
        const tensor& rhoUGrad,
        const vector& rhoEGrad,
        const vector& normal,
        const scalar& magSf
    ) const;

    void evaluateInletOutletFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rhoBound,
        const vector& UBound,
        const scalar& pBound,
        const scalar& rhoInter,
        const vector& rhoUInter,
        const scalar& rhoEInter,
        const vector& rhoGrad,
        const tensor& rhoUGrad,
        const vector& rhoEGrad,
        const vector& normal,
        const scalar& magSf
    ) const;

public:

    // Constructors - default
    navierStokes3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~navierStokes3rdSolver() {};

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
