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

#ifndef euler3rdSolver_H
#define euler3rdSolver_H

#include "solver.H"
#include "reconstruction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class euler3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class euler3rdSolver: public solver, public Reconstruction
{
    using solver::mesh_;
    using solver::rho_;
    using solver::rhoU_;
    using solver::rhoE_;

private:
    // Private Member Functions
    //- Disallow default bitwise copy construct
    euler3rdSolver(const euler3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const euler3rdSolver&);
    
    // Member Functions

    void evaluateVars
    (
        scalar& rho,
        vector& rhoU,
        scalar& rhoE,
        const scalar& rho_0,
        const vector& rhoU_0,
        const scalar& rhoE_0,
        const Mat9X5& coefs,
        const symmTensor& basisMean,
        const vector& rDeltaXYZ,
        const vector& delta
    ) const;

    void evaluateVars
    (
        scalar& rho,
        vector& rhoU,
        scalar& rhoE,
        const scalar& rho_0,
        const vector& rhoU_0,
        const scalar& rhoE_0,
        const vector& d1Rho,
        const tensor& d1RhoU,
        const vector& d1RhoE,
        const symmTensor& d2Rho,
        const symmTensor& d2RhoUx,
        const symmTensor& d2RhoUy,
        const symmTensor& d2RhoUz,
        const symmTensor& d2RhoE,
        const symmTensor& basisMean,
        const vector& rDeltaXYZ,
        const vector& delta
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
        const vector& normal,
        const scalar& magSf
    ) const;

public:

    // Constructors - default
    euler3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~euler3rdSolver() {};

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
