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

protected:

    void matrixVectorProduct
    (
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE,
        const scalarField& v_Rho,
        const vectorField& v_RhoU,
        const scalarField& v_RhoE,
        scalarField& y_Rho,
        vectorField& y_RhoU,
        scalarField& y_RhoE
    );

    void evaluateFlowRes
    (
        const volScalarField& rho,
        const volVectorField& U,
        const volScalarField& p,
        const volScalarField& T,
        const volVectorField& rhoGrad,
        const volTensorField& UGrad,
        const volVectorField& TGrad,
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

private:

    void limitGrad
    (
        const volScalarField& rho,
        const volVectorField& U,
        const volScalarField& T,
        const volVectorField& rhoGrad,
        const volTensorField& UGrad,
        const volVectorField& TGrad,
        volScalarField& rhoLimit,
        volVectorField& ULimit,
        volScalarField& TLimit
    );

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
