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
    solver

Description
    Encompasses variables and functions related to the flow solver.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef eulerSolver_H
#define eulerSolver_H

#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class eulerSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class eulerSolver: public solver
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    eulerSolver(const eulerSolver&);
    //- Disallow default bitwise assignment
    void operator=(const eulerSolver&);

public:
    // Constructors
    eulerSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~eulerSolver() {};

    virtual void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) = 0;

    void solveFlowLinearSystemByLUSGS
    (
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE
    );

    void solveFlowLinearSystemByGMRES
    (
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE
    );

    void solveFlowPseudoTimeSystem
    (
        const scalar& dt,
        const scalar& alpha,
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE,
        scalar& resL1
    );

    virtual void correctFields();

protected:

    virtual void evaluateMatrixLDU
    (
        const scalar& alpha,
        scalarField& L,
        scalarField& U,
        scalarField& D
    ) const;

    virtual void evaluateFlowResForGMRES
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

    virtual void matrixVectorProduct
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

private:

    void precondition
    (
        const scalar& alpha,
        const scalarField& L,
        const scalarField& U,
        const scalarField& D,
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE,
        scalarField& dRho,
        vectorField& dRhoU,
        scalarField& dRhoE
    ) const;

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
