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
    eulerSolver

Description
    Solver for Euler flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef eulerSolver_H
#define eulerSolver_H

#include <memory>
#include "fvCFD.H"
#include "emptyFvPatch.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "boundMinMax.H"
#include "dbnsFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Structure fluidProperties Declaration
\*---------------------------------------------------------------------------*/

struct fluidProperties
{
    word simulationType;

    scalar gamma,       /* Gas constant */
           AOA,         /* Angle of attack */
           Mach_inf,    /* Free-stream Mach number */
           T_inf,       /* Free-stream Temperature */
           p_inf;       /* Free-stream pressure */
};

/*---------------------------------------------------------------------------*\
                           Class eulerSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class eulerSolver
{
    // Private Data
    const fluidProperties& fluidProps_;
    //- Reference to mesh
    const fvMesh& mesh_;
    //- Normal
    surfaceVectorField normal_;
    //- Static density
    volScalarField& rho_;
    //- Velocity
    volVectorField& U_;
    //- Static pressure
    volScalarField& p_;
    //- Momentum
    volVectorField rhoU_;
    //- Energy
    volScalarField rhoE_;
    //- Static temperature
    volScalarField T_;
    //- Mach number
    volScalarField Ma_;
    //- Density grad
    volVectorField rhoGrad_;
    //- Velocity grad
    volTensorField UGrad_;
    //- pressure grad
    volVectorField pGrad_;
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
    //- Riemann solver
    std::unique_ptr<dbnsFlux> riemann_;

    // Private Member Functions
    //- Disallow default bitwise copy construct
    eulerSolver(const eulerSolver&);
    //- Disallow default bitwise assignment
    void operator=(const eulerSolver&);

public:

    // Constructors - default
    eulerSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~eulerSolver() {};

    // Member Functions
    const scalarField& rho()  const { return rho_.primitiveField(); }
    const vectorField& rhoU() const { return rhoU_.primitiveField(); }
    const scalarField& rhoE() const { return rhoE_.primitiveField(); }
    
    scalarField& rho()  { return rho_.primitiveFieldRef(); }
    vectorField& rhoU() { return rhoU_.primitiveFieldRef(); }
    scalarField& rhoE() { return rhoE_.primitiveFieldRef(); }

    void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

    void correctFields();

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
