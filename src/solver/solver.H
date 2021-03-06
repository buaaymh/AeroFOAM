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

#ifndef solver_H
#define solver_H

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
                           Class solver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class solver
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    solver(const solver&);
    //- Disallow default bitwise assignment
    void operator=(const solver&);

protected:
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
    //- Static temperature
    volScalarField T_;
    //- Momentum
    volVectorField rhoU_;
    //- Energy
    volScalarField rhoE_;
    //- Sound speed
    scalarField c_;
    //- Mach number
    volScalarField Ma_;
    //- Density grad
    volVectorField rhoGrad_;
    //- Velocity grad
    volTensorField UGrad_;
    //- temperature grad
    volVectorField TGrad_;
    //- Riemann solver
    std::unique_ptr<dbnsFlux> riemann_;
    //- Volume projections
    vectorField volProjections_;
    //- Local time steps
    scalarField localDtDv_;

public:
    // Constructors
    solver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~solver() {};

    // Member Functions
    const scalarField& rho()  const { return rho_.primitiveField(); }
    const vectorField& rhoU() const { return rhoU_.primitiveField(); }
    const scalarField& rhoE() const { return rhoE_.primitiveField(); }
    
    scalarField& rho()  { return rho_.primitiveFieldRef(); }
    vectorField& rhoU() { return rhoU_.primitiveFieldRef(); }
    scalarField& rhoE() { return rhoE_.primitiveFieldRef(); }
    
    const scalarField& LTS() const { return localDtDv_; }

    virtual void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) {};

    void solveFlowLinearSystem
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

    void solveFlowPseudoTimeSystem
    (
        const scalar& dt,
        const scalar& alpha,
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE,
        scalar& resL1
    );

    virtual void correctFields();

    virtual void updateLTS();

private:

    void volProjectionsInit();
    
    void evaluateMatrixLDU
    (
        scalarField& L,
        scalarField& U,
        scalarField& D
    ) const;

    void matrixVectorProduct
    (
        const scalarField& L,
        const scalarField& U,
        const scalarField& D,
        const scalarField& v_Rho,
        const vectorField& v_RhoU,
        const scalarField& v_RhoE,
        scalarField& y_Rho,
        vectorField& y_RhoU,
        scalarField& y_RhoE
    ) const;

    void precondition
    (
        const scalarField& L,
        const scalarField& U,
        const scalarField& D,
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) const;

};

/*---------------------------------------------------------------------------*\
                           Methods Declaration
\*---------------------------------------------------------------------------*/

void ApplyGivens
(
    const scalar s,
    const scalar c,
    scalar& h1,
    scalar& h2
);

scalar Sign(scalar x, scalar y);

void GenerateGivens
(
    scalar& dx,
    scalar& dy,
    scalar& s,
    scalar& c
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
