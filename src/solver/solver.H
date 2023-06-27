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

#include "emptyFvPatch.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "boundMinMax.H"
#include "dbnsFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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
    //- Riemann solver
    std::unique_ptr<dbnsFlux> riemann_;
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
    //- Mach number
    volScalarField Ma_;
    //- Momentum
    volVectorField rhoU_;
    //- Energy
    volScalarField rhoE_;
    //- Sound speed
    scalarField c_;
    //- Volume projections
    vectorField volProjections_;
    //- Local time steps
    scalarField localDtDv_;
    //- Continuity residual
    scalarField dRho_;
    //- Momentum residual
    vectorField dRhoU_;
    //- Energy residual
    scalarField dRhoE_;

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
    
    const scalarField& dRho()  const { return dRho_; }
    const vectorField& dRhoU() const { return dRhoU_; }
    const scalarField& dRhoE() const { return dRhoE_; }
    
    const scalarField& LTS() const { return localDtDv_; }

    void updateLTS()
    {
        scalar CFL = mesh_.solution().subDict("SOLVER").lookupOrDefault<scalar>("CFL", 1.0);
        localDtDv_ = CFL/(volProjections_&(cmptMag(U_.primitiveField())+c_*vector::one));
    }

    virtual void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) = 0;

    virtual void solveFlowLinearSystemByLUSGS
    (
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE
    ) = 0;

    virtual void solveFlowLinearSystemByGMRES
    (
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE
    ) = 0;

    virtual void solveFlowPseudoTimeSystem
    (
        const scalar& dt,
        const scalar& alpha,
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE
    ) = 0;

    virtual void correctFields();

protected:

    virtual void evaluateMatrixLDU
    (
        const scalar& alpha,
        scalarField& L,
        scalarField& U,
        scalarField& D
    ) const = 0;

    virtual void evaluateFlowResForGMRES
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) = 0;

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
    ) = 0;

private:

    void volProjectionsInit();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
