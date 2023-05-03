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
    //- Constant
    const scalar Ma_Re_;
    const scalar S_T_;
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
    //- Momentum
    volVectorField rhoU_;
    //- Energy
    volScalarField rhoE_;
    //- Mach number
    volScalarField Ma_;
    //- nuTilde
    volScalarField nuTilde_;
    //- laminar viscosity
    volScalarField laminarViscosity_;
    //- eddy viscosity
    volScalarField eddyViscosity_;
    //- Velocity grad
    volTensorField UGrad_;
    //- Velocity grad
    volVectorField TGrad_;
    // - nuTilde grad
    volVectorField nuTildeGrad_;
    //- Sound speed
    scalarField c_;
    //- delta between ij
    vectorField delta_;
    //- Volume projections
    vectorField volProjections_;
    //- Local time steps
    scalarField localDtDv_;
    //- nuMax
    scalarField nuMax_;
    //- hybrid delta
    scalarField deltaDES_;

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
    );

    void evaluateViscousFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );

    void solveFlowLinearSystemByLUSGS
    (
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE,
        scalar& resL2
    );

    void solveFlowLinearSystemByGMRES
    (
        const scalarField& resRho,
        const vectorField& resRhoU,
        const scalarField& resRhoE,
        scalar& resL2
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

    void correctFields();

    void solveTurbulence();

protected:

    virtual void evaluateMatrixLDU
    (
        const scalar& alpha,
        scalarField& L,
        scalarField& U,
        scalarField& D
    ) const;

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

    void correctPrimitiveBoundary
    (
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p,
        volScalarField& T
    );

    void correctViscousFields();

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

    void evaluateViscousFlux
    (
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const vector& U,
        const vector& TGrad,
        const tensor& UGrad,
        const scalar& laminarViscosity,
        const scalar& eddyViscosity,
        const vector& normal
    ) const;

private:

    void volProjectionsInit();

};


template<class value, class grad>
grad AverageGrad
(
    const value& value_L,
    const value& value_R,
    const grad&  grad_L,
    const grad&  grad_R,
    const vector& delta
)
{
    const scalar l = mag(delta);
    const vector t = delta/l;
    grad averageGrad = 0.5*(grad_L+grad_R);
    averageGrad -= ((averageGrad&t) - (value_R-value_L)/l)*t;
    return averageGrad;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
