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
    navierStokesSolver

Description
    2nd order solver for Euler flow.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef navierStokesSolver_H
#define navierStokesSolver_H

#include "eulerSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class navierStokesSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class navierStokesSolver: public eulerSolver
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    navierStokesSolver(const navierStokesSolver&);
    //- Disallow default bitwise assignment
    void operator=(const navierStokesSolver&);

protected:
    //- Constant
    const scalar Ma_Re_;
    const scalar S_T_;
    //- laminar viscosity
    volScalarField muLam_;
    //- Velocity grad
    volTensorField UGrad_;
    //- Velocity grad
    volVectorField TGrad_;
    //- delta between ij
    vectorField delta_;

public:
    // Constructors
    navierStokesSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~navierStokesSolver() {};

    // Member Functions

    virtual void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) = 0;

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
    
    void evaluateViscousFlux
    (
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const vector& U,
        const vector& TGrad,
        const tensor& UGrad,
        const scalar& muLam,
        const vector& normal
    ) const;

};

template<class Type>
Type muLam
(
    const Type& T,
    const scalar& S_T
)
{ 
    return pow(T, 1.5)*(1.0+S_T)/(T+S_T);
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
