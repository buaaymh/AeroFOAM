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

#ifndef turbulenceSolver_H
#define turbulenceSolver_H

#include "navierStokesSolver.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace SA
{
    static constexpr scalar k2 = 0.1681;
    static constexpr scalar sigma = 2.0/3.0;
    static constexpr scalar Cb1 = 0.1355;
    static constexpr scalar Cb2 = 0.622;
    static constexpr scalar Cv1_3 = 357.911;
    static constexpr scalar Cv2 = 0.7;
    static constexpr scalar Cv3 = 0.9;
    static constexpr scalar Cw1 = Cb1/k2+(1.0+Cb2)/sigma;
    static constexpr scalar Cw2 = 0.3;
    static constexpr scalar Cw3_6 = 64;
    static constexpr scalar Cn1 = 16.0;
    static constexpr scalar Ct3 = 1.2;
    static constexpr scalar Ct4 = 0.5;
    static constexpr scalar Cr1 = 1.0;
    static constexpr scalar Cr2 = 2.0;
    static constexpr scalar Cr3 = 1.0;
    static constexpr scalar constDES = 0.65;
}

class turbulenceSolver: public navierStokesSolver
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    turbulenceSolver(const turbulenceSolver&);
    //- Disallow default bitwise assignment
    void operator=(const turbulenceSolver&);

    protected:
        //- turbulence viscosity
        volScalarField& nuTilda_;
        //- laminar viscosity
        volScalarField nuLam_;
        //- turbulence viscosity
        volScalarField muTurb_;
        //- nuTilda grad
        volVectorField nuTildaGrad_;
        //- vorticity
        volScalarField Omega_;
        //- rotationFunction
        volScalarField fr1_;
        //- nuTilda residual
        scalarField dNuTilda_;
        //- wall distance
        scalarField dist_;
        //- DES distance
        scalarField maxDelta_;

public:
    // Constructors
    turbulenceSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p,
        volScalarField& nuTilda
    );

    // Destructor - default
    virtual ~turbulenceSolver() {};

    // Member Functions
    const scalarField& nuTilda()  const { return nuTilda_.primitiveFieldRef(); }

    scalarField& nuTilda()  { return nuTilda_.primitiveFieldRef(); }

    const scalarField& dNuTilda()  const { return dNuTilda_; }

    virtual void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) = 0;

    virtual void evaluateTurbRes
    (
        scalarField& resNuTilda
    );

    void solveTurbLinearSystemByLUSGS
    (
        const scalarField& resNuTilda
    );

    void solveTurbLinearSystemByGMRES
    (
        const scalarField& resNuTilda
    );
    
    virtual void solveTurbPseudoTimeSystem
    (
        const scalar& dt,
        const scalar& alpha,
        const scalarField& resNuTilda
    );
    
    virtual void correctTurbulenceFields();

    virtual void correctFields();

protected:

    virtual void evaluateMatrixLDU
    (
        const scalar& alpha,
        scalarField& L,
        scalarField& U,
        scalarField& D
    ) const;

    virtual void evaluateTurbLDU
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
        const scalar& muTurb,
        const vector& normal
    ) const;

    scalar evaluateSource
    (
        const scalar& d,
        const scalar& Omega,
        const scalar& Strain,
        const scalar& nuLam,
        const scalar& nuTilda,
        const vector& nuTildaGrad,
        const scalar& volume
    ) const;

    scalar evaluateSourceJacobian
    (
        const scalar& d,
        const scalar& Omega,
        const scalar& Strain,
        const scalar& nuLam,
        const scalar& nuTilda,
        const scalar& volume
    ) const;

private:

    void precondition
    (
        const scalarField& L,
        const scalarField& U,
        const scalarField& D,
        const scalarField& resNuTilda,
        scalarField& dNuTilda
    ) const;

    void matrixVectorProduct
    (
        const scalarField& L,
        const scalarField& U,
        const scalarField& D,
        const scalarField& v_NuTilda,
        scalarField& y_NuTilda
    ) const;
};

// Functions
template<class Type>
Type fv1(Type X)
{
    Type X3 = pow3(X);
    return X3/(X3 + SA::Cv1_3);
}

template<class Type>
Type nu(const Type& nuTilda, const Type& nuLam)
{
    Type nu;
    if (nuTilda > 0.0)
    {
        nu = nuLam + nuTilda;
    }
    else
    {
        const Type X3 = pow3(nuTilda/nuLam);
        const Type fn = (SA::Cn1 + X3)/(SA::Cn1 - X3);
        nu = nuLam + fn*nuTilda;
    }
    return nu;
}

template<class Type>
Type muTurb(const Type& nuTilda, const Type& rho, const Type& muLam)
{
    Type nuTildaCorrected = max(nuTilda, 0.0);
    Type X = nuTildaCorrected*rho/muLam;
    return rho*nuTildaCorrected*fv1(X);
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
