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
    AUSMplusUpFlux

Description
    AUSM+up flux class, computes numerical fluxes based on AUSM+ scheme. Uses primitive 
    variables as input and gives back conservative numerical fluxes.

    Meng-Sing Liou, "A sequel to AUSM, PartII: AUSM+ -up for all speeds"
    Journal of Computational Physics 214 (2006): 137-170

Author
    Sebastian Saegeler  All rights reserved.
    Rewritten by Jiri Furst, Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef AUSMplusUpFlux_H
#define AUSMplusUpFlux_H

#include "fvCFD.H"
#include "dbnsFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class AUSMplusUpFlux Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class AUSMplusUpFlux: public dbnsFlux
{
public:

    // Constructors
    AUSMplusUpFlux() {};
    AUSMplusUpFlux(const AUSMplusUpFlux &AUSMplusUpFlux);               // override default copy constructor
    AUSMplusUpFlux & operator = (const AUSMplusUpFlux &AUSMplusUpFlux); // and assignment operator

    // Destructor - default
    virtual ~AUSMplusUpFlux() {};

    // Member Functions
    virtual void evaluateFlux
    (
        scalar& rhoFlux,
        vector& rhoUFlux,
        scalar& rhoEFlux,
        const scalar& rho_L,
        const scalar& rho_R,
        const vector& U_L,
        const vector& U_R,
        const scalar& p_L,
        const scalar& p_R,
        const vector& normal
    ) const;

private:
    // Private Data

    //- Coefficients
    scalar beta_  = 0.125;
    scalar Kp_    = 0.25;
    scalar Ku_    = 0.75;
    scalar sigma_ = 1.0;
    scalar fa_    = 1.0;
    scalar alpha_ = 3.0/16.0*(5.0*sqr(fa_) - 4.0);
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
