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
    rotor

Description
    Encompasses variables and functions related to the flow solver.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef actuationSource_H
#define actuationSource_H

#include "rotor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class ActuationSource
{

public:
    //- Reference to mesh
    const fvMesh& mesh_;
    //- Density
    volScalarField& rho_;
    //- Velocity
    volVectorField& U_;
    
    volVectorField ForceSource_;

    //- second-invariant
    volScalarField Q_;

    std::unique_ptr<Blade> blade_;

    std::vector<Rotor> rotors_;

    // Constructors
    ActuationSource
    (
        volScalarField& rho,
        volVectorField& U
    );

    // Member Functions
    void addSourceTerms
    (
        scalar time,
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    );
    
    void write();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //