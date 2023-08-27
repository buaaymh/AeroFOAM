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
    Source

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef source_H
#define source_H

#include "solver.H"
#include "model.H"
#include "rotorALM.H"
#include "rotorACE.H"
#include "rotorADM.H"
#include "wingALM.H"
#include "wingACE.H"
#include "ALM2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

class Source
{

friend class Model;

private:
    //- Reference to mesh
    const fvMesh& mesh_;
    //- Density
    volScalarField& rho_;
    //- Velocity
    volVectorField& U_;
    //- ForceSource
    volVectorField force_;
    //- second-invariant
    volScalarField Q_;
    //- model
    std::vector<std::unique_ptr<Model>> models_;

public:
    // Constructors
    Source
    (
        volScalarField& rho,
        volVectorField& U
    );

    // Member Functions
    void updatePosition(scalar time);

    void evaluateForce(const solver* solver);

    const vectorField& force() const { return force_.primitiveField(); }
    
    void write();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //