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
    Model

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef model_H
#define model_H

#include "blade.H"
#include "definitions.H"
#include "element.H"
#include "KDTree.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

struct Section
{
    vector coord;
    labelList projectedCells;
    scalarList weights;
    vector x_unit, y_unit, z_unit;
};

class Model
{
    
protected:
    const word name_;
    label zoneI_;
    //- Reference to mesh
    const fvMesh& mesh_;
    //- Density
    const volScalarField& rho_;
    //- Velocity
    const volVectorField& U_;
    //- ForceSource
    volVectorField& force_;
    //- Blade
    std::unique_ptr<Blade> blade_;
    
    // Member Functions
    scalar getAngleOfAttack(scalar u, scalar w, scalar twist) const;
    
    scalar get3DGaussWeight(scalar d2, scalar eps) const;

public:
    // Constructors
    Model
    (
        const word& name,
        const volScalarField& rho,
        const volVectorField& U,
        volVectorField& force
    );

    // Destructor - default
    virtual ~Model() {};

    // Member Functions
    virtual void updatePosition(scalar time) {};

    virtual void evaluateForce() = 0;

    virtual void write() = 0;

    bool empty() const
    {
        if (mesh_.cellZones()[zoneI_].size() > 0) return false;
        else return true;
    }

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //