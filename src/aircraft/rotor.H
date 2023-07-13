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

#ifndef rotor_H
#define rotor_H

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
    label  adjCell;
    scalar y_value, eps;
    labelList projectedCells;
    vector x_unit, y_unit, z_unit;
};

struct SpanInfo
{
    scalar Fz;
    scalar Fx;
    scalar AOA;
    vector force;
};

class Rotor
{

public:
    //- Reference to mesh
    word name_;
    const fvMesh& mesh_;
    const Blade& blade_;
    label zoneI_;
    label nBlades_, nSections_;
    scalar t_current_, thrust_, torque_;
    vector origin_, rotate_;
    scalar frequence_, degOmega_, radOmega_;
    scalar sigma_;
    std::unique_ptr<KDTree> tree_;
    //- processor No of all actuator points
    labelList procNo_;
    //- coordinates of all actuator points
    vectorField coords_;
    //- force of all actuator points
    vectorField force_;
    
    std::map<size_t, Section> sections_;
    //- Gauss cell
    std::map<size_t, Cell> quads_;

    std::map<size_t, SpanInfo> spanInfo_;

    // Constructors
    Rotor
    (
        const word& name,
        const fvMesh& mesh,
        const Blade& blade
    );

    Rotor(const Rotor&) = delete;
    Rotor(Rotor&&) = default;

    // Destructor - default
    ~Rotor() noexcept = default;

    // Member Functions
    void updateSections(scalar time);

    SpanInfo getForce(scalar rho, vector velocity_rel, const Section& section);

    scalar getAngleOfAttack(scalar u, scalar w) const;

    scalar getProjectedWeight(scalar d2, scalar eps) const;

};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //