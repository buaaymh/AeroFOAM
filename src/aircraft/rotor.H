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

#include "airfoil.H"
#include "definitions.H"
#include "KDTree.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct Section
{
    label  adjCell;
    scalar y_value;
    vector x_unit, y_unit, z_unit;
};

class Rotor
{

public:
    //- Reference to mesh
    const fvMesh& mesh_;
    label zoneI_;
    scalar t_current_;
    vector origin_, rotate_;
    scalar degOmega_, radOmega_, twist_, chord_;
    scalar interiorRadius_, exteriorRadius_;
    label nBlades_, nPoints_, nLines_;
    scalar width_, dSpan_, sigma_;
    std::unique_ptr<Airfoil> airfoil_;
    std::unique_ptr<KDTree> tree_;
    //- processor No of all actuator points
    labelList procNo_;
    //- coordinates of all actuator points
    vectorField coords_;
    //- force of all actuator points
    vectorField force_;
    //- nearest cell
    std::map<size_t, Section> sections_;
    //- cell within radius
    std::map<size_t, labelList> projectedCells_;

    // Constructors
    Rotor
    (
        const word& name,
        const fvMesh& mesh
    );

    Rotor(const Rotor&) = delete;
    Rotor(Rotor&&) = default;

    // Destructor - default
    ~Rotor() noexcept = default;

    // Member Functions
    void updateSections(scalar time);

    vector getForce(scalar rho, vector velocity_rel, const Section& section) const;

    scalar getAngleOfAttack(scalar u, scalar w) const;

    scalar getProjectedWeight(scalar d2) const;

};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //