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
    advanced rotor ACE

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef rotorACE_H
#define rotorACE_H

#include "model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

class RotorACE: public Model
{
    bool isCorrected_;
    label  nBlades_, nSpans_;
    vector origin_, rotate_;
    scalar time_, refRho_, thrust_, torque_, twist_,
           frequence_, degOmega_, radOmega_,
           dSpan_, dSector_, eps_;
    std::unique_ptr<KDTree> tree_;
    //- section info of all actuator points
    scalarField sectionCount_;
    scalarField sectionWeight_;
    vectorField sectionForce_;
    scalarField sectionAOA_;
    scalarField sectionCz_;
    scalarField sectionUzDes_;
    scalarField sectionUzOpt_;
    //- section of all actuator points
    std::map<size_t, Section> sections_;

    scalar getACEGaussWeight(scalar r, scalar ps, scalar pn) const;

    std::pair<scalar, scalar> evaluateInducedVelocity
    (
        const scalarField& dG,
        scalar U_in,
        scalar epsDes,
        scalar epsOpt,
        label  sectionI
    ) const;

public:
    // Constructors
    RotorACE
    (
        const word& name,
        const volScalarField& rho,
        const volVectorField& U,
        volVectorField& force
    );

    // Member Functions
    void updatePosition(scalar time);

    void evaluateForce(const solver* solver);

    void write();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //