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
    wing ALM

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef wingALM_H
#define wingALM_H

#include "model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

class WingALM: public Model
{
    bool isCorrected_;
    word wingType_;
    label  nSpans_;
    vector origin_, rotate_, refU_;
    scalar refRho_, twist_, eps_, optimalPara_, dSpan_;
    std::unique_ptr<KDTree> tree_;
    //- section info of all actuator points
    scalarField sectionCount_;
    scalarField sectionWeight_;
    vectorField sectionForce_;
    scalarField sectionDwDu_;
    scalarField sectionUzDes_;
    scalarField sectionUzOpt_;
    //- section of all actuator points
    std::map<size_t, Section> sections_;

    std::pair<scalar, scalar> evaluateInducedVelocity
    (
        const scalarField& dG,
        scalar U_in,
        scalar epsDes,
        scalar epsOpt,
        label  sectionI
    ) const;
    
    scalar getGaussWeight(scalar r, scalar ps, scalar pn) const;

public:
    // Constructors
    WingALM
    (
        const word& name,
        const volScalarField& rho,
        const volVectorField& U,
        volVectorField& force
    );

    // Member Functions
    void evaluateForce(const solver* solver);
    
    void getConstCirculationForce(const solver* solver);
    
    void getEllipticallyLoadedForce(const solver* solver);

    void write();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //