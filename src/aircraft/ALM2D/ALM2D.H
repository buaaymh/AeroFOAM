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
    2D ALM

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef ALM2D_H
#define ALM2D_H

#include "model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

class ALM2D: public Model
{
    word simplingType_;
    bool isNormalized_, isSimplingP0_, isCorreForCd_;
    vector origin_, refU_;
    scalar refRho_, Cl_, Cd_, eps_;
    std::unique_ptr<KDTree> tree_;
    labelList projectedCells_;
    //- section info of all actuator points
    scalar sectionWeight_, sectionAOA_;
    vector sectionVelocity_, sectionForce_;

public:
    // Constructors
    ALM2D
    (
        const word& name,
        const volScalarField& rho,
        const volVectorField& U,
        volVectorField& force
    );

    // Member Functions
    void evaluateForce(const solver* solver);

    void getDefinedSamplingForce(const solver* solver);
    
    void getPointSamplingForce(const solver* solver);
    
    void getIntegralSamplingForce(const solver* solver);

    void write();

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //