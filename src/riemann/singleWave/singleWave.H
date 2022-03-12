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
    singleWave

Description
    Riemann solver for single wave fluxes.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef singleWave_H
#define singleWave_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class singleWave Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class singleWave
{
    // Private Data
    vector a_;

    // Private Member Functions
    //- Disallow default bitwise copy construct
    singleWave(const singleWave&);
    //- Disallow default bitwise assignment
    void operator=(const singleWave&);

public:

    // Constructors
    singleWave(const vector& a);

    // Destructor - default
    virtual ~singleWave() {};

    // Member Functions
    virtual void evaluateFlux
    (
        scalar& flux,
        const scalar& uLeft,
        const scalar& uRight,
        const vector& normal
    ) const;
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
