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
    Limiter

Description
    Slope limiter for 2nd solver.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef limiter_H
#define limiter_H

#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Limiter
{
    // Private Member Functions
    //- Reference to mesh
    const fvMesh& mesh_;
    //- Disallow default bitwise copy construct
    Limiter(const Limiter&);
    //- Disallow default bitwise assignment
    void operator=(const Limiter&);

    void deltaMaxMin
    (
        const volScalarField& rho,
        const volVectorField& U,
        const volScalarField& p,
        scalarField& rhoMin,
        scalarField& rhoMax,
        vectorField& UMin,
        vectorField& UMax,
        scalarField& pMin,
        scalarField& pMax
    ) const;

protected:
    //- Density limit
    volScalarField rhoLimit_;
    //- Velocity limit
    volVectorField ULimit_;
    //- pressure limit
    volScalarField pLimit_;

public:
    Limiter(const fvMesh& mesh_);

    void limiterVenkatakrishnan
    (
        const volScalarField& rho,
        const volVectorField& U,
        const volScalarField& p,
        const volVectorField& rhoGrad,
        const volTensorField& UGrad,
        const volVectorField& pGrad
    );
    
    void limiterBarthJespersen
    (
        const volScalarField& rho,
        const volVectorField& U,
        const volScalarField& p,
        const volVectorField& rhoGrad,
        const volTensorField& UGrad,
        const volVectorField& pGrad
    );

    // Destructor - default
    virtual ~Limiter() {};

};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //