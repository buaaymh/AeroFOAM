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
    Face

Description
    Face for 3rd order solver.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef element_H
#define element_H

#include "euler3rdSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class Face Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class Face
{
public:

    // Destructor - default
    virtual ~Face() {};

    // Member Functions
    virtual label size() const = 0;

    virtual scalar weight(label i) const = 0;
    
    virtual const vector& at(label i) const = 0;

};

class Triangle3: public Face
{
public:
    // Constructors - default
    Triangle3(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Triangle3() {};

    // Member Functions
    label size() const { return 3; };

    scalar weight(label i) const { return Tria3::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};

class Triangle4: public Face
{
public:
    // Constructors - default
    Triangle4(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Triangle4() {};

    // Member Functions
    label size() const { return 4; };

    scalar weight(label i) const { return Tria4::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};

class Quadrangle4: public Face
{
public:
    // Constructors - default
    Quadrangle4(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Quadrangle4() {};

    // Member Functions
    label size() const { return 4; };

    scalar weight(label i) const { return Quad4::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};

class Quadrangle9: public Face
{
public:
    // Constructors - default
    Quadrangle9(const fvMesh& mesh, const label& faceI);

    // Destructor - default
    virtual ~Quadrangle9() {};

    // Member Functions
    label size() const { return 9; };

    scalar weight(label i) const { return Quad9::w[i]; };
    
    const vector& at(label i) const { return quadPoints[i]; };

private:
    std::vector<vector> quadPoints;
};


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //