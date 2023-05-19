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

\*---------------------------------------------------------------------------*/

#include "boundMinMax.H"
#include "volFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::positiveCorrect
(
    volScalarField& vsf
)
{
    const dimensionedScalar vsf0("boundMin", dimless, 1e-6);
    label index = 0;
    scalar minVsf = GREAT;
    forAll(vsf, cellI)
    {
        if(vsf[cellI] < minVsf)
        {
            minVsf = vsf[cellI];
            index  = cellI;
        }
    }
    bool changed = false;
    if (minVsf < vsf0.value())
    {
        Pout << "----------------------------------------" << nl;
        Pout << "# Min " << vsf.name() << " = " << minVsf  << " at " << vsf.mesh().C()[index]<< endl;
        vsf.primitiveFieldRef() = max
        (
            max
            (
                vsf.primitiveField(),
                fvc::average(max(vsf, vsf0))().primitiveField()
              * pos0(-vsf.primitiveField())
            ),
            vsf0.value()
        );
        changed = true;
        Pout << "----------------------------------------" << nl;
    }
    return changed;
}

// ************************************************************************* //
