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

bool Foam::boundMinMax
(
    volScalarField& vsf,
    const dimensionedScalar& vsf0,
    const dimensionedScalar& vsf1
)
{
    scalar minVsf = min(vsf).value();
    scalar maxVsf = max(vsf).value();
    bool changed = false;

    if (minVsf < vsf0.value() || maxVsf > vsf1.value())
    {
        Info<< "bounding " << vsf.name()
            << ", min: " << minVsf
            << " max: " << maxVsf
            << " average: " << gAverage(vsf.internalField())
            << endl;
    }

    if (minVsf < vsf0.value())
    {
        vsf.primitiveFieldRef() = max
        (
            max
            (
                vsf.primitiveField(),
                fvc::average(max(vsf, vsf0))().primitiveField()
                *pos(vsf0.value() - vsf.primitiveField())
            ),
            vsf0.value()
        );
        Info<< "new min: " << gMin(vsf.internalField()) << endl;
        vsf.correctBoundaryConditions();
        changed = true;
    }

    if (maxVsf > vsf1.value())
    {
        vsf.primitiveFieldRef() = min
        (
            min
            (
                vsf.primitiveField(),
                fvc::average(min(vsf, vsf1))().primitiveField()
                *neg(vsf1.value() - vsf.primitiveField())
                // This is needed when all values are above max
                // HJ, 18/Apr/2009
              + pos(vsf1.value() - vsf.primitiveField())*vsf1.value()
            ),
            vsf1.value()
        );
        Info<< "new max: " << gMax(vsf.internalField()) << endl;
        vsf.correctBoundaryConditions();
        changed = true;
    }
    return changed;
}

// ************************************************************************* //
