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

#include "source.H"

Foam::Source::Source
(
    volScalarField& rho,
    volVectorField& U
)
:
    mesh_(rho.mesh()),
    rho_(rho),
    U_(U),
    force_
    (
        IOobject
        (
            "force",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, vector::zero)
    )
{
    forAll(mesh_.cellZones(), zoneI)
    {
        const word name = mesh_.cellZones()[zoneI].name();
        if (mesh_.solutionDict().isDict(name))
        {
            word model = mesh_.solutionDict().subDict(name).lookup<word>("model");
            if (model == "rotorALM")
            {
                models_.emplace_back(std::make_unique<RotorALM>(name, rho_, U_, force_));
            }
            else if (model == "rotorACE")
            {
                models_.emplace_back(std::make_unique<RotorACE>(name, rho_, U_, force_));
            }
            else if (model == "rotorADM")
            {
                models_.emplace_back(std::make_unique<RotorADM>(name, rho_, U_, force_));
            }
            else if (model == "rotorADE")
            {
                models_.emplace_back(std::make_unique<RotorADE>(name, rho_, U_, force_));
            }
            else if (model == "wingALM")
            {
                models_.emplace_back(std::make_unique<WingALM>(name, rho_, U_, force_));
            }
            else if (model == "wingACE")
            {
                models_.emplace_back(std::make_unique<WingACE>(name, rho_, U_, force_));
            }
            else if (model == "ALM2D")
            {
                models_.emplace_back(std::make_unique<ALM2D>(name, rho_, U_, force_));
            }
            else
            {
                Info << "Error in model type" << nl
                     << "(" << nl
                     << " rotorALM" << nl
                     << " rotorACE" << nl
                     << " rotorADM" << nl
                     << " wingALM"  << nl
                     << " wingACE"  << nl
                     << " ALM2D"  << nl
                     << ")" << nl
                     << endl;
            }
        }
    }
}

void Foam::Source::updatePosition(scalar time)
{
    for (auto& model : models_)
    {
        model->updatePosition(time);
    }
}

void Foam::Source::evaluateForce(const solver* solver)
{
    force_.primitiveFieldRef() = vector::zero; 
    for (auto& model : models_)
    { 
        model->evaluateForce(solver);
    }
}

void Foam::Source::write()
{
    if (mesh_.time().outputTime())
    {
        force_.primitiveFieldRef() /= mesh_.V();
    }
    for (auto& model : models_) { model->write(); }
    Info << "----------------------------------------" << endl;
}

