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

#include "actuationSource.H"

Foam::ActuationSource::ActuationSource
(
    volScalarField& rho,
    volVectorField& U
)
:
    mesh_(rho.mesh()),
    rho_(rho),
    U_(U),
    rhoGrad_
    (
        IOobject
        (
            "rhoGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(dimless/dimLength, vector::zero)
    ),
    UGrad_
    (
        IOobject
        (
            "UGrad",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedTensor(dimless/dimLength, tensor::zero)
    ),
    ForceSource_
    (
        IOobject
        (
            "ForceSource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector(dimless, vector::zero)
    ),
    EnergySource_
    (
        IOobject
        (
            "EnergySource",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0.0)
    ),
    Q_
    (
        IOobject
        (
            "Q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless/dimArea, 0.0)
    )
{
    blade_ = std::make_unique<Rectangle<SC1095>>(mesh_);
    forAll(mesh_.cellZones(), zoneI)
    {
        const word name = mesh_.cellZones()[zoneI].name();
        if (mesh_.solution().isDict(name))
        {
            Pout << "# Install a rotor in Zone " << name << endl;
            rotors_.emplace_back(name, mesh_, *blade_);
        }
    }
}

void Foam::ActuationSource::addSourceTerms
(
    scalar time,
    scalarField& resRho,
    vectorField& resRhoU,
    scalarField& resRhoE
)
{
    ForceSource_ = vector::zero;
    EnergySource_ = 0.0;
    rhoGrad_ = fvc::grad(rho_);
    UGrad_ = fvc::grad(U_);
    for (auto& rotor : rotors_)
    {
        if (mag(time - rotor.t_current_) > 1e-10) rotor.updateSections(time);
        rotor.force_ = vectorField(rotor.procNo_.size(), vector::zero);
        rotor.energy_ = scalarField(rotor.procNo_.size(), 0.0);
        rotor.thrust_ = 0; rotor.torque_ = 0;
        for (const auto& [pointI, section] : rotor.sections_)
        {
            if (rotor.procNo_[pointI] == Pstream::myProcNo())
            {
                const label i = section.adjCell;
                const vector delta = rotor.coords_[pointI] - mesh_.C()[i];
                const scalar rho = rho_[i] + (rhoGrad_[i]&delta);
                const vector U   = U_[i]   + (UGrad_[i]&delta);
                rotor.force_[pointI] = rotor.getForce(rho, U, section);
                rotor.energy_[pointI] = rotor.force_[pointI]&U;
            }
        }
        rotor.force_ = returnReduce(rotor.force_, sumOp<vectorField>());
        rotor.energy_ = returnReduce(rotor.energy_, sumOp<scalarField>());
        reduce(rotor.thrust_, sumOp<scalar>());
        reduce(rotor.torque_, sumOp<scalar>());
        for (const auto& [pointI, section] : rotor.sections_)
        {
            forAll(section.projectedCells, cellI)
            {
                label i = section.projectedCells[cellI];
                const Cell& quad = rotor.quads_[i];
                for (label gaussI = 0; gaussI < quad.size(); gaussI++)
                {
                    const scalar d2 = magSqr(rotor.coords_[pointI] - quad.at(gaussI));
                    const scalar weight = rotor.getProjectedWeight(d2, section.eps);
                    ForceSource_[i]  += rotor.force_[pointI] *weight*quad.weight(gaussI);
                    EnergySource_[i] += rotor.energy_[pointI]*weight*quad.weight(gaussI);
                }
            }
        }
    }
    resRhoU += ForceSource_.primitiveField();
    resRhoE += EnergySource_.primitiveField();
}

void Foam::ActuationSource::write()
{
    if (mesh_.time().outputTime())
    {
        ForceSource_.primitiveFieldRef()  /= mesh_.V();
        EnergySource_.primitiveFieldRef() /= mesh_.V();
        Q_ = 0.5*(sqr(tr(UGrad_)) - tr(UGrad_&UGrad_));
        for (const auto& rotor : rotors_)
        {
            scalar CT = 2.0*rotor.thrust_/(sqr(rotor.radOmega_*blade_->maxRadius())*sqr(blade_->maxRadius())*constant::mathematical::pi);
            scalar CM = 2.0*rotor.torque_/(sqr(rotor.radOmega_*blade_->maxRadius())*pow3(blade_->maxRadius())*constant::mathematical::pi);
            Info << "# ------ " << rotor.name_ << " ------ #" << nl
                 << "# CT   [-] = " << setprecision(4) << CT << nl
                 << "# CM   [-] = " << setprecision(4) << CM << endl;
        }
        Info << "----------------------------------------" << nl;
    }
}

