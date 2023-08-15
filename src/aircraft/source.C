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
            // else if (model == "ACM")
            // {
            //     models_.emplace_back(std::make_unique<ACM>(name, *this));
            // }
            // else if (model == "ADM")
            // {
            //     models_.emplace_back(std::make_unique<ADM>(name, *this));
            // }
            // else if (model == "FIX")
            // {
            //     models_.emplace_back(std::make_unique<FIX>(name, *this));
            // }
            else
            {
                Info << "Error in model type" << nl
                     << "(" << nl
                     << " rotorALM" << nl
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

void Foam::Source::evaluateForce()
{
    force_ = vector::zero; 
    for (auto& model : models_)
    { 
        model->evaluateForce();
    }
}

void Foam::Source::write()
{
    force_.primitiveFieldRef() /= mesh_.V();
    volTensorField UGrad = fvc::grad(U_);
    Q_ = 0.5*(sqr(tr(UGrad)) - tr(UGrad&UGrad));
    for (auto& model : models_) { model->write(); }
    Info << "----------------------------------------" << endl;
}

// void Foam::Source::addSourceTerms
// (
//     scalar time,
//     scalarField& resRho,
//     vectorField& resRhoU,
//     scalarField& resRhoE
// )
// {
//     ForceSource_ = vector::zero;
//     volVectorField rhoGrad = fvc::grad(rho_);
//     volTensorField UGrad = fvc::grad(U_);
//     for (auto& rotor : rotors_)
//     {
//         if (mag(time - rotor.t_current_) > 1e-10) rotor.updateSections(time);
//         rotor.force_ = vectorField(rotor.procNo_.size(), vector::zero);
//         rotor.thrust_ = 0; rotor.torque_ = 0;
//         rotor.spanInfo_.clear();
//         for (const auto& [pointI, section] : rotor.sections_)
//         {
//             if (rotor.procNo_[pointI] == Pstream::myProcNo())
//             {
//                 const label i = section.adjCell;
//                 const vector delta = rotor.coords_[pointI] - mesh_.C()[i];
//                 const scalar rho = rho_[i] + (rhoGrad[i]&delta);
//                 const vector U   = U_[i]   + (UGrad[i]&delta);
//                 auto span_info = rotor.getForce(rho, U, section);
//                 if (pointI < rotor.nSpans_) rotor.spanInfo_[pointI] = span_info;
//                 rotor.force_[pointI] = span_info.force;
//             }
//         }
//         rotor.force_ = returnReduce(rotor.force_, sumOp<vectorField>());
//         reduce(rotor.thrust_, sumOp<scalar>());
//         reduce(rotor.torque_, sumOp<scalar>());
//         for (const auto& [pointI, section] : rotor.sections_)
//         {
//             forAll(section.projectedCells, cellI)
//             {
//                 label i = section.projectedCells[cellI];
//                 const scalar d2 = magSqr(rotor.coords_[pointI] - mesh_.C()[i]);
//                 const scalar weight = rotor.getProjectedWeight(d2, section.eps);
//                 ForceSource_[i] += rotor.force_[pointI]*weight*mesh_.V()[i];
//             }
            
//         }
//     }
//     resRhoU += ForceSource_.primitiveField();
// }

// void Foam::Source::write()
// {
    // if (mesh_.time().outputTime())
    // {
    //     ForceSource_.primitiveFieldRef() /= mesh_.V();
    //     volTensorField UGrad = fvc::grad(U_);
    //     Q_ = 0.5*(sqr(tr(UGrad)) - tr(UGrad&UGrad));
    //     for (const auto& rotor : rotors_)
    //     {
    //         scalarField Cl(rotor.nSpans_,0.0);
    //         for (const auto& [pointI, span_info] : rotor.spanInfo_)
    //         {
    //             if (rotor.procNo_[pointI] == Pstream::myProcNo())
    //             {
    //                 Cl[pointI] = span_info.Cl;
    //             }
    //         }
    //         reduce(Cl, sumOp<scalarField>());
    //         if (Pstream::master())
    //         {
    //             fileName outputDir = mesh_.time().timePath();
    //             mkDir(outputDir);
    //             // File pointer to direct the output to
    //             autoPtr<OFstream> outputFilePtr;
    //             // Open the file in the newly created directory
    //             outputFilePtr.reset(new OFstream(outputDir/"spanInfo.dat"));
    //             outputFilePtr() << "#r/R" << tab << "Cl" << endl;
    //             forAll(Cl, spanI)
    //             {
    //                 scalar r_R = (blade_->minRadius()+spanI*rotor.dSpan_+0.5*rotor.dSpan_)/blade_->maxRadius();
    //                 outputFilePtr() << r_R << tab << Cl[spanI] << endl;
    //             }
    //         }
    //     }
    // }
    // for (const auto& rotor : rotors_)
    // {
    //     scalar CT = rotor.thrust_/(sqr(rotor.radOmega_*blade_->maxRadius())*sqr(blade_->maxRadius())*constant::mathematical::pi);
    //     scalar CM = rotor.torque_/(sqr(rotor.radOmega_*blade_->maxRadius())*pow3(blade_->maxRadius())*constant::mathematical::pi);
    //     Info << "# ------ " << rotor.name_ << " ------ #" << nl
    //             << "# CT   [-] = " << setprecision(4) << CT << nl
    //             << "# CM   [-] = " << setprecision(4) << CM << endl;
    // }
    // Info << "----------------------------------------" << nl;
// }

