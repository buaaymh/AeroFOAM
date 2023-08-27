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

#include "rotorACE.H"

Foam::RotorACE::RotorACE
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    volVectorField& force
)
:
    Model(name, rho, U, force)
{
    Info << "Install rotor ACE model in Zone " << name << endl;
    origin_  = mesh_.solutionDict().subDict(name).lookup<vector>("origin");
    rotate_  = mesh_.solutionDict().subDict(name).lookup<vector>("rotate");
    nBlades_ = mesh_.solutionDict().subDict(name).lookup<label>("nBlades");
    nSpans_ = mesh_.solutionDict().subDict(name).lookup<label>("nActuatorPoints");
    frequence_ = mesh_.solutionDict().subDict(name).lookup<scalar>("frequence");
    refRho_ = mesh_.solutionDict().subDict(name).lookup<scalar>("referenceDensity");
    twist_ = mesh_.solutionDict().subDict(name).lookup<scalar>("twist");
    dGrid_ = mesh_.solutionDict().subDict(name).lookup<scalar>("gridSize");
    nMax_ = mesh_.solutionDict().subDict(name).lookup<scalar>("nMax");
    dSpan_ = (blade_->maxRadius() - blade_->minRadius()) / nSpans_;
    dSector_ = 360.0/nBlades_;
    degOmega_ = frequence_ * 360.0;
    radOmega_ = frequence_ * 2 * constant::mathematical::pi;
    // Section
    sectionWeight_ = scalarField(nSpans_*nBlades_);
    sectionForce_ = vectorField(nSpans_*nBlades_);
    sectionAOA_ = scalarField(nSpans_*nBlades_);
    sectionCz_ = scalarField(nSpans_*nBlades_);
    // Build KDTree
    std::vector<std::vector<scalar>> points;
    points.reserve(mesh_.cellZones()[zoneI_].size());
    labelList cellListOfZone = mesh_.cellZones()[zoneI_];
    cellListOfZone.append(-1);
    forAll(cellListOfZone, cellI)
    {   
        vector point = mesh_.C()[cellListOfZone[cellI]];
        if (cellListOfZone[cellI] >= 0)
        {
            points.push_back(std::vector<scalar>{point[0], point[1], point[2]});
        }
    }
    tree_ = std::make_unique<KDTree>(points);
}

void Foam::RotorACE::updatePosition(scalar time)
{
    time_ = time;
    sectionWeight_ = 0;
    vector x_unit(1, 0, 0), y_unit(0, 1, 0), z_unit(0, 0, 1);
    // Rotate for Frame
    rotateX(y_unit, z_unit, rotate_.x());
    rotateY(x_unit, z_unit, rotate_.y());
    rotateZ(x_unit, y_unit, rotate_.z());
    // Rotate for Actuator Lines
    rotateZ(x_unit, y_unit, degOmega_*time_);
    sections_.clear();
    for (label bladeI = 0; bladeI < nBlades_; bladeI++)
    {
        scalar y_value = blade_->minRadius() + 0.5*dSpan_;
        for (label pointI = 0; pointI < nSpans_; pointI++)
        {
            vector coordTmp = origin_ + y_value * y_unit;
            std::vector<scalar> coord{coordTmp[0], coordTmp[1], coordTmp[2]};
            const scalar eps = gaussRadius(y_value);
            const scalar searchRadius = sqrt(sqr(3*eps) + sqr(dSpan_));
            auto neighIds = tree_->neighborhood_indices(coord, searchRadius);
            if (!neighIds.empty())
            {
                const size_t sectionI = nSpans_*bladeI + pointI;
                std::vector<label> cellIs; cellIs.reserve(neighIds.size());
                std::vector<scalar> weights; weights.reserve(neighIds.size());
                for (auto& cellI : neighIds)
                {
                    cellI = mesh_.cellZones()[zoneI_][cellI];
                    const scalar r = mag((mesh_.C()[cellI]-origin_)&y_unit);
                    if (r >= blade_->maxRadius() || r <= blade_->minRadius()) continue;
                    scalar ps = (r-y_value)/dSpan_;
                    scalar pn = sqrt(magSqr(mesh_.C()[cellI]-origin_)-sqr(r));
                    if ((pointI == 0) && (ps < 0)) ps = 0;
                    else if ((pointI == nSpans_-1) && (ps > 0)) ps = 0;
                    scalar weight = getACEGaussWeight(r, mag(ps), pn);
                    if (weight > 1e-8)
                    {
                        weight *= mesh_.V()[cellI];
                        cellIs.push_back(cellI);
                        weights.push_back(weight);
                        sectionWeight_[sectionI] += weight;
                    }
                }
                sections_[sectionI] = Section();
                sections_[sectionI].coord = coordTmp;
                sections_[sectionI].projectedCells = labelList(cellIs.begin(), cellIs.end());
                sections_[sectionI].weights = scalarList(weights.begin(), weights.end());
                sections_[sectionI].x_unit = x_unit;
                sections_[sectionI].y_unit = y_unit;
                sections_[sectionI].z_unit = z_unit;
            }
            y_value += dSpan_;
        }
        rotateZ(x_unit, y_unit, dSector_);
    }
    sectionWeight_ = returnReduce(sectionWeight_, sumOp<scalarField>());
}

void Foam::RotorACE::evaluateForce(const solver* solver)
{
    // Sample rho and U
    vectorField sectionU(nSpans_*nBlades_, vector::zero);
    for (const auto& [sectionI, section] : sections_)
    {
        forAll(section.projectedCells, cellI)
        {
            label i = section.projectedCells[cellI];
            sectionU[sectionI] += section.weights[cellI] * U_[i];
        }
    }
    sectionU = returnReduce(sectionU, sumOp<vectorField>())/sectionWeight_;
    // Evaluate force on section
    sectionForce_ = vector::zero;
    sectionAOA_ = 0;
    sectionCz_ = 0;
    for (const auto& [sectionI, section] : sections_)
    {
        scalar r = mag(section.coord - origin_);
        vector velocity_rel = sectionU[sectionI] + r*radOmega_*section.x_unit;
        vector x_unit = section.x_unit;
        if (degOmega_ < 0) x_unit = -section.x_unit;
        scalar u = velocity_rel&x_unit;
        scalar w = velocity_rel&section.z_unit;
        scalar U_in = sqrt(sqr(u) + sqr(w));
        scalar twist = twist_ + blade_->twist(r);
        sectionAOA_[sectionI] = getAngleOfAttack(u, w, twist);
        auto [Cl, Cd] = blade_->Cl_Cd(U_in, r, sectionAOA_[sectionI]);
        auto [cos, sin] = cosSin(sectionAOA_[sectionI] - twist); // angle of inflow
        sectionCz_[sectionI] = Cl * cos + Cd * sin;
        scalar Cx = Cd * cos - Cl * sin;
        sectionForce_[sectionI]  = sectionCz_[sectionI]*section.z_unit + Cx*x_unit;
        sectionForce_[sectionI] *= -0.5*refRho_*sqr(U_in)*blade_->chord(r)*dSpan_;
        forAll(section.projectedCells, cellI)
        {
            label i = section.projectedCells[cellI];
            force_[i] += section.weights[cellI]*sectionForce_[sectionI]/sectionWeight_[sectionI];
        }
    }
}

void Foam::RotorACE::write()
{
    scalarField sectionCount(nSpans_*nBlades_, 0);
    for (const auto& [sectionI, section] : sections_) sectionCount[sectionI] = 1.0;
    sectionCount = returnReduce(sectionCount, sumOp<scalarField>());
    scalar thrust = 0, torque = 0;
    for (const auto& [sectionI, section] : sections_)
    { 
        scalar r = mag(section.coord - origin_);
        thrust -= (sectionForce_[sectionI]&section.z_unit)  /sectionCount[sectionI];
        torque -= (sectionForce_[sectionI]&section.x_unit)*r/sectionCount[sectionI];
    }
    thrust = returnReduce(thrust, sumOp<scalar>());
    torque = returnReduce(torque, sumOp<scalar>());
    scalar CT = thrust/(sqr(radOmega_*blade_->maxRadius())*sqr(blade_->maxRadius()) *constant::mathematical::pi);
    scalar CM = torque/(sqr(radOmega_*blade_->maxRadius())*pow3(blade_->maxRadius())*constant::mathematical::pi);
    Info << "# ------ " << name_ << " ------ #" << nl
         << "# CT   [-] = " << setprecision(4) << CT << nl
         << "# CM   [-] = " << setprecision(4) << CM << endl;

    if (mesh_.time().outputTime())
    {
        sectionAOA_ = returnReduce(sectionAOA_, sumOp<scalarField>());
        sectionCz_  = returnReduce(sectionCz_, sumOp<scalarField>());
        if (Pstream::master())
        {
            sectionAOA_ /= sectionCount;
            sectionCz_  /= sectionCount;
            fileName outputDir = mesh_.time().timePath();
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr<OFstream> outputFilePtr;
            // Open the file in the newly created directory
            outputFilePtr.reset(new OFstream(outputDir/"sectionInfo.dat"));
            outputFilePtr() << "#r/R" << tab << "Cl" << tab << "AOA" << endl;
            for (label pointI = 0; pointI < nSpans_; pointI++)
            {
                scalar r_R = (blade_->minRadius()+(pointI+0.5)*dSpan_)/blade_->maxRadius();
                outputFilePtr() << r_R << tab
                                << sectionCz_[pointI] << tab
                                << sectionAOA_[pointI] << endl;
            }
        }
    }
}

scalar Foam::RotorACE::gaussRadius(scalar r) const
{
    scalar eps = nMax_*dGrid_;
    eps *= sqrt(1.0-sqr((2*r-blade_->maxRadius())/blade_->maxRadius()));
    return max(eps, dGrid_);
}

scalar Foam::RotorACE::getACEGaussWeight
(
    scalar r,
    scalar ps,
    scalar pn
) const
{
    scalar eps = gaussRadius(r);
    if (pn > 3*eps || ps >= 1) return 0.0;
    return (1-ps)*Foam::exp(-sqr(pn/eps))/(sqr(eps)*dSpan_*constant::mathematical::pi);
}   