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

#include "rotorADM.H"

Foam::RotorADM::RotorADM
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    volVectorField& force
)
:
    Model(name, rho, U, force)
{
    Info << "Install rotor ADM model in Zone " << name << endl;
    isCorrected_ = mesh_.solutionDict().subDict(name).lookup<Switch>("isCorrected");
    origin_  = mesh_.solutionDict().subDict(name).lookup<vector>("origin");
    rotate_  = mesh_.solutionDict().subDict(name).lookup<vector>("rotate");
    nBlades_ = mesh_.solutionDict().subDict(name).lookup<label>("nBlades");
    nSectors_ = mesh_.solutionDict().subDict(name).lookup<label>("nSectors");
    nSpans_ = mesh_.solutionDict().subDict(name).lookup<label>("nActuatorPoints");
    frequence_ = mesh_.solutionDict().subDict(name).lookup<scalar>("frequence");
    refRho_ = mesh_.solutionDict().subDict(name).lookup<scalar>("referenceDensity");
    twist_ = mesh_.solutionDict().subDict(name).lookup<scalar>("twist");
    eps_ = mesh_.solutionDict().subDict(name).lookup<scalar>("gaussRadius");
    dSpan_ = (blade_->maxRadius() - blade_->minRadius()) / nSpans_;
    dSector_ = 360.0/nSectors_;
    sigma_ = scalar(nBlades_)/scalar(nSectors_);
    degOmega_ = frequence_ * 360.0;
    radOmega_ = frequence_ * 2 * constant::mathematical::pi;
    // Section
    sectionCount_ = scalarField(nSpans_*nSectors_, 0);
    sectionWeight_ = scalarField(nSpans_*nSectors_, 0);
    sectionForce_ = vectorField(nSpans_*nSectors_);
    sectionAOA_ = scalarField(nSpans_*nSectors_);
    sectionCz_ = scalarField(nSpans_*nSectors_);
    sectionUzDes_ = scalarField(nSpans_*nSectors_, 0);
    sectionUzOpt_ = scalarField(nSpans_*nSectors_, 0);
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
    // Initialization for Actuator Disks
    vector x_unit(1, 0, 0), y_unit(0, 1, 0), z_unit(0, 0, 1);
    // Rotate for Frame
    rotateX(y_unit, z_unit, rotate_.x());
    rotateY(x_unit, z_unit, rotate_.y());
    rotateZ(x_unit, y_unit, rotate_.z());
    for (label sectorI = 0; sectorI < nSectors_; sectorI++)
    {
        scalar y_value = blade_->minRadius() + 0.5*dSpan_;
        for (label pointI = 0; pointI < nSpans_; pointI++)
        {
            vector coordTmp = origin_ + y_value * y_unit;
            std::vector<scalar> coord{coordTmp[0], coordTmp[1], coordTmp[2]};
            const scalar searchRadius = sqrt(sqr(3*eps_)+sqr(dSpan_));
            auto neighIds = tree_->neighborhood_indices(coord, searchRadius);
            if (!neighIds.empty())
            {
                const size_t sectionI = nSpans_*sectorI + pointI;
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
                    scalar weight = getGaussWeight(mag(ps), pn);
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
                sectionCount_[sectionI] = 1;
            }
            y_value += dSpan_;
        }
        rotateZ(x_unit, y_unit, dSector_);
    }
    sectionCount_  = returnReduce(sectionCount_,  sumOp<scalarField>());
    sectionWeight_ = returnReduce(sectionWeight_, sumOp<scalarField>());
}

void Foam::RotorADM::evaluateForce(const solver* solver)
{
    // Sample rho and U
    vectorField sectionU(nSpans_*nSectors_, vector::zero);
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
    scalarField G(nSpans_*nSectors_, 0);
    scalarField U_in(nSpans_*nSectors_, 0);
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
        scalar w = (velocity_rel&section.z_unit) - sectionUzDes_[sectionI];
        U_in[sectionI] = sqrt(sqr(u) + sqr(w));
        w += sectionUzOpt_[sectionI];
        scalar twist = twist_ + blade_->twist(r);
        sectionAOA_[sectionI] = getAngleOfAttack(u, w, twist);
        auto [Cl, Cd] = blade_->Cl_Cd(U_in[sectionI], r, sectionAOA_[sectionI]);
        auto [cos, sin] = cosSin(sectionAOA_[sectionI] - twist); // angle of inflow
        sectionCz_[sectionI] = Cl * cos + Cd * sin;
        scalar Cx = Cd * cos - Cl * sin;
        sectionForce_[sectionI]  = sectionCz_[sectionI]*section.z_unit + Cx*x_unit;
        sectionForce_[sectionI] *= -0.5*refRho_*(sqr(u)+sqr(w))*blade_->chord(r)*dSpan_*sigma_;
        G[sectionI] = 0.5*sectionCz_[sectionI]*sqr(U_in[sectionI])*blade_->chord(r); // G for correction
        forAll(section.projectedCells, cellI)
        {
            label i = section.projectedCells[cellI];
            force_[i] += section.weights[cellI]*sectionForce_[sectionI]/sectionWeight_[sectionI];
        }
    }
    // Correction
    if (isCorrected_)
    {
        scalarField dG((nSpans_+2)*nSectors_, 0);
        G = returnReduce(G, sumOp<scalarField>())/sectionCount_;
        for (label sectorI = 0; sectorI < nSectors_; sectorI++)
        {
            label head_G  = sectorI*nSpans_;
            label head_dG = sectorI*(nSpans_+2);
            dG[head_dG]   = 2*G[head_G];
            dG[head_dG+1] = G[head_G+1] - G[head_G];
            dG[head_dG+nSpans_+1] = -2*G[head_G+nSpans_-1];
            dG[head_dG+nSpans_]   = G[head_G+nSpans_-1] - G[head_G+nSpans_-2];
            for (label i = 2; i < nSpans_; i++)
                dG[head_dG+i] = 0.5*(G[head_G+i]-G[head_G+i-2]);
        }
        for (const auto& [sectionI, section] : sections_)
        {
            scalar r = mag(section.coord - origin_);
            scalar epsOpt = 0.25*blade_->chord(r);
            auto [UzDes, UzOpt] = evaluateInducedVelocity(dG, U_in[sectionI], eps_, epsOpt, sectionI);
            sectionUzDes_[sectionI] = 0.1*UzDes + 0.9*sectionUzDes_[sectionI];
            sectionUzOpt_[sectionI] = 0.1*UzOpt + 0.9*sectionUzOpt_[sectionI];
        }
    }
}

void Foam::RotorADM::write()
{
    scalar thrust = 0, torque = 0;
    for (const auto& [sectionI, section] : sections_)
    { 
        scalar r = mag(section.coord - origin_);
        thrust -= (sectionForce_[sectionI]&section.z_unit)  /sectionCount_[sectionI];
        torque -= (sectionForce_[sectionI]&section.x_unit)*r/sectionCount_[sectionI];
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
            sectionAOA_ /= sectionCount_;
            sectionCz_  /= sectionCount_;
            fileName outputDir = mesh_.time().timePath();
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr<OFstream> outputFilePtr;
            // Open the file in the newly created directory
            outputFilePtr.reset(new OFstream(outputDir/(name_+".dat")));
            outputFilePtr() << "# CT   [-] = " << setprecision(4) << CT << nl
                            << "# CM   [-] = " << setprecision(4) << CM << endl;
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

scalar Foam::RotorADM::getGaussWeight
(
    scalar ps,
    scalar pn
) const
{
    if (pn > 3*eps_ || ps >= 1) return 0.0;
    return (1-ps)*Foam::exp(-sqr(pn/eps_))/(sqr(eps_)*dSpan_*constant::mathematical::pi);
}

std::pair<scalar, scalar> Foam::RotorADM::evaluateInducedVelocity
(
    const scalarField& dG,
    scalar U_in,
    scalar epsDes,
    scalar epsOpt,
    label  sectionI
) const
{
    label sectorI = sectionI/nSpans_;
    label spanI   = sectionI%nSpans_;
    scalar UyDes = 0, UyOpt = 0;
    for (label i = 0; i < nSpans_+2; i++)
    {
        if (i == spanI+1) continue;
        scalar dy = dSpan_*(spanI+1-i);
        scalar temp = dG[i+(nSpans_+2)*sectorI]/(4*constant::mathematical::pi*dy);
        UyDes += temp*(1-Foam::exp(-sqr(dy/epsDes)));
        UyOpt += temp*(1-Foam::exp(-sqr(dy/epsOpt)));
    }
    UyDes /= -U_in;
    UyOpt /= -U_in;
    UyDes *= sigma_;
    return {UyDes, UyOpt};
}