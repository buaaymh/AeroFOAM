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

#include "stabilator.H"

Foam::Stabilator::Stabilator
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    volVectorField& force
)
:
    Model(name, rho, U, force)
{
    Info << "Install stabilator model in Zone " << name << endl;
    isModified_ = mesh_.solutionDict().subDict(name).lookup<Switch>("isModified");
    isCorrected_ = mesh_.solutionDict().subDict(name).lookup<Switch>("isCorrected");
    origin_  = mesh_.solutionDict().subDict(name).lookup<vector>("origin");
    rotate_  = mesh_.solutionDict().subDict(name).lookup<vector>("rotate");
    nSpans_ = mesh_.solutionDict().subDict(name).lookup<label>("nActuatorPoints");
    refRho_ = mesh_.solutionDict().subDict(name).lookup<scalar>("referenceDensity");
    refU_ = mesh_.solutionDict().subDict(name).lookup<vector>("referenceVelocity");
    twist_ = mesh_.solutionDict().subDict(name).lookup<scalar>("twist");
    eps_ = mesh_.solutionDict().subDict(name).lookup<scalar>("gaussRadius");
    optimalPara_ = mesh_.solutionDict().subDict(name).lookup<scalar>("optimalPara");
    dSpan_ = (blade_->maxRadius() - blade_->minRadius()) / nSpans_;
    // Section
    sectionCount_ = scalarField(nSpans_, 0);
    sectionForce_ = vectorField(nSpans_);
    sectionAOA_   = scalarField(nSpans_);
    sectionUzDes_ = scalarField(nSpans_, 0);
    sectionUzOpt_ = scalarField(nSpans_, 0);
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
    // Build SphereCell
    sectionWeight_ = scalarField(nSpans_, 0);
    vector x_unit(1, 0, 0), y_unit(0, 1, 0), z_unit(0, 0, 1);
    // Rotate for Frame
    rotateX(y_unit, z_unit, rotate_.x());
    rotateY(x_unit, z_unit, rotate_.y());
    rotateZ(x_unit, y_unit, rotate_.z());
    scalar y_value = blade_->minRadius() + 0.5*dSpan_;
    for (label sectionI = 0; sectionI < nSpans_; sectionI++)
    {
        vector coordTmp = origin_ + y_value * y_unit;
        std::vector<scalar> coord{coordTmp[0], coordTmp[1], coordTmp[2]};
        const scalar searchRadius = sqrt(sqr(3*eps_) + sqr(dSpan_));
        auto neighIds = tree_->neighborhood_indices(coord, searchRadius);
        if (!neighIds.empty())
        {
            std::vector<label> cellIs; cellIs.reserve(neighIds.size());
            std::vector<scalar> weights; weights.reserve(neighIds.size());
            for (auto& cellI : neighIds)
            {
                cellI = mesh_.cellZones()[zoneI_][cellI];
                scalar r = (mesh_.C()[cellI]-origin_)&y_unit;
                if (r >= blade_->maxRadius() || r <= blade_->minRadius()) continue;
                scalar ps = (r-y_value)/dSpan_;
                scalar pn = sqrt(magSqr(mesh_.C()[cellI]-origin_)-sqr(r));
                if ((sectionI == 0) && (ps < 0)) ps = 0;
                else if ((sectionI == nSpans_-1) && (ps > 0)) ps = 0;
                scalar weight = getGaussWeight(r, mag(ps), pn);
                if (weight > 1e-6)
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
    sectionCount_  = returnReduce(sectionCount_,  sumOp<scalarField>());
    sectionWeight_ = returnReduce(sectionWeight_, sumOp<scalarField>());
    Info << sectionWeight_ << endl;
}

void Foam::Stabilator::evaluateForce(const solver* solver)
{
    // Sample Velocity
    vectorField sectionU(nSpans_, vector::zero);
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
    scalarField G(nSpans_, 0);
    scalarField U_in(nSpans_, 0);
    sectionForce_ = vector::zero;
    sectionAOA_  = 0;
    for (const auto& [sectionI, section] : sections_)
    {
        scalar r = (section.coord - origin_)&section.y_unit;
        scalar u = (sectionU[sectionI]&section.x_unit);
        scalar w = (sectionU[sectionI]&section.z_unit) - sectionUzDes_[sectionI] + sectionUzOpt_[sectionI];
        U_in[sectionI] = sqrt(sqr(u) + sqr(w) + SMALL);
        scalar twist = twist_ + blade_->twist(r);
        sectionAOA_[sectionI] = getAngleOfAttack(u, w, twist);
        auto [Cl, Cd] = blade_->Cl_Cd(U_in[sectionI], r, sectionAOA_[sectionI]);
        auto [cos, sin] = cosSin(sectionAOA_[sectionI] - twist); // angle of inflow
        scalar Cz = Cl * cos + Cd * sin;
        scalar Cx = Cd * cos - Cl * sin;
        sectionForce_[sectionI]  = Cz*section.z_unit + Cx*section.x_unit;
        sectionForce_[sectionI] *= -0.5*refRho_*magSqr(U_in[sectionI])*blade_->chord(r)*dSpan_;
        G[sectionI] = 0.5*Cz*magSqr(U_in[sectionI])*blade_->chord(r); // G for correction
        forAll(section.projectedCells, cellI)
        {
            label i = section.projectedCells[cellI];
            force_[i] += section.weights[cellI] * sectionForce_[sectionI] / sectionWeight_[sectionI];
        }
    }
    // Correction
    if (isCorrected_)
    {
        scalarField dG(nSpans_+1, 0);
        G = returnReduce(G, sumOp<scalarField>())/sectionCount_;
        dG[0] = 2*G[0];
        dG[nSpans_] = -2*G[nSpans_-1];
        for (label sectionI = 1; sectionI < nSpans_; sectionI++)
            dG[sectionI] = G[sectionI]-G[sectionI-1];
        for (const auto& [sectionI, section] : sections_)
        {
            scalar r = (section.coord - origin_)&section.y_unit;
            scalar epsOpt = optimalPara_*blade_->chord(r);
            auto [UzDes, UzOpt] = evaluateInducedVelocity(dG, U_in[sectionI], eps_, epsOpt, sectionI);
            sectionUzDes_[sectionI] = 0.1*UzDes + 0.9*sectionUzDes_[sectionI];
            sectionUzOpt_[sectionI] = 0.1*UzOpt + 0.9*sectionUzOpt_[sectionI];
        }
    }
}

void Foam::Stabilator::write()
{
    scalar lift = 0, drag = 0;
    for (const auto& [sectionI, section] : sections_)
    { 
        lift -= (sectionForce_[sectionI]&section.z_unit)/sectionCount_[sectionI];
        drag -= (sectionForce_[sectionI]&section.x_unit)/sectionCount_[sectionI];
    }
    lift = returnReduce(lift, sumOp<scalar>());
    drag = returnReduce(drag, sumOp<scalar>());
    const scalar area = sqr(blade_->maxRadius()-blade_->minRadius())/blade_->aspectRatio();
    scalar Cl = lift/(0.5*refRho_*magSqr(refU_)*area);
    scalar Cd = drag/(0.5*refRho_*magSqr(refU_)*area);
    Info << "# " << name_
         << ", AOA[deg] = " << setprecision(4) << twist_
         << ", Cl[-] = " << setprecision(4) << Cl
         << ", Cd[-] = " << setprecision(4) << Cd << endl;

    if (isModified_) twist_ -= Cl/0.2;
    
    if (mesh_.time().outputTime())
    {
        sectionAOA_ = returnReduce(sectionAOA_, sumOp<scalarField>());
        if (Pstream::master())
        {
            sectionAOA_ /= sectionCount_;
            fileName outputDir = mesh_.time().timePath();
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr<OFstream> outputFilePtr;
            // Open the file in the newly created directory
            outputFilePtr.reset(new OFstream(outputDir/"sectionInfo.dat"));
            outputFilePtr() << "#r/R" << tab << "AOA" << endl;
            for (label pointI = 0; pointI < nSpans_; pointI++)
            {
                scalar r_R = (blade_->minRadius()+(pointI+0.5)*dSpan_)/(blade_->maxRadius()-blade_->minRadius());
                outputFilePtr() << r_R << tab << sectionAOA_[pointI] << endl;
            }
        }
    }
}

std::pair<scalar, scalar> Foam::Stabilator::evaluateInducedVelocity
(
    const scalarField& dG,
    scalar U_in,
    scalar epsDes,
    scalar epsOpt,
    label  sectionI
) const
{
    scalar UyDes = 0, UyOpt = 0;
    for (label i = 0; i < nSpans_+1; i++)
    {
        scalar dy = dSpan_*(sectionI-i+0.5);
        scalar temp = dG[i]/(4*constant::mathematical::pi*dy);
        UyDes += temp*(1-Foam::exp(-sqr(dy/epsDes)));
        UyOpt += temp*(1-Foam::exp(-sqr(dy/epsOpt)));
    }
    UyDes /= -U_in;
    UyOpt /= -U_in;
    return {UyDes, UyOpt};
}

scalar Foam::Stabilator::getGaussWeight
(
    scalar r,
    scalar ps,
    scalar pn
) const
{
    if (pn > 3*eps_ || ps >= 1) return 0.0;
    return (1-ps)*Foam::exp(-sqr(pn/eps_))/(sqr(eps_)*dSpan_*constant::mathematical::pi);
}