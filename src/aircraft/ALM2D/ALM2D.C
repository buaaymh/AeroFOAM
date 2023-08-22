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

#include "ALM2D.H"

Foam::ALM2D::ALM2D
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    volVectorField& force
)
:
    Model(name, rho, U, force)
{
    simplingType_ = mesh_.solutionDict().subDict(name).lookup<word>("simplingType");
    isNormalized_ = mesh_.solutionDict().subDict(name).lookup<Switch>("isNormalized");
    isSimplingP0_ = mesh_.solutionDict().subDict(name).lookup<Switch>("isSimplingP0");
    isCorreForCd_ = mesh_.solutionDict().subDict(name).lookup<Switch>("isCorreForCd");
    origin_  = mesh_.solutionDict().subDict(name).lookup<vector>("origin");
    refRho_ = mesh_.solutionDict().subDict(name).lookupOrDefault<scalar>("referenceDensity", 1.0);
    refU_ = mesh_.solutionDict().subDict(name).lookup<vector>("referenceVelocity");
    Cl_ = mesh_.solutionDict().subDict(name).lookup<scalar>("Cl");
    Cd_ = mesh_.solutionDict().subDict(name).lookup<scalar>("Cd");
    eps_ = mesh_.solutionDict().subDict(name).lookup<scalar>("eps");
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
    sectionWeight_ = 0;
    std::vector<scalar> coord{origin_[0], origin_[1], origin_[2]};
    auto neighIds = tree_->neighborhood_indices(coord, 3*eps_);
    if (!neighIds.empty())
    {
        for (auto& cellI : neighIds) cellI = mesh_.cellZones()[zoneI_][cellI];
        projectedCells_ = labelList(neighIds.begin(), neighIds.end());
    }
    forAll(projectedCells_, cellI)
    {
        label globalCellI = projectedCells_[cellI];
        scalar r2 = magSqr(mesh_.C()[globalCellI] - origin_);
        sectionWeight_ += mesh_.V()[globalCellI]*Foam::exp(-r2/sqr(eps_))/(sqr(eps_)*constant::mathematical::pi);
    }
    if (isNormalized_) sectionWeight_ = returnReduce(sectionWeight_, sumOp<scalar>());
    else sectionWeight_ = 1.0;
}

void Foam::ALM2D::evaluateForce(const solver* solver)
{
    if (simplingType_   == "defined")  getDefinedSamplingForce(solver);
    else if (simplingType_ == "point")  getPointSamplingForce(solver);
    else if (simplingType_ == "integral") getIntegralSamplingForce(solver);
    else
    {
        Info << "Error in simpling type" << nl
             << "(" << nl
             << " defined" << nl
             << " point" << nl
             << " integral" << nl
             << ")" << nl
             << endl;
    }
}

void Foam::ALM2D::getDefinedSamplingForce(const solver* solver)
{
    sectionU_ = refU_;
    sectionAOA_ = getAngleOfAttack(sectionU_.x(), sectionU_.y(), 0);
    scalar refAOA = getAngleOfAttack(refU_.x(), refU_.y(), 0);
    auto [cos, sin] = cosSin(refAOA); // angle of inflow
    scalar Cx = Cd_ * cos - Cl_ * sin;
    scalar Cy = Cl_ * cos + Cd_ * sin;
    sectionForce_ = -0.5*refRho_*magSqr(refU_)*blade_->chord(0)*vector(Cx, Cy, 0);
    forAll(projectedCells_, cellLocalI)
    {
        label cellGlobalI = projectedCells_[cellLocalI];
        scalar r2 = magSqr(mesh_.C()[cellGlobalI] - origin_);
        scalar weight = mesh_.V()[cellGlobalI]*Foam::exp(-r2/sqr(eps_))/(sqr(eps_)*constant::mathematical::pi);
        force_[cellGlobalI] += sectionForce_ * weight / sectionWeight_;
    }
}

void Foam::ALM2D::getPointSamplingForce(const solver* solver)
{
    std::vector<scalar> coord{origin_[0], origin_[1], origin_[2]};
    auto adjId = tree_->nearest_index(coord);
    label simplingCellI = mesh_.cellZones()[zoneI_][adjId];
    scalar d2 = magSqr(mesh_.C()[simplingCellI] - origin_);
    scalar minDist2 = returnReduce(d2, minOp<scalar>());
    label procNo = 0;
    if (d2 == minDist2) procNo = Pstream::myProcNo();
    procNo = returnReduce(d2, maxOp<label>());
    sectionU_ = vector::zero;
    if (procNo == Pstream::myProcNo())
    {
        if (isSimplingP0_) sectionU_ = U_[simplingCellI];
        else sectionU_ = solver->sampleVelocity(simplingCellI, origin_);
    }
    sectionU_ = returnReduce(sectionU_, sumOp<vector>());
    // Correction for Cd
    if (isCorreForCd_) sectionU_ /= (1.0 - Cd_*blade_->chord(0)/(4*sqrt(constant::mathematical::pi)*eps_));
    sectionAOA_ = getAngleOfAttack(sectionU_.x(), sectionU_.y(), 0);
    scalar refAOA = getAngleOfAttack(refU_.x(), refU_.y(), 0);
    auto [cos, sin] = cosSin(refAOA); // angle of inflow
    scalar Cx = Cd_ * cos - Cl_ * sin;
    scalar Cy = Cl_ * cos + Cd_ * sin;
    sectionForce_ = -0.5*refRho_*magSqr(refU_)*blade_->chord(0)*vector(Cx, Cy, 0);
    forAll(projectedCells_, cellLocalI)
    {
        label cellGlobalI = projectedCells_[cellLocalI];
        scalar r2 = magSqr(mesh_.C()[cellGlobalI] - origin_);
        scalar weight = mesh_.V()[cellGlobalI]*Foam::exp(-r2/sqr(eps_))/(sqr(eps_)*constant::mathematical::pi);
        force_[cellGlobalI] += sectionForce_ * weight / sectionWeight_;
    }
}

void Foam::ALM2D::getIntegralSamplingForce(const solver* solver)
{
    sectionU_ = vector::zero;
    forAll(projectedCells_, cellLocalI)
    {
        label cellGlobalI = projectedCells_[cellLocalI];
        scalar r2 = magSqr(mesh_.C()[cellGlobalI] - origin_);
        scalar weight = mesh_.V()[cellGlobalI]*Foam::exp(-r2/sqr(eps_))/(sqr(eps_)*constant::mathematical::pi);
        sectionU_ += U_[cellGlobalI] * weight;
    }
    sectionU_ = returnReduce(sectionU_, sumOp<vector>())/sectionWeight_;
    // Correction for Cd
    if (isCorreForCd_) sectionU_ /= (1.0 - Cd_*blade_->chord(0)/(4*sqrt(constant::mathematical::pi)*eps_));
    sectionAOA_ = getAngleOfAttack(sectionU_.x(), sectionU_.y(), 0);
    scalar refAOA = getAngleOfAttack(refU_.x(), refU_.y(), 0);
    auto [cos, sin] = cosSin(refAOA); // angle of inflow
    scalar Cx = Cd_ * cos - Cl_ * sin;
    scalar Cy = Cl_ * cos + Cd_ * sin;
    sectionForce_ = -0.5*refRho_*magSqr(refU_)*blade_->chord(0)*vector(Cx, Cy, 0);
    forAll(projectedCells_, cellLocalI)
    {
        label cellGlobalI = projectedCells_[cellLocalI];
        scalar r2 = magSqr(mesh_.C()[cellGlobalI] - origin_);
        scalar weight = mesh_.V()[cellGlobalI]*Foam::exp(-r2/sqr(eps_))/(sqr(eps_)*constant::mathematical::pi);
        force_[cellGlobalI] += sectionForce_ * weight / sectionWeight_;
    }
}

void Foam::ALM2D::write()
{
    scalar refAOA = getAngleOfAttack(refU_.x(), refU_.y(), 0);
    scalar error_AOA = mag(sectionAOA_ - refAOA);
    scalar error_velocity = 100*(mag(sectionU_)-mag(refU_))/mag(refU_);
    Info << "# simpling weight        [-] = " << sectionWeight_  << endl;
    Info << "# simpling AOA           [-] = " << sectionAOA_  << endl;
    Info << "# error in AOA           [-] = " << error_AOA << endl;
    Info << "# simpling velocity      [-] = " << sectionU_  << endl;
    Info << "# error in velocity      [%] = " << error_velocity << endl;
}

