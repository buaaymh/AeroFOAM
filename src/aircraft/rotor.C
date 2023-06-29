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

#include "rotor.H"

Foam::Rotor::Rotor
(
    const word& name,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    t_current_(0.0)
{
    origin_  = mesh_.solution().subDict(name).lookup<vector>("origin");
    rotate_  = mesh_.solution().subDict(name).lookup<vector>("rotate");
    twist_   = mesh_.solution().subDict(name).lookup<scalar>("twist");
    chord_   = mesh_.solution().subDict(name).lookup<scalar>("chord");
    width_   = mesh_.solution().subDict(name).lookup<scalar>("GaussianWidth");
    nBlades_ = mesh_.solution().subDict(name).lookup<label>("nBlades");
    nPoints_ = mesh_.solution().subDict(name).lookup<label>("nPoints");
    nLines_  = mesh_.solution().subDict(name).lookup<label>("nLines");
    interiorRadius_ = mesh_.solution().subDict(name).lookup<scalar>("interiorRadius");
    exteriorRadius_ = mesh_.solution().subDict(name).lookup<scalar>("exteriorRadius");
    dSpan_ = (exteriorRadius_-interiorRadius_)/nPoints_;
    sigma_ = scalar(nBlades_)/scalar(nLines_);
    word airfoil = mesh_.solution().subDict(name).lookup<word>("airfoil");
    if (airfoil == "SC1095") airfoil_ = std::make_unique<Airfoil2D<SC1095>>();
    scalar frequence = mesh_.solution().subDict(name).lookup<scalar>("frequence");
    degOmega_ = frequence * 360.0;
    radOmega_ = frequence * 2 * constant::mathematical::pi;
    coords_ = vectorField(nPoints_*nLines_);
    scalarField minDist2Local(nPoints_*nLines_, GREAT);
    // Initialize Sections
    procNo_ = labelList(nPoints_*nLines_, 0);
    force_ = vectorField(nPoints_*nLines_, vector::zero);
        
    // Build KDTree
    zoneI_ = mesh_.cellZones().findZoneID(name);
    std::vector<std::vector<scalar>> points;
    points.reserve(mesh_.cellZones()[zoneI_].size());
    labelList cellListOfZone = mesh.cellZones()[zoneI_];
    cellListOfZone.append(-1);
    forAll(cellListOfZone, cellI)
    {   
        vector point = mesh_.C()[cellListOfZone[cellI]];
        if (cellListOfZone[cellI] >= 0)
            points.push_back(std::vector<scalar>{point[0], point[1], point[2]});
    }
    tree_ = std::make_unique<KDTree>(points);

    scalar dAngle = 360.0/nLines_;
    vector x_unit(1, 0, 0), y_unit(0, 1, 0), z_unit(0, 0, 1);
    rotateX(y_unit, z_unit, rotate_.x());
    rotateY(x_unit, z_unit, rotate_.y());
    rotateZ(x_unit, y_unit, rotate_.z());
    label i = 0;
    for (label lineI = 0; lineI < nLines_; lineI++)
    {
        scalar y_value = interiorRadius_ + 0.5 * dSpan_;
        for (label pointI = 0; pointI < nPoints_; pointI++)
        {
            coords_[i] = origin_ + y_value * y_unit;
            std::vector<scalar> coord{coords_[i][0], coords_[i][1], coords_[i][2]};
            auto neighIds = tree_->neighborhood_indices(coord, 2.63*width_);
            if (!neighIds.empty())
            {
                for (auto& id : neighIds) id = mesh_.cellZones()[zoneI_][id];
                projectedCells_[i] = labelList(neighIds.begin(), neighIds.end());
                auto adjId = tree_->nearest_index(coord);
                sections_[i] = Section();
                sections_[i].adjCell = mesh_.cellZones()[zoneI_][adjId];
                sections_[i].y_value = y_value;
                sections_[i].x_unit = x_unit; sections_[i].y_unit = y_unit; sections_[i].z_unit = z_unit;
                minDist2Local[i] = magSqr(mesh_.C()[sections_[i].adjCell] - coords_[i]);
            }
            y_value += dSpan_;
            i++;
        }
        rotateZ(x_unit, y_unit, dAngle);
    }
    scalarField minDist2Global = returnReduce(minDist2Local, minOp<scalarField>());
    forAll(minDist2Local, i)
    {
        if (minDist2Global[i] == minDist2Local[i]) procNo_[i] = Pstream::myProcNo();
    }
    procNo_ = returnReduce(procNo_, maxOp<labelList>());
}

Foam::vector Foam::Rotor::getForce
(
    scalar rho,
    vector velocity_rel,
    const Section& section
) const
{
    velocity_rel -= section.y_value * radOmega_ * (-section.x_unit);
    vector x_unit = section.x_unit;
    if (degOmega_ < 0) x_unit = -section.x_unit;
    scalar u = velocity_rel&x_unit;
    scalar w = velocity_rel&section.z_unit;
    scalar deg = getAngleOfAttack(u, w);
    scalar c_l = airfoil_->cLift(deg), c_d = airfoil_->cDrag(deg);
    deg -= twist_;  // angle of inflow
    auto con_sin = cosSin(deg);
    scalar c_z = c_l * con_sin.first + c_d * con_sin.second;
    scalar c_x = c_d * con_sin.first - c_l * con_sin.second;
    vector force = c_z * section.z_unit + c_x * x_unit;
    force *= -0.5 * rho * (sqr(u) + sqr(w)) * (chord_*dSpan_*sigma_);
    return force;
}

Foam::scalar Foam::Rotor::getAngleOfAttack
(
    scalar u,
    scalar w
) const
{
    auto deg = rad2deg(Foam::atan(w / u));  // [-90, 90]
    if (u < 0) {
      if (w > 0) {
        deg += 180/* [-90, 0] -> [90, 180] */;
      } else {
        deg -= 180/* [0, 90] -> [-180, -90] */;
      }
    }
    // deg := angle of inflow
    deg += twist_;
    // deg := angle of attack
    if (deg < -180 || 180 < deg) {
      deg += (deg < 0 ? 360 : -360);
    }
    // deg in [-180, 180]
    return deg;
}

scalar Foam::Rotor::getProjectedWeight
(
    scalar d2
) const
{
    return Foam::exp(-d2/sqr(width_))/(pow3(width_)*sqrt(pow3(constant::mathematical::pi)));
}

void Foam::Rotor::updateSections(scalar time)
{
    t_current_ = time;
    projectedCells_.clear();
    sections_.clear();
    procNo_ = labelList(nPoints_*nLines_, 0);
    force_ = vectorField(nPoints_*nLines_, vector::zero);
    scalar dAngle = 360.0/nLines_;
    vector x_unit(1, 0, 0), y_unit(0, 1, 0), z_unit(0, 0, 1);
    rotateX(y_unit, z_unit, rotate_.x());
    rotateY(x_unit, z_unit, rotate_.y());
    rotateZ(x_unit, y_unit, rotate_.z());
    rotateZ(x_unit, y_unit, degOmega_*t_current_); // Rotate for Actuator Lines
    scalarField minDist2Local(nPoints_*nLines_, GREAT);
    label i = 0;
    for (label lineI = 0; lineI < nLines_; lineI++)
    {
        scalar y_value = interiorRadius_ + 0.5 * dSpan_;
        for (label pointI = 0; pointI < nPoints_; pointI++)
        {
            coords_[i] = origin_ + y_value * y_unit;
            std::vector<scalar> coord{coords_[i][0], coords_[i][1], coords_[i][2]};
            auto neighIds = tree_->neighborhood_indices(coord, 2.63*width_);
            if (!neighIds.empty())
            {
                for (auto& id : neighIds) id = mesh_.cellZones()[zoneI_][id];
                projectedCells_[i] = labelList(neighIds.begin(), neighIds.end());
                auto adjId = tree_->nearest_index(coord);
                sections_[i] = Section();
                sections_[i].adjCell = mesh_.cellZones()[zoneI_][adjId];
                sections_[i].y_value = y_value;
                sections_[i].x_unit = x_unit; sections_[i].y_unit = y_unit; sections_[i].z_unit = z_unit;
                minDist2Local[i] = magSqr(mesh_.C()[sections_[i].adjCell] - coords_[i]);
            }
            y_value += dSpan_;
            i++;
        }
        rotateZ(x_unit, y_unit, dAngle);
    }
    scalarField minDist2Global = returnReduce(minDist2Local, minOp<scalarField>());
    forAll(minDist2Local, i)
    {
        if (minDist2Global[i] == minDist2Local[i]) procNo_[i] = Pstream::myProcNo();
    }
    procNo_ = returnReduce(procNo_, maxOp<labelList>());
}