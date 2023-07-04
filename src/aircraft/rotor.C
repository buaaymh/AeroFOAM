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
    const fvMesh& mesh,
    const Blade& blade
)
:
    name_(name),
    mesh_(mesh),
    blade_(blade),
    t_current_(0.0)
{
    origin_  = mesh_.solution().subDict(name).lookup<vector>("origin");
    rotate_  = mesh_.solution().subDict(name).lookup<vector>("rotate");
    nBlades_ = mesh_.solution().subDict(name).lookup<label>("nBlades");
    nSections_  = mesh_.solution().subDict(name).lookup<label>("nSections");
    frequence_ = mesh_.solution().subDict(name).lookup<scalar>("frequence");
    sigma_ = scalar(nBlades_)/scalar(nSections_);
    degOmega_ = frequence_ * 360.0;
    radOmega_ = frequence_ * 2 * constant::mathematical::pi;
    coords_ = vectorField(blade_.numOfSpans()*nSections_);
    scalarField minDist2Local(blade_.numOfSpans()*nSections_, GREAT);
    // Initialize Sections
    procNo_ = labelList(blade_.numOfSpans()*nSections_, 0);
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
        {
            points.push_back(std::vector<scalar>{point[0], point[1], point[2]});
            quads_[cellListOfZone[cellI]] = Cell(mesh_, cellListOfZone[cellI]);
        }
    }
    tree_ = std::make_unique<KDTree>(points);

    scalar dAngle = 360.0/nSections_;
    vector x_unit(1, 0, 0), y_unit(0, 1, 0), z_unit(0, 0, 1);
    rotateX(y_unit, z_unit, rotate_.x());
    rotateY(x_unit, z_unit, rotate_.y());
    rotateZ(x_unit, y_unit, rotate_.z());
    label i = 0;
    for (label lineI = 0; lineI < nSections_; lineI++)
    {
        scalar y_value = blade_.minRadius() + 0.5 * blade_.dSpan();
        for (label pointI = 0; pointI < blade_.numOfSpans(); pointI++)
        {
            coords_[i] = origin_ + y_value * y_unit;
            std::vector<scalar> coord{coords_[i][0], coords_[i][1], coords_[i][2]};
            auto neighIds = tree_->neighborhood_indices(coord, 2.15*blade_.eps0());
            if (!neighIds.empty())
            {
                for (auto& id : neighIds) id = mesh_.cellZones()[zoneI_][id];
                auto adjId = tree_->nearest_index(coord);
                sections_[i] = Section();
                sections_[i].adjCell = mesh_.cellZones()[zoneI_][adjId];
                sections_[i].y_value = y_value;
                sections_[i].eps = blade_.GaussianRadius(y_value);
                sections_[i].projectedCells = labelList(neighIds.begin(), neighIds.end());
                sections_[i].x_unit = x_unit; sections_[i].y_unit = y_unit; sections_[i].z_unit = z_unit;
                minDist2Local[i] = magSqr(mesh_.C()[sections_[i].adjCell] - coords_[i]);
            }
            y_value += blade_.dSpan();
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
)
{
    velocity_rel -= section.y_value * radOmega_ * (-section.x_unit);
    vector x_unit = section.x_unit;
    if (degOmega_ < 0) x_unit = -section.x_unit;
    scalar u = velocity_rel&x_unit;
    scalar w = velocity_rel&section.z_unit;
    scalar deg = getAngleOfAttack(u, w);
    scalar c_l = blade_.Cl(section.y_value, deg), c_d = blade_.Cd(section.y_value, deg);
    deg -= blade_.twist();  // angle of inflow
    auto con_sin = cosSin(deg);
    scalar c_z = c_l * con_sin.first + c_d * con_sin.second;
    scalar c_x = c_d * con_sin.first - c_l * con_sin.second;
    vector force = c_z * section.z_unit + c_x * x_unit;
    force *= -0.5 * rho * (sqr(u) + sqr(w)) * (blade_.chord()*blade_.dSpan()*sigma_);
    thrust_ -= force&section.z_unit;
    torque_ -= (force&x_unit)*section.y_value;
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
    deg += blade_.twist();
    // deg := angle of attack
    if (deg < -180 || 180 < deg) {
      deg += (deg < 0 ? 360 : -360);
    }
    // deg in [-180, 180]
    return deg;
}

scalar Foam::Rotor::getProjectedWeight
(
    scalar d2,
    scalar eps
) const
{
    return Foam::exp(-d2/sqr(eps))/(pow3(eps)*sqrt(pow3(constant::mathematical::pi)));
}

void Foam::Rotor::updateSections(scalar time)
{
    t_current_ = time;
    sections_.clear();
    procNo_ = labelList(blade_.numOfSpans()*nSections_, 0);
    scalar dAngle = 360.0/nSections_;
    vector x_unit(1, 0, 0), y_unit(0, 1, 0), z_unit(0, 0, 1);
    rotateX(y_unit, z_unit, rotate_.x());
    rotateY(x_unit, z_unit, rotate_.y());
    rotateZ(x_unit, y_unit, rotate_.z());
    rotateZ(x_unit, y_unit, degOmega_*t_current_); // Rotate for Actuator Lines
    scalarField minDist2Local(blade_.numOfSpans()*nSections_, GREAT);
    label i = 0;
    for (label lineI = 0; lineI < nSections_; lineI++)
    {
        scalar y_value = blade_.minRadius() + 0.5 * blade_.dSpan();
        for (label pointI = 0; pointI < blade_.numOfSpans(); pointI++)
        {
            coords_[i] = origin_ + y_value * y_unit;
            std::vector<scalar> coord{coords_[i][0], coords_[i][1], coords_[i][2]};
            auto neighIds = tree_->neighborhood_indices(coord, 2.15*blade_.eps0());
            if (!neighIds.empty())
            {
                for (auto& id : neighIds) id = mesh_.cellZones()[zoneI_][id];
                auto adjId = tree_->nearest_index(coord);
                sections_[i] = Section();
                sections_[i].adjCell = mesh_.cellZones()[zoneI_][adjId];
                sections_[i].y_value = y_value;
                sections_[i].eps = blade_.GaussianRadius(y_value);
                sections_[i].projectedCells = labelList(neighIds.begin(), neighIds.end());
                sections_[i].x_unit = x_unit; sections_[i].y_unit = y_unit; sections_[i].z_unit = z_unit;
                minDist2Local[i] = magSqr(mesh_.C()[sections_[i].adjCell] - coords_[i]);
            }
            y_value += blade_.dSpan();
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