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

#include "wingALM.H"

Foam::WingALM::WingALM
(
    const word& name,
    const volScalarField& rho,
    const volVectorField& U,
    volVectorField& force
)
:
    Model(name, rho, U, force)
{
    Info << "Install wing ALM model in Zone " << name << endl;
    wingType_ = mesh_.solutionDict().subDict(name).lookup<word>("blade");
    origin_  = mesh_.solutionDict().subDict(name).lookup<vector>("origin");
    rotate_  = mesh_.solutionDict().subDict(name).lookup<vector>("rotate");
    nSpans_ = mesh_.solutionDict().subDict(name).lookup<label>("nActuatorPoints");
    refRho_ = mesh_.solutionDict().subDict(name).lookup<scalar>("referenceDensity");
    refU_ = mesh_.solutionDict().subDict(name).lookup<vector>("referenceVelocity");
    twist_ = mesh_.solutionDict().subDict(name).lookup<scalar>("twist");
    eps_ = mesh_.solutionDict().subDict(name).lookup<scalar>("kernelSize");
    dSpan_ = (blade_->maxRadius() - blade_->minRadius()) / nSpans_;
    sectionForce_ = vectorField(nSpans_);
    sectionDwDu_  = scalarField(nSpans_);
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
        auto neighIds = tree_->neighborhood_indices(coord, 3*eps_);
        if (!neighIds.empty())
        {
            std::vector<label> cellIs; cellIs.reserve(neighIds.size());
            std::vector<scalar> weights; weights.reserve(neighIds.size());
            for (auto& cellI : neighIds)
            {
                cellI = mesh_.cellZones()[zoneI_][cellI];
                scalar r = mag((mesh_.C()[cellI]-origin_)&y_unit);
                if (r >= blade_->maxRadius()) continue;
                const scalar d2 = magSqr(coordTmp - mesh_.C()[cellI]);
                const scalar weight = get3DGaussWeight(d2, eps_) * mesh_.V()[cellI];
                cellIs.push_back(cellI);
                weights.push_back(weight);
                sectionWeight_[sectionI] += weight;
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
    sectionWeight_ = returnReduce(sectionWeight_, sumOp<scalarField>());
}

void Foam::WingALM::evaluateForce(const solver* solver)
{
    if (wingType_ == "RectangularWing")  getConstCirculationForce(solver);
    else if (wingType_ == "EllipticWing")  getEllipticallyLoadedForce(solver);
    else
    {
        Info << "Error in wing type" << nl
             << "(" << nl
             << " RectangularWing" << nl
             << " EllipticWing" << nl
             << ")" << nl
             << endl;
    }
}
    
void Foam::WingALM::getConstCirculationForce(const solver* solver)
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
    sectionForce_ = vector::zero;
    sectionDwDu_  = 0;
    for (const auto& [sectionI, section] : sections_)
    {
        scalar r = (section.coord - origin_)&section.y_unit;
        scalar u = refU_&section.x_unit;
        scalar w = sectionU[sectionI]&section.z_unit;
        auto [Cl, Cd] = blade_->Cl_Cd(0, r, 0);
        sectionDwDu_[sectionI] = w / u;
        scalar sectionCir = 0.5*mag(refU_)*blade_->chord(r)*Cl;
        // angle of priori inflow 
        u = mag(refU_);
        w = (blade_->maxRadius()+r)/(sqr((blade_->maxRadius()+r)) + sqr(0.1*blade_->chord(r)))
          + (blade_->maxRadius()-r)/(sqr((blade_->maxRadius()-r)) + sqr(0.1*blade_->chord(r)));
        w *= -sectionCir/(4*constant::mathematical::pi);
        auto [cos, sin] = cosSin(getAngleOfAttack(u, w, 0));
        scalar Cz = Cl * cos + Cd * sin;
        scalar Cx = Cd * cos - Cl * sin;
        sectionForce_[sectionI]  = Cz*section.z_unit + Cx*section.x_unit;
        sectionForce_[sectionI] *= -0.5*refRho_*magSqr(refU_)*blade_->chord(r)*dSpan_;
        forAll(section.projectedCells, cellI)
        {
            label i = section.projectedCells[cellI];
            force_[i] += section.weights[cellI] * sectionForce_[sectionI] / sectionWeight_[sectionI];
        }
    }
}

void Foam::WingALM::getEllipticallyLoadedForce(const solver* solver)
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
    sectionForce_ = vector::zero;
    sectionDwDu_  = 0;
    for (const auto& [sectionI, section] : sections_)
    {
        scalar r = (section.coord - origin_)&section.y_unit;
        scalar u = refU_&section.x_unit;
        scalar w = sectionU[sectionI]&section.z_unit;
        sectionDwDu_[sectionI] = w / u;
        scalar twist = twist_ + blade_->twist(r);
        scalar AOA = getAngleOfAttack(u, w, twist);
        auto [Cl, Cd] = blade_->Cl_Cd(0, r, AOA);
        auto [cos, sin] = cosSin(AOA - twist); // angle of inflow
        scalar Cz = Cl * cos + Cd * sin;
        scalar Cx = Cd * cos - Cl * sin;
        sectionForce_[sectionI]  = Cz*section.z_unit + Cx*section.x_unit;
        sectionForce_[sectionI] *= -0.5*refRho_*magSqr(refU_)*blade_->chord(r)*dSpan_;
        forAll(section.projectedCells, cellI)
        {
            label i = section.projectedCells[cellI];
            force_[i] += section.weights[cellI] * sectionForce_[sectionI] / sectionWeight_[sectionI];
        }
    }
}

void Foam::WingALM::write()
{
    scalarField sectionCount(nSpans_, 0);
    for (const auto& [sectionI, section] : sections_) sectionCount[sectionI] = 1.0;
    sectionCount = returnReduce(sectionCount, sumOp<scalarField>());
    scalar lift = 0, drag = 0;
    for (const auto& [sectionI, section] : sections_)
    { 
        lift += - (sectionForce_[sectionI]&section.z_unit)/sectionCount[sectionI];
        drag += - (sectionForce_[sectionI]&section.x_unit)/sectionCount[sectionI];
    }
    lift = returnReduce(lift, sumOp<scalar>());
    drag = returnReduce(drag, sumOp<scalar>());
    const scalar area = sqr(blade_->maxRadius()-blade_->minRadius())/blade_->aspectRatio();
    scalar Cl = lift/(0.5*refRho_*magSqr(refU_)*area);
    scalar Cd = drag/(0.5*refRho_*magSqr(refU_)*area);
    Info << "# ------ " << name_ << " ------ #" << nl
         << "# Cl   [-] = " << setprecision(4) << Cl << nl
         << "# Cd   [-] = " << setprecision(4) << Cd << endl;

    if (mesh_.time().outputTime())
    {
        sectionDwDu_ = returnReduce(sectionDwDu_, sumOp<scalarField>());
        if (Pstream::master())
        {
            sectionDwDu_ /= sectionCount;
            fileName outputDir = mesh_.time().timePath();
            mkDir(outputDir);
            // File pointer to direct the output to
            autoPtr<OFstream> outputFilePtr;
            // Open the file in the newly created directory
            outputFilePtr.reset(new OFstream(outputDir/"sectionInfo.dat"));
            outputFilePtr() << "#r/R" << tab << "DwDu" << endl;
            for (label pointI = 0; pointI < nSpans_; pointI++)
            {
                scalar r_R = (blade_->minRadius()+(pointI+0.5)*dSpan_)/(blade_->maxRadius()-blade_->minRadius());
                outputFilePtr() << r_R << tab << sectionDwDu_[pointI] << endl;
            }
        }
    }
}
