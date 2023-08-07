/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file isn't part of OpenFOAM.

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

Application
    unsteady2ndFoam

Description
    Density-based compressible flow solver based on cell-centered 2nd schemes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "solver.H"
#include "euler2ndSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run())
    {
        const scalarField dt_dv(runTime.deltaT().value()/mesh.V().field());

        scalarField rho_0(solver.rho());
        vectorField rhoU_0(solver.rhoU());
        scalarField rhoE_0(solver.rhoE());

        // Stage 1
        solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
        if (fluidProps.withSourceTerm)
        {
            actuationSource->addSourceTerms(runTime.value(), resRho, resRhoU, resRhoE);
        }
        solver.rho()  = rho_0  + resRho  * dt_dv;
        solver.rhoU() = rhoU_0 + resRhoU * dt_dv;
        solver.rhoE() = rhoE_0 + resRhoE * dt_dv;
        solver.correctFields();

        // Stage 2
        solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
        if (fluidProps.withSourceTerm)
        {
            actuationSource->addSourceTerms(runTime.value()+runTime.deltaT().value()/3.0, resRho, resRhoU, resRhoE);
        }
        solver.rho()  = 0.75 * rho_0  + 0.25 * (solver.rho()  + resRho  * dt_dv);
        solver.rhoU() = 0.75 * rhoU_0 + 0.25 * (solver.rhoU() + resRhoU * dt_dv);
        solver.rhoE() = 0.75 * rhoE_0 + 0.25 * (solver.rhoE() + resRhoE * dt_dv);
        solver.correctFields();

        // Stage 3
        solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
        if (fluidProps.withSourceTerm)
        {
            actuationSource->addSourceTerms(runTime.value()+2.0*runTime.deltaT().value()/3.0, resRho, resRhoU, resRhoE);
        }
        solver.rho()  = 1.0/3 * rho_0  + 2.0/3 * (solver.rho()  + resRho  * dt_dv);
        solver.rhoU() = 1.0/3 * rhoU_0 + 2.0/3 * (solver.rhoU() + resRhoU * dt_dv);
        solver.rhoE() = 1.0/3 * rhoE_0 + 2.0/3 * (solver.rhoE() + resRhoE * dt_dv);
        solver.correctFields();

        runTime++;
        Info<< "Time = " << runTime.value() << " s" << nl;
        if (fluidProps.withSourceTerm) actuationSource->write();
        runTime.write();
	
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

