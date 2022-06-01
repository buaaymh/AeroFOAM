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

    const scalar tolerance = mesh.solutionDict().subDict("LUSGS").lookupOrDefault<scalar>("tolerance", 1e-6);
    label step = 0;

    while (runTime.run())
    {
        scalarField rho_0(solver.rho());
        vectorField rhoU_0(solver.rhoU());
        scalarField rhoE_0(solver.rhoE());

        // Stage 1
        solver.updateLTS();
        solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
        solver.rho()  = rho_0  + resRho  * solver.LTS();
        solver.rhoU() = rhoU_0 + resRhoU * solver.LTS();
        solver.rhoE() = rhoE_0 + resRhoE * solver.LTS();
        solver.correctFields();

        const scalar L1Res = gSum(mag(resRho));
        outputFilePtr() << step++ << tab << runTime.elapsedCpuTime() << tab << L1Res << endl;
        if (L1Res < tolerance) break;

        // Stage 2
        solver.updateLTS();
        solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
        solver.rho()  = 0.75 * rho_0  + 0.25 * (solver.rho()  + resRho  * solver.LTS());
        solver.rhoU() = 0.75 * rhoU_0 + 0.25 * (solver.rhoU() + resRhoU * solver.LTS());
        solver.rhoE() = 0.75 * rhoE_0 + 0.25 * (solver.rhoE() + resRhoE * solver.LTS());
        solver.correctFields();

        // Stage 3
        solver.updateLTS();
        solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
        solver.rho()  = 1.0/3 * rho_0  + 2.0/3 * (solver.rho()  + resRho  * solver.LTS());
        solver.rhoU() = 1.0/3 * rhoU_0 + 2.0/3 * (solver.rhoU() + resRhoU * solver.LTS());
        solver.rhoE() = 1.0/3 * rhoE_0 + 2.0/3 * (solver.rhoE() + resRhoE * solver.LTS());
        solver.correctFields();

        runTime++;
        Info<< "Time = " << runTime.value() << " s" << nl;
        runTime.write();
	
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

