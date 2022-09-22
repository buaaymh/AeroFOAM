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
    steadyTurbFoam

Description
    Density-based compressible flow solver for TurbSA euqations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "solver.H"
#include "turb3rdSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const scalar tolerance = mesh.solutionDict().subDict("SOLVER").lookupOrDefault<scalar>("tolerance", 1e-6);
    
    while (runTime.run())
    {
        runTime++;
        Info<< "Time = " << runTime.value() << " s" << nl;

        solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
        // solver.evaluateTurbRes(resNuTilde, L_turb, U_turb, D_turb);
        const scalar L2ResRho  = Foam::sqrt(gSumSqr(resRho));
        const scalar L2ResRhoU = Foam::sqrt(gSum(magSqr(resRhoU)));
        const scalar L2ResRhoE = Foam::sqrt(gSumSqr(resRhoE));
        Info << "Residual for rho  = " << L2ResRho << endl;
        Info << "Residual for rhoU = " << L2ResRhoU << endl;
        Info << "Residual for rhoE = " << L2ResRhoE << endl;

        if (Pstream::master())
        {
            outputFilePtr() << runTime.value() << tab
                            << runTime.elapsedCpuTime() << tab
                            << resRho << endl;
        }
        if ((L2ResRho < tolerance) && (L2ResRhoU < tolerance) && (L2ResRhoE < tolerance)) break;
        solver.solveFlowLinearSystem(resRho, resRhoU, resRhoE);
        // solver.solveTurbLinearSystem(resNuTilde, L_turb, U_turb, D_turb);
        solver.correctFields();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //

