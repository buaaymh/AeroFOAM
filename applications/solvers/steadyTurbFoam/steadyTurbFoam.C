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
#include "turbulenceSolver.H"
#include "turbulence2ndSolver.H"
#include "turbulence3rdSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const scalar tolerance = mesh.solution().subDict("SOLVER").lookupOrDefault<scalar>("tolerance", 1e-6);
    const scalar nCells    = scalar(returnReduce(mesh.nCells(), sumOp<label>()));
    const word method = mesh.solution().subDict("SOLVER").lookupOrDefault<word>("method", "LUSGS");
    const scalar CFL  = mesh.solution().subDict("SOLVER").lookupOrDefault<scalar>("CFL", 1.0);

    while (runTime.run())
    {
        runTime++;
        Info << "========================================" << nl;
        Info << "# " << method << " # Iteration Step = " << runTime.value() << endl;
        Info << "========================================" << nl;
        solver->evaluateFlowRes(resRho, resRhoU, resRhoE);
        solver->evaluateTurbRes(resNuTilda);
        Info << "# Local Courant          [-] = " << CFL << endl;
        Info << "----------------------------------------" << nl;

        if (method == "LUSGS")
        {
            solver->solveTurbLinearSystemByLUSGS(resNuTilda);
            solver->solveFlowLinearSystemByLUSGS(resRho, resRhoU, resRhoE);
        }
        else
        {
            solver->solveTurbLinearSystemByGMRES(resNuTilda);
            solver->solveFlowLinearSystemByGMRES(resRho, resRhoU, resRhoE);
        }
        solver->correctFields();

        scalar resRho  = gSum(mag(solver->dRho()))/nCells;
        scalar resRhoU = gSum(mag(solver->dRhoU()))/nCells;
        scalar resRhoE = gSum(mag(solver->dRhoE()))/nCells;
        scalar resNuTilda = gSum(mag(solver->dNuTilda()))/nCells;
        Info << "# Continuity residual    [-] = " << resRho  << endl;
        Info << "# Momentum   residual    [-] = " << resRhoU << endl;
        Info << "# Energy     residual    [-] = " << resRhoE << endl;
        Info << "# Turbulence residual    [-] = " << resNuTilda << endl;
        Info << "----------------------------------------" << nl;

        if (Pstream::master())
        {
            outputFilePtr() << runTime.value() << tab
                            << runTime.elapsedCpuTime() << tab
                            << resRho << tab
                            << resRhoU << tab
                            << resRhoE << tab
                            << resNuTilda << endl;
        }
        if (resRho < tolerance) runTime.writeAndEnd();
        else runTime.write();
	
        Info << "# ExecutionTime          [s] = " << runTime.elapsedCpuTime()  << nl
             << "# ClockTime              [s] = " << runTime.elapsedClockTime() << nl
             << "----------------------------------------" << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //

