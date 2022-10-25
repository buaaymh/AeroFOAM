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
    Density-based compressible flow solver based on cell-centered 3rd schemes.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "solver.H"
#include "euler2ndSolver.H"
#include "eulerPrimVar3rdSolver.H"

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

        solver->evaluateFlowRes(resRho, resRhoU, resRhoE);
        const scalar L2Res = Foam::sqrt(gSumSqr(resRho));
        Info << "Residual for density = " << L2Res << endl;
        if (Pstream::master())
        {
            outputFilePtr() << runTime.value() << tab
                            << runTime.elapsedCpuTime() << tab
                            << L2Res << endl;
        }
        if (L2Res < tolerance) break;

        solver->solveFlowLinearSystem(resRho, resRhoU, resRhoE);
        solver->correctFields();
        runTime.write();
	
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //

