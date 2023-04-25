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
#include "euler3rdSolver.H"
#include "turb2ndSolver.H"

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
    scalar L2Res_0 = 0;
    
    while (runTime.run())
    {
        runTime++;
        Info << "Time = " << runTime.value() << " s" << nl << endl;;
        solver->evaluateFlowRes(resRho, resRhoU, resRhoE);
        const scalar L2ResRho  = Foam::sqrt(gSumSqr(resRho));
        const scalar L2ResRhoU = Foam::sqrt(gSum(magSqr(resRhoU)));
        const scalar L2ResRhoE = Foam::sqrt(gSumSqr(resRhoE));
        const scalar L2Res     = Foam::sqrt(sqr(L2ResRho) + sqr(L2ResRhoU) + sqr(L2ResRhoE));
        if (L2Res_0 < SMALL) { L2Res_0 = L2Res; }
        const scalar relRes = L2Res/L2Res_0;
        
        Info << "# FGMRES solving for rho  = " << L2ResRho << endl;
        Info << "# FGMRES solving for rhoU = " << L2ResRhoU << endl;
        Info << "# FGMRES solving for rhoE = " << L2ResRhoE << endl;
        Info << "# Relative residual = " << relRes << endl;

        if (relRes > 1e-1) solver->updateCFL(20);
        else if ((relRes < 1e-1) && (relRes > 1e-2)) solver->updateCFL(100);
        else if ((relRes < 1e-2) && (relRes > 1e-4)) solver->updateCFL(500);
        else solver->updateCFL(5000);

        if (Pstream::master())
        {
            outputFilePtr() << runTime.value() << tab
                            << runTime.elapsedCpuTime() << tab
                            << relRes << endl;
        }
        if (fluidProps.simulationType == "SATurb") solver->solveTurbulence();
        solver->solveFlowLinearSystem(resRho, resRhoU, resRhoE);
        solver->correctFields();

        if (relRes < tolerance) runTime.writeAndEnd();
        else runTime.write();
	
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //

