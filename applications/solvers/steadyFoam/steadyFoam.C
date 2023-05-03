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
    scalar L2ResRho = GREAT;

    while (runTime.run())
    {
        runTime++;
        Info << "========================================" << nl;
        Info << "# Iteration step    = " << runTime.value() << endl;
        Info << "========================================" << nl;
        solver->evaluateFlowRes(resRho, resRhoU, resRhoE);
        if (fluidProps.simulationType == "SATurb") solver->solveTurbulence();
        if (L2ResRho > 1)
        {
            Info << "# Impilict method   = LUSGS" << endl;
            solver->solveFlowLinearSystemByLUSGS(resRho, resRhoU, resRhoE, L2ResRho);
            Info << "# L2 rho residual   = " << L2ResRho << endl;
        }
        else
        {
            Info << "# Impilict method   = GMRES" << endl;
            solver->solveFlowLinearSystemByGMRES(resRho, resRhoU, resRhoE, L2ResRho);
            Info << "# L2 rho residual   = " << L2ResRho << endl;
        }
        solver->correctFields();
        
        if (Pstream::master())
        {
            outputFilePtr() << runTime.value() << tab
                            << runTime.elapsedCpuTime() << tab
                            << L2ResRho << endl;
        }
        if (L2ResRho < tolerance) runTime.writeAndEnd();
        else runTime.write();
	
        Info<< "# ExecutionTime     = " << runTime.elapsedCpuTime() << " s" << nl
            << "# ClockTime         = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //

