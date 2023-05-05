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
    const scalar nCells    = scalar(returnReduce(mesh.nCells(), sumOp<label>()));
    scalar initialResRho  = SMALL;
    scalar initialResRhoU = SMALL;
    scalar initialResRhoE = SMALL;

    while (runTime.run())
    {
        const word method = mesh.solution().subDict("SOLVER").lookupOrDefault<word>("method", "LUSGS");
        const scalar CFL  = mesh.solution().subDict("SOLVER").lookupOrDefault<scalar>("CFL", 1.0);
        runTime++;
        Info << "========================================" << nl;
        Info << "# " << method << " # Iteration Step = " << runTime.value() << endl;
        Info << "========================================" << nl;
        solver->evaluateFlowRes(resRho, resRhoU, resRhoE);
        Info << "# Local Courant          [-] = " << CFL << endl;
        Info << "----------------------------------------" << nl;
        scalar absoluteResRho  = Foam::sqrt(gSumSqr(resRho/mesh.V())/nCells);
        scalar absoluteResRhoU = Foam::sqrt(gSum(magSqr(resRhoU)/mesh.V())/nCells);
        scalar absoluteResRhoE = Foam::sqrt(gSumSqr(resRhoE/mesh.V())/nCells);

        if (runTime.value() <= 5)
        {
            initialResRho  = max(initialResRho,  absoluteResRho);
            initialResRhoU = max(initialResRhoU, absoluteResRhoU);
            initialResRhoE = max(initialResRhoE, absoluteResRhoE);
        }
        scalar relativeResRho  = absoluteResRho/initialResRho;
        scalar relativeResRhoU = absoluteResRhoU/initialResRhoU;
        scalar relativeResRhoE = absoluteResRhoE/initialResRhoE;
        Info << "# Continuity relativeRes [-] = " << relativeResRho  << endl;
        Info << "# Momentum   relativeRes [-] = " << relativeResRhoU << endl;
        Info << "# Energy     relativeRes [-] = " << relativeResRhoE << endl;
        if (fluidProps.simulationType == "SATurb") solver->solveTurbulence();
        Info << "----------------------------------------" << nl;

        if (method == "LUSGS") solver->solveFlowLinearSystemByLUSGS(resRho, resRhoU, resRhoE);
        else solver->solveFlowLinearSystemByGMRES(resRho, resRhoU, resRhoE);
        solver->correctFields();
        
        if (Pstream::master())
        {
            outputFilePtr() << runTime.value() << tab
                            << runTime.elapsedCpuTime() << tab
                            << relativeResRho << tab
                            << relativeResRhoU << tab
                            << relativeResRhoE << endl;
        }
        if (relativeResRho < tolerance) runTime.writeAndEnd();
        else runTime.write();
	
        Info << "# ExecutionTime          [s] = " << runTime.elapsedCpuTime()  << nl
             << "# ClockTime              [s] = " << runTime.elapsedClockTime() << nl
             << "----------------------------------------" << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //

