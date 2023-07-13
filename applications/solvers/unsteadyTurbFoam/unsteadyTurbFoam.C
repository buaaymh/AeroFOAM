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

    const scalar k11_22 = (3.0-Foam::sqrt(3.0))/6.0;
    const scalar k21    = 1.0/Foam::sqrt(3.0);
    const scalar nCells = scalar(returnReduce(mesh.nCells(), sumOp<label>()));
    const label innerIter = mesh.solution().subDict("SOLVER").lookupOrDefault<label>("innerIter", 20);
    const scalar tolerance = mesh.solution().subDict("SOLVER").lookupOrDefault<scalar>("tolerance", 1e-6);
    const scalar relTol = mesh.solution().subDict("SOLVER").lookupOrDefault<scalar>("relTol", 0.001);
    const scalar dt = runTime.deltaT().value();
    const scalarField dt_dv(dt/mesh.V().field());
    
    while (runTime.run())
    {
        runTime++;
        Info<< "Time = " << runTime.value() << " s" << nl << endl;

        scalarField rho_0(solver->rho());
        vectorField rhoU_0(solver->rhoU());
        scalarField rhoE_0(solver->rhoE());
        scalarField nuTilda_0(solver->nuTilda());

        scalarField resRho_2(resRho_1);
        vectorField resRhoU_2(resRhoU_1);
        scalarField resRhoE_2(resRhoE_1);
        scalarField resNuTilda_2(resNuTilda_1);

        scalar L1_deltaRho_0 = 0.0;
        scalar L1_deltaRho   = 0.0;
        
        // Stage 1
        label count = 0;
        for (int i = 0; i < innerIter; ++i)
        {
            count++;
            solver->evaluateFlowRes(resRho_1, resRhoU_1, resRhoE_1);
            solver->evaluateTurbRes(resNuTilda_1);
            if (fluidProps.withSourceTerm)
            {
                actuationSource->addSourceTerms(runTime.value(), resRho_1, resRhoU_1, resRhoE_1);
            }
            scalarField pseudoResRho((rho_0   - solver->rho()) /dt_dv + k11_22*resRho_1);
            vectorField pseudoResRhoU((rhoU_0 - solver->rhoU())/dt_dv + k11_22*resRhoU_1);
            scalarField pseudoResRhoE((rhoE_0 - solver->rhoE())/dt_dv + k11_22*resRhoE_1);
            scalarField pseudoResNuTilda((nuTilda_0 - solver->nuTilda())/dt_dv + k11_22* resNuTilda_1);
            if (i == 0)
            {
                solver->solveTurbPseudoTimeSystem(dt, k11_22, pseudoResNuTilda);
                solver->solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE);
                solver->correctFields();
                L1_deltaRho_0 = gSum(mag(solver->dRho()))/nCells;
                Info << "# L1(dRho_0) = " << setprecision(4) << L1_deltaRho_0 << endl;
                if (L1_deltaRho_0 < tolerance)  { L1_deltaRho = L1_deltaRho_0; break; }
            }
            else
            {
                solver->solveTurbPseudoTimeSystem(dt, k11_22, pseudoResNuTilda);
                solver->solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE);
                solver->correctFields();
                L1_deltaRho = gSum(mag(solver->dRho()))/nCells;
                Info << "# L1(dRho) [-] = "            << setprecision(4) << L1_deltaRho
                     << ", L1(dRho)/L1(dRho_0) [-] = " << setprecision(4) << L1_deltaRho/L1_deltaRho_0 << endl;
                if (L1_deltaRho/L1_deltaRho_0 < relTol || L1_deltaRho < tolerance)  break;
            }
        }
        Info << "----------------------------------------" << nl;

        // Stage 2
        count = 0;
        for (int i = 0; i < innerIter; ++i)
        {
            count++;
            solver->evaluateFlowRes(resRho_2, resRhoU_2, resRhoE_2);
            solver->evaluateTurbRes(resNuTilda_2);
            if (fluidProps.withSourceTerm)
            {
                actuationSource->addSourceTerms(runTime.value(), resRho_2, resRhoU_2, resRhoE_2);
            }
            scalarField pseudoResRho((rho_0   - solver->rho())/dt_dv  + k11_22*resRho_2  + k21*resRho_1);
            vectorField pseudoResRhoU((rhoU_0 - solver->rhoU())/dt_dv + k11_22*resRhoU_2 + k21*resRhoU_1);
            scalarField pseudoResRhoE((rhoE_0 - solver->rhoE())/dt_dv + k11_22*resRhoE_2 + k21*resRhoE_1);
            scalarField pseudoResNuTilda((nuTilda_0 - solver->nuTilda())/dt_dv + k11_22*resNuTilda_2 + k21*resNuTilda_1);
            if (i == 0)
            {
                solver->solveTurbPseudoTimeSystem(dt, k11_22, pseudoResNuTilda);
                solver->solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE);
                solver->correctFields();
                L1_deltaRho_0 = gSum(mag(solver->dRho()))/nCells;
                Info << "# L1(dRho_0) = " << setprecision(4) << L1_deltaRho_0 << endl;
                if (L1_deltaRho_0 < tolerance)  { L1_deltaRho = L1_deltaRho_0; break; }
            }
            else
            {
                solver->solveTurbPseudoTimeSystem(dt, k11_22, pseudoResNuTilda);
                solver->solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE);
                solver->correctFields();
                L1_deltaRho = gSum(mag(solver->dRho()))/nCells;
                Info << "# L1(dRho) [-] = "            << setprecision(4) << L1_deltaRho
                     << ", L1(dRho)/L1(dRho_0) [-] = " << setprecision(4) << L1_deltaRho/L1_deltaRho_0 << endl;
                if (L1_deltaRho/L1_deltaRho_0 < relTol || L1_deltaRho < tolerance)  break;
            }
        }
        Info << "----------------------------------------" << nl;

        solver->rho()  = rho_0  + 0.5 * dt_dv * (resRho_1 + resRho_2);
        solver->rhoU() = rhoU_0 + 0.5 * dt_dv * (resRhoU_1 + resRhoU_2);
        solver->rhoE() = rhoE_0 + 0.5 * dt_dv * (resRhoE_1 + resRhoE_2);
        solver->nuTilda() = nuTilda_0 + 0.5 * dt_dv * (resNuTilda_1+ resNuTilda_2);
        
        solver->correctFields();
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

