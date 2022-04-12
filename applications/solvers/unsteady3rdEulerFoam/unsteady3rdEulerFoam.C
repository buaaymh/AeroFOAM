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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "createFields.H"
    // #include "createControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const scalar k11_22 = (3.0 - Sqrt(3.0)) / 6.0;
    const scalar k21    = 1.0 / Sqrt(3.0);

    while (runTime.run())
    {
        runTime++;
        Info<< "Time = " << runTime.value() << " s" << nl;

        scalar dt = runTime.deltaT();

        volScalarField rho_0(solver.rho());
        volVectorField rhoU_0(solver.rhoU());
        volScalarField rhoE_0(solver.rhoE());

        volScalarField resRho_1(resRho);
        volVectorField resRhoU_1(resRhoU);
        volScalarField resRhoE_1(resRhoE);

        scalar L1_deltaRho_0 = 0.0;
        scalar L1_deltaRho   = 0.0;
        
        // Stage 1
        for (int i = 0; i < innerIter; ++i)
        {
            solver.evaluateFlowRes(resRho, resRhoU, resRhoE);
            volScalarField pseudoResRho((rho_0 - solver.rho())/dt + k11_22* resRho);
            volVectorField pseudoResRhoU((rhoU_0 - solver.rhoU())/dt + k11_22*resRhoU);
            volScalarField pseudoResRhoE((rhoE_0 - solver.rhoE())/dt + k11_22*resRhoE);
            if (i == 0)
            {  
                solver.solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE, L1_deltaRho_0);
            }
            else
            {
                solver.solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE, L1_deltaRho);
                if (L1_deltaRho/L1_deltaRho_0 < inner_eps) break;
            }
        }
        // Stage 2
        for (int i = 0; i < innerIter; ++i)
        {
            solver.evaluateFlowRes(resRho_1, resRhoU_1, resRhoE_1);
            volScalarField pseudoResRho((rho_0 - solver.rho())/dt + k11_22*(resRho+resRho_1));
            volVectorField pseudoResRhoU((rhoU_0 - solver.rhoU())/dt + k11_22*(resRhoU+resRhoU_1));
            volScalarField pseudoResRhoE((rhoE_0 - solver.rhoE())/dt + k11_22*(resRhoE+resRhoE_1));
            if (i == 0)
            {  
                solver.solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE, L1_deltaRho_0);
            }
            else
            {
                solver.solveFlowPseudoTimeSystem(dt, k11_22, pseudoResRho, pseudoResRhoU, pseudoResRhoE, L1_deltaRho);
                if (L1_deltaRho/L1_deltaRho_0 < inner_eps) break;
            }
        }

        solver.rho()  = rho_0  + 0.5 * dt * (resRho + resRho_1);
        solver.rhoU() = rhoU_0 + 0.5 * dt * (resRhoU + resRhoU_1);
        solver.rhoE() = rhoE_0 + 0.5 * dt * (resRhoE + resRhoE_1);

        solver.correctField();

        runTime.write();
	
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //

