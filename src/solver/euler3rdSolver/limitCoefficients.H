void Foam::euler3rdSolver::limitCoefficients()
{
    /* Detect trouble cells */
    const scalar IS = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("IS", 1.0);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    std::vector<bool> isTrouble(mesh_.nCells(), false); 
    label count = 0;
    #include "troubleIndicator.H"
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "Trouble cell [%] = " << percent;

    /* Limiting d2 coefficients of polynomials */
    const Switch characteristic = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("characteristic", false);
    const word limiter = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<word>("limiter", "Slope1");
    if (limiter == "Slope1")
    {
        #include "slope1Limiter.H"
    }
    else if (limiter == "WBAP")
    {
        #include "WBAPLimiter.H"
    }
    else
    {
        Pout << "# Error in limiter type #" << endl;
    }

    /* Preserve positive coefficients of polynomials */
    const Switch positive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    if (positive)
    {
        count = 0;
        volVectorField pGrad = fvc::grad(p_);
        forAll(mesh_.C(), cellI)
        {
            const scalar p_min   = 0.5*p_[cellI] - cmptSum(cmptMag(cmptDivide(pGrad[cellI], rDeltaXYZ_[cellI])));
            // if (p_min < 0.0)
            // {
            //     count++;
            //     coefs_[cellI] = Mat9X5::Zero();
            // }
            if (p_min < 0.0)
            {
                count++;
                Col5X1 deltaW = coefs_[cellI].array().abs().matrix().transpose()*polyNorm_[cellI];
                // Limit rho
                scalar limiter = min(1.0, rho_[cellI]/max(deltaW(0), SMALL));
                coefs_[cellI].col(0) *= limiter;
                deltaW(0) *= limiter;
                // Limit rhoE
                limiter = min(1.0, rhoE_[cellI]/max(deltaW(4), SMALL));
                coefs_[cellI].col(4) *= limiter;
                deltaW(4) *= limiter;
                // Limit rhoU
                scalar rhoU = sqrt(2.0*(rhoE_[cellI]+deltaW(4))*(rho_[cellI]+deltaW(0)));
                scalar rhoU_max = mag(cmptMag(rhoU_[cellI]) + vector(deltaW(1),deltaW(2),deltaW(3)));
                scalar temp = mag(rhoU_[cellI]);
                limiter = min(1.0, mag(rhoU-temp)/max(rhoU_max-temp, SMALL));
                coefs_[cellI].block<9,3>(0,1) *= limiter;
            }
        }
        forAll(mesh_.boundary(), patchI)
        { 
            const word type = mesh_.boundary()[patchI].type();
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            if (type == "symmetryPlane" || type == "wall" || type == "symmetry")
            {
                forAll(bfaceCells, faceI)
                {
                    const label cellI = bfaceCells[faceI];
                    Col5X1 deltaW = coefs_[cellI].array().abs().matrix().transpose()*polyNorm_[cellI];
                    // Limit rhoU
                    scalar rhoU = sqrt(2.0*(rhoE_[cellI]+deltaW(4))*(rho_[cellI]+deltaW(0)));
                    scalar rhoU_max = mag(cmptMag(rhoU_[cellI]) + vector(deltaW(1),deltaW(2),deltaW(3)));
                    scalar temp = mag(rhoU_[cellI]);
                    scalar limiter = min(1.0, mag(rhoU-temp)/max(rhoU_max-temp, SMALL));
                    if (limiter < 1.0)
                    {
                        count++;
                        coefs_[cellI].block<9,3>(0,1) *= limiter;
                    }
                }
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}
