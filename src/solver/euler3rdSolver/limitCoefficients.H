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
        forAll(mesh_.C(), cellI)
        {
            vector limiter = vector(1, 1, 1);
            Col5X1 deltaW = coefs_[cellI].array().abs().matrix().transpose()*polyNorm_[cellI];
            // Limit rho
            if (rho_[cellI]-deltaW(0) < 1e-6)
            {
                limiter[0] = rho_[cellI]/deltaW(0);
                deltaW(0) *= limiter[0];
                coefs_[cellI].col(0) *= limiter[0];
            }
            // Limit rhoE
            if (rhoE_[cellI]-deltaW(4) < 1e-6)
            {
                limiter[2] = rhoE_[cellI]/deltaW(4);
                deltaW(4) *= limiter[2];
                coefs_[cellI].col(4) *= limiter[2];
            }
            // Limit rhoU
            scalar rhoU_magTop = sqrt(2.0*(rho_[cellI]-deltaW(0))*(rhoE_[cellI]-deltaW(4))+1e-6);
            scalar rhoU_magMax = mag(cmptMag(rhoU_[cellI]) + vector(deltaW(1),deltaW(2),deltaW(3)));
            scalar rhoU_magAvg = mag(rhoU_[cellI]);
            if (rhoU_magTop < rhoU_magAvg)
            {
                limiter[1] = 0;
                coefs_[cellI].block<9,3>(0,1) *= limiter[1];
            }
            else if (rhoU_magMax >= rhoU_magTop)
            {
                limiter[1] = min(1.0, mag(rhoU_magTop-rhoU_magAvg)/max(rhoU_magMax-rhoU_magAvg, SMALL));
                coefs_[cellI].block<9,3>(0,1) *= limiter[1];
            }
            if (cmptMin(limiter) < 1.0)  count++;
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}
