void Foam::euler3rdSolver::limitCoefficients()
{
    /* Detect trouble cells */
    const scalar IS = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("IS", 1.0);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    label count = 0;
    #include "troubleIndicator.H"
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "Trouble cell [%] = " << percent;

    /* Limiting d2 coefficients of polynomials */
    const Switch characteristic = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("characteristic", false);
    const word limiter = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<word>("limiter", "Slope1");
    if (limiter == "Slope")
    {
        #include "slopeLimiter.H"
    }
    else if (limiter == "WBAP")
    {
        #include "WBAPLimiter.H"
    }
    else if (limiter == "New")
    {
        #include "NewWBAPLimiter.H"
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
            if (!isTrouble_[cellI]) continue;
            scalar factor = 1.0;
            scalar delta = coefs_[cellI].col(0).array().abs().matrix().dot(polyNorm_[cellI]);
            if (rho_[cellI]  < delta) factor = min(factor, rho_[cellI]/delta);
            delta = coefs_[cellI].col(4).array().abs().matrix().dot(polyNorm_[cellI]);
            if (rhoE_[cellI] < delta) factor = min(factor, rhoE_[cellI]/delta);
            if (factor < 1.0)
            {
                coefs_[cellI] *= factor;
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}
