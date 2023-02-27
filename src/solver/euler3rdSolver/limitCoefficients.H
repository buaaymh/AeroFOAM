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
            scalar rho_min = rho_[cellI] - deltaW(0);
            if (rho_min < 1e-6)
            {
                rho_min = max(rho_min, 1e-6);
                limiter[0] = (rho_[cellI]-rho_min)/max(deltaW(0), SMALL);
                coefs_[cellI].col(0) *= limiter[0];
            }
            // Limit rhoE
            scalar rhoE_min = rhoE_[cellI]-deltaW(4);
            if (rhoE_min < 1e-6)
            {
                rhoE_min = max(rhoE_min, 1e-6);
                limiter[2] = (rhoE_[cellI]-rhoE_min)/max(deltaW(4), SMALL);
                coefs_[cellI].col(4) *= limiter[2];
            }
            // Limit rhoU, rhoE
            scalar rhoU_top = sqrt(2.0*rho_min*rhoE_min);
            scalar rhoU_max = mag(cmptMag(rhoU_[cellI]) + vector(deltaW(1),deltaW(2),deltaW(3)));
            scalar rhoU_avg = mag(rhoU_[cellI]);
            if (rhoU_top < rhoU_avg)
            {
                limiter[1] = 0;
                limiter[2] = 0;
                coefs_[cellI].block<9,3>(0,1) *= limiter[1];
                coefs_[cellI].col(4) *= limiter[2];
            }
            else if (rhoU_max >= rhoU_top)
            {
                limiter[1] = min(1.0, mag(rhoU_top-rhoU_avg)/max(rhoU_max-rhoU_avg, SMALL));
                coefs_[cellI].block<9,3>(0,1) *= limiter[1];
            }
            if (cmptMin(limiter) < 1.0)
            {
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}
