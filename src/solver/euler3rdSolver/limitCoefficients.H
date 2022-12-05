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
    else if (limiter == "Slope2")
    {
        #include "slope2Limiter.H"
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
            if (p_min < 0.0)
            {
                count++;
                coefs_[cellI] = Mat9X5::Zero();
            }
            else
            {
                const Col9X1 polyTemp
                (
                    1, 1, 1,
                    0.5*(1+basisMean_[cellI].xx()), 1+basisMean_[cellI].xy(),  1+basisMean_[cellI].xz(),
                    0.5*(1+basisMean_[cellI].yy()), 1+basisMean_[cellI].yz(),
                    0.5*(1+basisMean_[cellI].zz())
                );
                const scalar deltaRho  = coefs_[cellI].col(0).dot(polyTemp);
                const scalar deltaRhoE = coefs_[cellI].col(4).dot(polyTemp);
                if (min(rho_[cellI]-deltaRho, rhoE_[cellI]-deltaRhoE) < 0)
                {
                    count++;
                    const scalar temp = min(rho_[cellI]/deltaRho, rhoE_[cellI]/deltaRhoE);
                    coefs_[cellI] *= temp;
                }
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}
