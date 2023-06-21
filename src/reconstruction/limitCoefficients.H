void Foam::Reconstruction::limitCoefficients()
{
    /* Detect trouble cells */
    const scalar IS = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<scalar>("IS", 1.0);
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());
    vectorField delta_(mesh_.delta());
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    label count = 0;
    #include "troubleIndicator.H"
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "# Troubled cell          [%] = " << percent << endl;

    /* Limiting d2 coefficients of polynomials */
    const Switch characteristic = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<Switch>("characteristic", false);
    const word limiter = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<word>("limiter", "Slope");
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
    const Switch positive = mesh_.schemes().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    if (positive)
    {
        count = 0;
        forAll(mesh_.C(), cellI)
        {
            scalar factor = 1.0;
            Col5X1 delta = coefs_[cellI].array().abs().matrix().transpose()*polyNorm_[cellI];
            if (rho_[cellI]  < 2*delta(0)) factor = min(factor, 0.5*rho_[cellI]/delta(0));
            if (rhoE_[cellI] < 2*delta(4)) factor = min(factor, 0.5*rhoE_[cellI]/delta(4));
            scalar dRho  = -delta(0);
            vector dRhoU = vector(delta(1),delta(2),delta(3));
            scalar dRhoE = -delta(4);
            scalar rhoNew  = rho_[cellI] + dRho*factor;
            vector rhoUNew = cmptMag(rhoU_[cellI]) + dRhoU*factor;
            scalar rhoENew = rhoE_[cellI]+ dRhoE*factor;
            scalar pNew = (Gamma-1.0)*(rhoENew - 0.5*magSqr(rhoUNew)/rhoNew);
            if (pNew <= SMALL)
                factor = min(factor, solveForPressure(rho_[cellI], cmptMag(rhoU_[cellI]), rhoE_[cellI],
                                                      dRho, dRhoU, dRhoE, SMALL, factor));
            if (factor < 1.0)
            {
                coefs_[cellI] *= factor;
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << "# Negative cell          [%] = " << percent << endl;
    }
}
