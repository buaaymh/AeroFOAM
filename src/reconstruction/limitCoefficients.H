void Foam::Reconstruction::limitCoefficients()
{
    /* Detect trouble cells */
    surfaceVectorField normal(mesh_.Sf()/mesh_.magSf());
    vectorField delta_(mesh_.delta());
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    parameter_.nTroubled = 0;
    #include "troubleIndicator.H"
    /* Limiting d2 coefficients of polynomials */
    if (parameter_.limiter == "Slope")
    {
        #include "slopeLimiter.H"
    }
    else if (parameter_.limiter == "WBAP")
    {
        #include "WBAPLimiter.H"
    }
    else if (parameter_.limiter == "New")
    {
        #include "NewWBAPLimiter.H"
    }
    else
    {
        Pout << "# Error in limiter type #" << endl;
    }

    /* Preserve positive coefficients of polynomials */
    if (parameter_.positive)
    {
        parameter_.nNegative = 0;
        forAll(mesh_.C(), cellI)
        {
            scalar factor = 1.0;
            Col5X1 delta = coefs_[cellI].array().abs().matrix().transpose()*polyNorm_[cellI];
            if (0.99*rho_[cellI]  < delta(0)) factor = min(factor, 0.99*rho_[cellI]/delta(0));
            if (0.99*rhoE_[cellI] < delta(4)) factor = min(factor, 0.99*rhoE_[cellI]/delta(4));
            if (factor < 1.0)
            {
                coefs_[cellI] *= factor;
                parameter_.nNegative++;
                IS_[cellI] = 100;
            }
        }
    }
    Info << "# Troubled [%] = " << setprecision(4) << 100*scalar(returnReduce(parameter_.nTroubled, sumOp<label>()))/parameter_.nCells
         << ", Negative [%] = " << setprecision(4) << 100*scalar(returnReduce(parameter_.nNegative, sumOp<label>()))/parameter_.nCells << endl;
}
