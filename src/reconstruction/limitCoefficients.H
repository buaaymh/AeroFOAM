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
                coefs_[cellI] *= 0;
                parameter_.nNegative++;
                IS_[cellI] = 100;
            }
        }
    }
    Info << "# Troubled [%] = " << setprecision(4) << 100*scalar(returnReduce(parameter_.nTroubled, sumOp<label>()))/parameter_.nCells
         << ", Negative [%] = " << setprecision(4) << 100*scalar(returnReduce(parameter_.nNegative, sumOp<label>()))/parameter_.nCells << endl;
}
