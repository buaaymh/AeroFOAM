void Foam::eulerPrimVar3rdSolver::limitCoefficients()
{
    /* Detect trouble cells */
    const scalar IS = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("IS", 1.0);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField delta = mesh_.delta();
    scalarField rhoE_max(rhoE_.primitiveField());
    scalarField rhoE_mid(mesh_.nCells());
    scalarField error(mesh_.nCells(), 0.0);
    std::vector<bool> isTrouble(mesh_.nCells(), false);
    forAll(mesh_.cells(), cellI)
    {
        scalar rho; vector U; scalar p;
        evaluateVars(rho, U, p, rho_[cellI], U_[cellI], p_[cellI], rhoGrad_[cellI], UGrad_[cellI], pGrad_[cellI],
                     d2Var_[cellI], basisMean_[cellI], rDeltaXYZ_[cellI], vector::zero);
        rhoE_mid[cellI] = p/(fluidProps_.gamma-1.0) + 0.5*rho*magSqr(U);
    }
    forAll(owner, faceI)
    {
        scalar rho; vector U; scalar p;
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rhoE_max[i] = max(rhoE_max[i], rhoE_[j]);
        rhoE_max[j] = max(rhoE_max[j], rhoE_[i]);
        evaluateVars(rho, U, p, rho_[i], U_[i], p_[i], rhoGrad_[i], UGrad_[i], pGrad_[i],
                     d2Var_[i], basisMean_[i], rDeltaXYZ_[i], delta[faceI]);
        error[j] += mag(p/(fluidProps_.gamma-1.0)+0.5*rho*magSqr(U) - rhoE_mid[j]);
        evaluateVars(rho, U, p, rho_[j], U_[j], p_[j], rhoGrad_[j], UGrad_[j], pGrad_[j],
                     d2Var_[j], basisMean_[j], rDeltaXYZ_[j], -delta[faceI]);
        error[i] += mag(p/(fluidProps_.gamma-1.0)+0.5*rho*magSqr(U) - rhoE_mid[i]);
    }
    label count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (error[cellI]/(N_h_[cellI]*rhoE_max[cellI]) > IS)
        {
            isTrouble[cellI] = true;
            count++;
        }
    }
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << ", Trouble cell [%] = " << percent;

    /* Limiting d2 coefficients of polynomials */
    const scalar np = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("np", 10);
    std::vector<Mat5X6> d2VarLimited(count);
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<Mat5X6> d2VarList;
        d2VarList.reserve(mesh_.cellCells()[cellI].size()+1);
        d2VarList.emplace_back(d2Var_[cellI]);
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                d2VarList.emplace_back(d2Var_[j]);
            }
        }
        d2VarLimited[count++] = WBAP_L2_Limiter<6>(d2VarList, np);
    }

    /* Update d2 coefficients of polynomials */
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (isTrouble[cellI])
        {
            d2Var_[cellI] = d2VarLimited[count];
            count++;
        }
    }

    /* Update d1Var */
    forAll(mesh_.C(), cellI)
    {
        d1Var_[cellI] << rhoGrad_[cellI][0], rhoGrad_[cellI][1], rhoGrad_[cellI][2],
                         UGrad_[cellI][0],   UGrad_[cellI][1],   UGrad_[cellI][2],
                         UGrad_[cellI][3],   UGrad_[cellI][4],   UGrad_[cellI][5],
                         UGrad_[cellI][6],   UGrad_[cellI][7],   UGrad_[cellI][8],
                         pGrad_[cellI][0],   pGrad_[cellI][1],   pGrad_[cellI][2];
    }

    /* Limiting d1 coefficients of polynomials */
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<Mat5X3> d1VarList;
        d1VarList.reserve(mesh_.cellCells()[cellI].size()+1);
        d1VarList.emplace_back(d1Var_[cellI]);
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                const vector delta_ij = mesh_.C()[cellI] - mesh_.C()[j];
                d1VarList.emplace_back(evaluateGrad(delta_ij, d1Var_[j], d2Var_[j], rDeltaXYZ_[j]));
            }
        }
        Mat5X3 d1VarLimited = WBAP_L2_Limiter<3>(d1VarList, np);
        rhoGrad_[cellI] = vector(d1VarLimited(0,0), d1VarLimited(0,1), d1VarLimited(0,2));
        UGrad_[cellI]   = tensor(d1VarLimited(1,0), d1VarLimited(1,1), d1VarLimited(1,2),
                                 d1VarLimited(2,0), d1VarLimited(2,1), d1VarLimited(2,2),
                                 d1VarLimited(3,0), d1VarLimited(3,1), d1VarLimited(3,2));
        pGrad_[cellI]   = vector(d1VarLimited(4,0), d1VarLimited(4,1), d1VarLimited(4,2));
    }

    /* Update d2 coefficients of polynomials */
    forAll(mesh_.C(), cellI)
    {
        d2Rho_[cellI] = symmTensor(d2Var_[cellI](0,0), d2Var_[cellI](0,1), d2Var_[cellI](0,2),
                                   d2Var_[cellI](0,3), d2Var_[cellI](0,4), d2Var_[cellI](0,5));
        d2Ux_[cellI]  = symmTensor(d2Var_[cellI](1,0), d2Var_[cellI](1,1), d2Var_[cellI](1,2),
                                   d2Var_[cellI](1,3), d2Var_[cellI](1,4), d2Var_[cellI](1,5));
        d2Uy_[cellI]  = symmTensor(d2Var_[cellI](2,0), d2Var_[cellI](2,1), d2Var_[cellI](2,2),
                                   d2Var_[cellI](2,3), d2Var_[cellI](2,4), d2Var_[cellI](2,5));
        d2Uz_[cellI]  = symmTensor(d2Var_[cellI](3,0), d2Var_[cellI](3,1), d2Var_[cellI](3,2),
                                   d2Var_[cellI](3,3), d2Var_[cellI](3,4), d2Var_[cellI](3,5));
        d2P_[cellI]   = symmTensor(d2Var_[cellI](4,0), d2Var_[cellI](4,1), d2Var_[cellI](4,2),
                                   d2Var_[cellI](4,3), d2Var_[cellI](4,4), d2Var_[cellI](4,5));
    }

    /* Preserve positive coefficients of polynomials */
    const Switch positive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    if (positive)
    {
        count = 0;
        forAll(mesh_.C(), cellI)
        {
            if (isP0Cell_[cellI]) continue;
            const symmTensor basisFuncs
            (
                0.5*(1+basisMean_[cellI].xx()), 1+basisMean_[cellI].xy(),  1+basisMean_[cellI].xz(),
                                           0.5*(1+basisMean_[cellI].yy()), 1+basisMean_[cellI].yz(),
                                                                      0.5*(1+basisMean_[cellI].zz())
            );
            const scalar deltaRho = cmptSum(cmptMag(cmptDivide(rhoGrad_[cellI], rDeltaXYZ_[cellI])))
                                  + cmptSum(cmptMag(cmptMultiply(d2Rho_[cellI], basisFuncs)));
            const scalar deltaP   = cmptSum(cmptMag(cmptDivide(pGrad_[cellI], rDeltaXYZ_[cellI])))
                                  + cmptSum(cmptMag(cmptMultiply(d2P_[cellI], basisFuncs)));
            scalar temp = min(rho_[cellI]-deltaRho, p_[cellI]-deltaP);
            if (temp < 0)
            {
                temp = min(rho_[cellI]/deltaRho, p_[cellI]/deltaP);
                rhoGrad_[cellI] *= temp;
                UGrad_[cellI] *= temp;
                pGrad_[cellI] *= temp;
                d2Rho_[cellI] *= temp;
                d2Ux_[cellI] *= temp;
                d2Uy_[cellI] *= temp;
                d2Uz_[cellI] *= temp;
                d2P_[cellI] *= temp;
                d2Var_[cellI] *= temp;
                count++;
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}