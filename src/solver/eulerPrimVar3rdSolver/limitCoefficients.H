void Foam::eulerPrimVar3rdSolver::limitCoefficients()
{
    /* Detect trouble cells */
    const scalar IS = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("IS", 1.0);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField delta = mesh_.delta();
    scalarField rho_max(rho_.primitiveField());
    scalarField rho_mid(mesh_.nCells());
    scalarField error(mesh_.nCells(), 0.0);
    std::vector<bool> isTrouble(mesh_.nCells(), false);
    forAll(mesh_.cells(), cellI)
    {
        d2Rho_[cellI] = symmTensor(d2Var_[cellI](0,0), d2Var_[cellI](0,1), d2Var_[cellI](0,2),
                                   d2Var_[cellI](0,3), d2Var_[cellI](0,4), d2Var_[cellI](0,5));
        rho_mid[cellI] = rho_[cellI]+scalarPoly(vector::zero, rhoGrad_[cellI], d2Rho_[cellI], rDeltaXYZ_[cellI], basisMean_[cellI]);
    }
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        error[j] += mag(rho_[i]+scalarPoly( delta[faceI], rhoGrad_[i], d2Rho_[i], rDeltaXYZ_[i], basisMean_[i])-rho_mid[j]);
        error[i] += mag(rho_[j]+scalarPoly(-delta[faceI], rhoGrad_[j], d2Rho_[j], rDeltaXYZ_[j], basisMean_[j])-rho_mid[i]);
    }
    label count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (error[cellI]/(N_h_[cellI]*rho_max[cellI]) > IS)
        {
            isTrouble[cellI] = true;
            count++;
        }
    }
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "Trouble cell [%] = " << percent << endl;

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
        d2VarLimited[count++] = WBAP_L2_Limiter<5, 6>(d2VarList, np);
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
                         TGrad_[cellI][0],   TGrad_[cellI][1],   TGrad_[cellI][2];
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
        Mat5X3 d1VarLimited = WBAP_L2_Limiter<5, 3>(d1VarList, np);
        rhoGrad_[cellI] = vector(d1VarLimited(0,0), d1VarLimited(0,1), d1VarLimited(0,2));
        UGrad_[cellI]   = tensor(d1VarLimited(1,0), d1VarLimited(1,1), d1VarLimited(1,2),
                                 d1VarLimited(2,0), d1VarLimited(2,1), d1VarLimited(2,2),
                                 d1VarLimited(3,0), d1VarLimited(3,1), d1VarLimited(3,2));
        TGrad_[cellI]   = vector(d1VarLimited(4,0), d1VarLimited(4,1), d1VarLimited(4,2));
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
        d2T_[cellI]   = symmTensor(d2Var_[cellI](4,0), d2Var_[cellI](4,1), d2Var_[cellI](4,2),
                                   d2Var_[cellI](4,3), d2Var_[cellI](4,4), d2Var_[cellI](4,5));
    }

    /* Preserve positive coefficients of polynomials */
    const Switch positive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    if (positive)
    {
        forAll(mesh_.C(), cellI)
        {
            if (isP0Cell_[cellI]) continue;
            const symmTensor basisFuncs
            (
                0.5*(1-basisMean_[cellI].xx()), 1-basisMean_[cellI].xy(),  1-basisMean_[cellI].xz(),
                                           0.5*(1-basisMean_[cellI].yy()), 1-basisMean_[cellI].yz(),
                                                                      0.5*(1-basisMean_[cellI].yy())
            );
            const scalar deltaRho = cmptSum(cmptDivide(cmptMag(rhoGrad_[cellI]), rDeltaXYZ_[cellI]))
                                  + cmptSum(cmptMultiply(cmptMag(d2Rho_[cellI]), basisFuncs));
            const scalar deltaT   = cmptSum(cmptDivide(cmptMag(TGrad_[cellI]), rDeltaXYZ_[cellI]))
                                  + cmptSum(cmptMultiply(cmptMag(d2T_[cellI]), basisFuncs));
            scalar temp = min(rho_[cellI]-deltaRho, T_[cellI]-deltaT);
            if (temp < 0)
            {
                temp = min(rho_[cellI]/deltaRho, T_[cellI]/deltaT);
                rhoGrad_[cellI] *= temp;
                UGrad_[cellI] *= temp;
                TGrad_[cellI] *= temp;
                d2Rho_[cellI] *= temp;
                d2Ux_[cellI] *= temp;
                d2Uy_[cellI] *= temp;
                d2Uz_[cellI] *= temp;
                d2T_[cellI] *= temp;
            }
        }
    }
}