void Foam::euler3rdSolver::limitCoefficients()
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
        rho_mid[cellI] = rho_[cellI] + coefs_[cellI].col(0).dot(polynomialDn0(vector::zero, rDeltaXYZ_[cellI], basisMean_[cellI]));
    }
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rho_max[i] = max(rho_max[i], rho_[j]);
        rho_max[j] = max(rho_max[j], rho_[i]);
        error[j] += mag(rho_[i] + coefs_[i].col(0).dot(polynomialDn0( delta[faceI], rDeltaXYZ_[i], basisMean_[i])) - rho_mid[j]);
        error[i] += mag(rho_[j] + coefs_[j].col(0).dot(polynomialDn0(-delta[faceI], rDeltaXYZ_[j], basisMean_[j])) - rho_mid[i]);
    }
    label count = 0;
    forAll(mesh_.cells(), cellI)
    {
        IS_[cellI] = error[cellI]/(N_h_[cellI]*rho_max[cellI]);
        if (IS_[cellI] > IS)
        {
            isTrouble[cellI] = true;
            count++;
        }
    }
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "Trouble cell [%] = " << percent;

    /* Limiting d2 coefficients of polynomials */
    const Switch characteristic = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("characteristic", false);
    const scalar np = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("np", 10);
    std::vector<Mat6X5> d2CoefsLimited(count);
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<vector> normalList;
        std::vector<Mat5X6> d2CoefsList;
        normalList.reserve(mesh_.cellCells()[cellI].size());
        d2CoefsList.reserve(mesh_.cellCells()[cellI].size()+1);
        d2CoefsList.emplace_back(coefs_[cellI].block<6,5>(3,0).transpose());
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                normalList.emplace_back(normal_[faceI]);
                d2CoefsList.emplace_back(transformD2Coefs(coefs_[j].block<6,5>(3,0), rDeltaXYZ_[cellI], rDeltaXYZ_[j]).transpose());
            }
        }
        if (characteristic)
        {
            d2CoefsLimited[count++] = Eigen_WBAP_L2_Limiter<6>(normalList, d2CoefsList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], fluidProps_.gamma, np).transpose();
        }
        else
        {
            d2CoefsLimited[count++] = WBAP_L2_Limiter<6>(d2CoefsList, np).transpose();
        }
    }
    /* Update d2 coefficients of polynomials */
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (isTrouble[cellI])
        {
            coefs_[cellI].block<6,5>(3,0) = d2CoefsLimited[count++];
        }
    }
    /* Limiting d1 coefficients of polynomials */
    std::vector<Mat3X5> d1CoefsLimited(count);
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<vector> normalList;
        std::vector<Mat5X3> d1CoefsList;
        normalList.reserve(mesh_.cellCells()[cellI].size());
        d1CoefsList.reserve(mesh_.cellCells()[cellI].size()+1);
        d1CoefsList.emplace_back(coefs_[cellI].block<3,5>(0,0).transpose());
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                normalList.emplace_back(normal_[faceI]);
                label j = owner[faceI];
                vector d_ij = delta[faceI];
                if (j == cellI)
                {
                    j = neighbour[faceI];
                    d_ij *= -1;
                }
                d1CoefsList.emplace_back(transformD1Coefs(coefs_[j], d_ij, rDeltaXYZ_[cellI], rDeltaXYZ_[j]).transpose());
            }
        }
        if (characteristic)
        {
            d1CoefsLimited[count++] = Eigen_WBAP_L2_Limiter<3>(normalList, d1CoefsList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], fluidProps_.gamma, np).transpose();
        }
        else
        {
            d1CoefsLimited[count++] = WBAP_L2_Limiter<3>(d1CoefsList, np).transpose();
        }
    }
    /* Update d1 coefficients of polynomials */
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (isTrouble[cellI])
        {
            coefs_[cellI].block<3,5>(0,0) = d1CoefsLimited[count++];
        }
    }

    /* Preserve positive coefficients of polynomials */
    const Switch positive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    if (positive)
    {
        count = 0;
        forAll(mesh_.C(), cellI)
        {
            const Col9X1 polyTemp
            (
                1, 1, 1,
                0.5*(1+basisMean_[cellI].xx()), 1+basisMean_[cellI].xy(),  1+basisMean_[cellI].xz(),
                0.5*(1+basisMean_[cellI].yy()), 1+basisMean_[cellI].yz(),
                0.5*(1+basisMean_[cellI].zz())
            );
            const Col5X1 deltaVars = coefs_[cellI].transpose().array().abs().matrix()*polyTemp;
            const scalar rho_min  = rho_[cellI] - deltaVars(0);
            const vector rhoU_max = cmptMag(rhoU_[cellI]) + vector(deltaVars(1), deltaVars(2), deltaVars(3));
            const scalar rhoE_min = rhoE_[cellI] - deltaVars(4);
            const scalar p_min = rhoE_min - 0.5*magSqr(rhoU_max)/rho_min;
            if (min(min(rho_min, rhoE_min), p_min) < 0)
            {
                count++;
                coefs_[cellI] = Mat9X5::Zero();
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}
