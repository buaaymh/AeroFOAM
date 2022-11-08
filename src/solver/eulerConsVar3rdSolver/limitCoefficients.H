void Foam::eulerConsVar3rdSolver::limitCoefficients()
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
        d2RhoE_[cellI] = symmTensor(d2Var_[cellI](4,0), d2Var_[cellI](4,1), d2Var_[cellI](4,2),
                                    d2Var_[cellI](4,3), d2Var_[cellI](4,4), d2Var_[cellI](4,5));
        rhoE_mid[cellI] = rhoE_[cellI]+scalarPoly(vector::zero, rhoEGrad_[cellI], d2RhoE_[cellI], rDeltaXYZ_[cellI], basisMean_[cellI]);
    }
    forAll(owner, faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rhoE_max[i] = max(rhoE_max[i], rhoE_[j]);
        rhoE_max[j] = max(rhoE_max[j], rhoE_[i]);
        error[j] += mag(rhoE_[i]+scalarPoly( delta[faceI], rhoEGrad_[i], d2RhoE_[i], rDeltaXYZ_[i], basisMean_[i])-rhoE_mid[j]);
        error[i] += mag(rhoE_[j]+scalarPoly(-delta[faceI], rhoEGrad_[j], d2RhoE_[j], rDeltaXYZ_[j], basisMean_[j])-rhoE_mid[i]);
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
    const Switch characteristic = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("characteristic", false);
    const scalar np = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("np", 10);
    std::vector<Mat5X6> d2VarLimited(count);
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<vector> normalList;
        std::vector<Mat5X6> d2VarList;
        normalList.reserve(mesh_.cellCells()[cellI].size());
        d2VarList.reserve(mesh_.cellCells()[cellI].size()+1);
        d2VarList.emplace_back(d2Var_[cellI]);
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                normalList.emplace_back(normal_[faceI]);
                d2VarList.emplace_back(d2Var_[j]);
            }
        }
        if (characteristic)
        {
            d2VarLimited[count++] = Eigen_WBAP_L2_Limiter<6>(normalList, d2VarList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], fluidProps_.gamma, np);
        }
        else
        {
            d2VarLimited[count++] = WBAP_L2_Limiter<6>(d2VarList, np);
        }
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
        d1Var_[cellI] << rhoGrad_[cellI][0],  rhoGrad_[cellI][1],  rhoGrad_[cellI][2],
                         rhoUGrad_[cellI][0], rhoUGrad_[cellI][1], rhoUGrad_[cellI][2],
                         rhoUGrad_[cellI][3], rhoUGrad_[cellI][4], rhoUGrad_[cellI][5],
                         rhoUGrad_[cellI][6], rhoUGrad_[cellI][7], rhoUGrad_[cellI][8],
                         rhoEGrad_[cellI][0], rhoEGrad_[cellI][1], rhoEGrad_[cellI][2];
    }

    /* Limiting d1 coefficients of polynomials */
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<vector> normalList;
        std::vector<Mat5X3> d1VarList;
        normalList.reserve(mesh_.cellCells()[cellI].size());
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
                normalList.emplace_back(normal_[faceI]);
                d1VarList.emplace_back(evaluateGrad(delta_ij, d1Var_[j], d2Var_[j], rDeltaXYZ_[j]));
            }
        }
        Mat5X3 d1VarLimited;
        if (characteristic)
        {
            d1VarLimited = Eigen_WBAP_L2_Limiter<3>(normalList, d1VarList, rho_[cellI], rhoU_[cellI], rhoE_[cellI], fluidProps_.gamma, np);
        }
        else
        {
            d1VarLimited = WBAP_L2_Limiter<3>(d1VarList, np);
        }
        rhoGrad_[cellI]  = vector(d1VarLimited(0,0), d1VarLimited(0,1), d1VarLimited(0,2));
        rhoUGrad_[cellI] = tensor(d1VarLimited(1,0), d1VarLimited(1,1), d1VarLimited(1,2),
                                  d1VarLimited(2,0), d1VarLimited(2,1), d1VarLimited(2,2),
                                  d1VarLimited(3,0), d1VarLimited(3,1), d1VarLimited(3,2));
        rhoEGrad_[cellI] = vector(d1VarLimited(4,0), d1VarLimited(4,1), d1VarLimited(4,2));
    }

    /* Update d2 coefficients of polynomials */
    forAll(mesh_.C(), cellI)
    {
        d2Rho_[cellI]   = symmTensor(d2Var_[cellI](0,0), d2Var_[cellI](0,1), d2Var_[cellI](0,2),
                                     d2Var_[cellI](0,3), d2Var_[cellI](0,4), d2Var_[cellI](0,5));
        d2RhoUx_[cellI] = symmTensor(d2Var_[cellI](1,0), d2Var_[cellI](1,1), d2Var_[cellI](1,2),
                                     d2Var_[cellI](1,3), d2Var_[cellI](1,4), d2Var_[cellI](1,5));
        d2RhoUy_[cellI] = symmTensor(d2Var_[cellI](2,0), d2Var_[cellI](2,1), d2Var_[cellI](2,2),
                                     d2Var_[cellI](2,3), d2Var_[cellI](2,4), d2Var_[cellI](2,5));
        d2RhoUz_[cellI] = symmTensor(d2Var_[cellI](3,0), d2Var_[cellI](3,1), d2Var_[cellI](3,2),
                                     d2Var_[cellI](3,3), d2Var_[cellI](3,4), d2Var_[cellI](3,5));
        d2RhoE_[cellI]  = symmTensor(d2Var_[cellI](4,0), d2Var_[cellI](4,1), d2Var_[cellI](4,2),
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
            const scalar deltaRho  = cmptSum(cmptMag(cmptDivide(rhoGrad_[cellI], rDeltaXYZ_[cellI])))
                                   + cmptSum(cmptMag(cmptMultiply(d2Rho_[cellI], basisFuncs)));
            const scalar deltaRhoE = cmptSum(cmptMag(cmptDivide(rhoEGrad_[cellI], rDeltaXYZ_[cellI])))
                                   + cmptSum(cmptMag(cmptMultiply(d2RhoE_[cellI], basisFuncs)));
            scalar temp = min(rho_[cellI]-deltaRho, rhoE_[cellI]-deltaRhoE);
            if (temp > 0)
            {
                const scalar deltaRhoUx = cmptSum(cmptMag(cmptDivide(rhoUGrad_[cellI].x(), rDeltaXYZ_[cellI])))
                                        + cmptSum(cmptMag(cmptMultiply(d2RhoUx_[cellI], basisFuncs)));
                const scalar deltaRhoUy = cmptSum(cmptMag(cmptDivide(rhoUGrad_[cellI].y(), rDeltaXYZ_[cellI])))
                                        + cmptSum(cmptMag(cmptMultiply(d2RhoUy_[cellI], basisFuncs)));
                const scalar deltaRhoUz = cmptSum(cmptMag(cmptDivide(rhoUGrad_[cellI].z(), rDeltaXYZ_[cellI])))
                                        + cmptSum(cmptMag(cmptMultiply(d2RhoUz_[cellI], basisFuncs)));
                vector rhoU_max = cmptMag(rhoU_[cellI]) + vector(deltaRhoUx, deltaRhoUy, deltaRhoUz);
                temp = rhoE_[cellI]-deltaRhoE - 0.5*magSqr(rhoU_max)/(rho_[cellI]-deltaRho);
            }
            if (temp < 0)
            {
                count++;
                rhoGrad_[cellI]  = vector::zero;
                rhoUGrad_[cellI] = tensor::zero;
                rhoEGrad_[cellI] = vector::zero;
                d2Rho_[cellI]    = symmTensor::zero;
                d2RhoUx_[cellI] = symmTensor::zero;
                d2RhoUy_[cellI] = symmTensor::zero;
                d2RhoUz_[cellI] = symmTensor::zero;
                d2RhoE_[cellI]  = symmTensor::zero;
                d2Var_[cellI] = Mat5X6::Zero();
            }
        }
        const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
        Info << ", Negative cell [%] = " << percent;
    }
    Info << endl;
}

void Foam::evaluateEigenMatrix
(
    Mat5X5& L,
    Mat5X5& R,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const vector& normal,
    const scalar& gamma
)
{
    scalar rRho = 1.0/rho;
    vector U = rhoU*rRho;
    scalar ek = 0.5*magSqr(U);
    scalar p = (rhoE-rho*ek)*(gamma-1);
    scalar a = sqrt(p*gamma*rRho);
    scalar h = (rhoE+p)*rRho;
    scalar Vn = U&normal;
    // build Mat5x5 R
    vector an = a*normal;
    R.col(0) << 1, U.x()-an.x(), U.y()-an.y(), U.z()-an.z(), h-a*Vn;
    R.col(1) << 1, U.x()       , U.y()       , U.z()       , ek;
    R.col(2) << 1, U.x()+an.x(), U.y()+an.y(), U.z()+an.z(), h+a*Vn;
    // build Mat5x5 L
    scalar b1 = (gamma-1.0)/sqr(a), b2 = b1*ek, ra = 1.0/a;
    vector n_a = normal*ra, b1xU = b1*U;
    L.row(0) << (b2+Vn*ra), -(b1xU.x()+n_a.x()), -(b1xU.y()+n_a.y()), -(b1xU.z()+n_a.z()), b1;
    L.row(1) << -b2+1     ,   b1xU.x()         ,   b1xU.y()         ,   b1xU.z()     ,    -b1;
    L.row(2) << (b2-Vn*ra), -(b1xU.x()-n_a.x()), -(b1xU.y()-n_a.y()), -(b1xU.z()-n_a.z()), b1;
    L.row(0) *= 0.5;
    L.row(2) *= 0.5;
    if (mag(normal.x()) > 1e-8)
    {
        R.col(3) << 0,  normal.y(), -normal.x(), 0, U.x()*normal.y()-U.y()*normal.x();
        R.col(4) << 0, -normal.z(),  0, normal.x(), U.z()*normal.x()-U.x()*normal.z();
        scalar rnx = 1.0/normal.x();
        L.row(3) <<  (U.y()-Vn*normal.y())*rnx,  normal.y(),  (sqr(normal.y())-1.0)*rnx, normal.y()*normal.z()*rnx, 0;
        L.row(4) << -(U.z()-Vn*normal.z())*rnx, -normal.z(), -normal.y()*normal.z()*rnx, (1.0-sqr(normal.z()))*rnx, 0;
    }
    else if (mag(normal.y()) > 1e-8)
    {
        R.col(3) << 0, normal.y(), -normal.x(), 0, U.x()*normal.y()-U.y()*normal.x();
        R.col(4) << 0, 0, normal.z(), -normal.y(), U.y()*normal.z()-U.z()*normal.y();
        scalar rny = 1.0/normal.y();
        L.row(3) << -(U.x()-Vn*normal.x())*rny, (1.0-sqr(normal.x()))*rny, -normal.x(), normal.x()*normal.z()*rny, 0;
        L.row(4) <<  (U.z()-Vn*normal.z())*rny, normal.x()*normal.z()*rny,  normal.z(), (sqr(normal.z())-1.0)*rny, 0;
    }
    else
    {
        R.col(3) << 0, -normal.z(), 0, normal.x(), U.z()*normal.x()-U.x()*normal.z();
        R.col(4) << 0, 0, normal.z(), -normal.y(), U.y()*normal.z()-U.z()*normal.y();
        scalar rnz = 1.0/normal.z();
        L.row(3) <<  (U.x()-Vn*normal.x())*rnz,  (sqr(normal.x())-1.0)*rnz,  normal.x()*normal.y()*rnz,  normal.x(), 0;
        L.row(4) << -(U.y()-Vn*normal.y())*rnz, -normal.x()*normal.y()*rnz,  (1.0-sqr(normal.y()))*rnz, -normal.y(), 0;
    }
}

template<int cols>
Eigen::Matrix<scalar, 5, cols> Foam::Eigen_WBAP_L2_Limiter
(
    const std::vector<vector>& normalList,
    const std::vector<Eigen::Matrix<scalar, 5, cols>>& matList,
    const scalar& rho,
    const vector& rhoU,
    const scalar& rhoE,
    const scalar& gamma,
    const scalar& value
)
{
    const label n = matList.size();
    std::vector<Eigen::Matrix<scalar, 5, cols>> list;
    list.reserve(normalList.size());
    Mat5X5 L, R;
    for (const auto& normal : normalList)
    {
        evaluateEigenMatrix(L, R, rho, rhoU, rhoE, normal, gamma); 
        std::vector<Eigen::Matrix<scalar, 5, cols>> tempList;
        tempList.reserve(n);
        for (const auto& mat : matList)
        {
            tempList.emplace_back(L*mat);
        }
        list.emplace_back(R*WBAP_L2_Limiter<cols>(tempList, value));
    }
    return WBAP_L2_Limiter<cols>(list, 1.0);
}