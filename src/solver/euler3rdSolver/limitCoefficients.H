void Foam::euler3rdSolver::limitCoefficients()
{
    // detect trouble cells
    const scalar IS = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("IS", 1.0);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    trouble_ = std::vector<bool>(mesh_.nCells(), false);
    const scalar one_gamm1 = 1.0/(fluidProps_.gamma-1.0);
    label count = 0;
    scalarField rhoE_max = rhoE_.primitiveField();
    forAll(mesh_.owner(), faceI)
    {
        const label i = owner[faceI];
        const label j = neighbour[faceI];
        rhoE_max[i] = max(rhoE_max[i], rhoE_[j]);
        rhoE_max[j] = max(rhoE_max[j], rhoE_[i]);
    }
    forAll(mesh_.cells(), cellI)
    {
        scalar error = 0.0;
        const labelList& cellFaces = mesh_.cells()[cellI];
        const scalar rho_i = rho_[cellI] + scalarPoly(vector::zero, rhoGrad_[cellI], d2Rho_[cellI], rLengthScale_[cellI], basisConst_[cellI]);
        const vector U_i = U_[cellI] + vectorPoly(vector::zero, UGrad_[cellI], d2Ux_[cellI], d2Uy_[cellI], d2Uz_[cellI], rLengthScale_[cellI], basisConst_[cellI]);
        const scalar p_i = p_[cellI] + scalarPoly(vector::zero, pGrad_[cellI], d2P_[cellI], rLengthScale_[cellI], basisConst_[cellI]);
        const scalar rhoE_i = p_i*one_gamm1 + 0.5*rho_i*magSqr(U_i);
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                const vector delta = mesh_.C()[cellI] - mesh_.C()[j];
                const scalar rho_j = rho_[j] + scalarPoly(delta, rhoGrad_[j], d2Rho_[j], rLengthScale_[j], basisConst_[j]);
                const vector U_j = U_[j] + vectorPoly(delta, UGrad_[j], d2Ux_[j], d2Uy_[j], d2Uz_[j], rLengthScale_[j], basisConst_[j]);
                const scalar p_j = p_[j] + scalarPoly(delta, pGrad_[j], d2P_[j], rLengthScale_[j], basisConst_[j]);
                error += mag(p_j*one_gamm1+0.5*rho_j*magSqr(U_j)-rhoE_i);
            }
        }
        if (error/(N_h_[cellI]*rhoE_max[cellI]) > IS)
        {
            trouble_[cellI] = true;
            count++;
        }
    }
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))
                            /scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "Trouble cell [%] = " << percent << endl;

    // WBAP limiting
    const scalar np = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("np", 10);
    d2RhoLimited_ = d2Rho_;
    d2UxLimited_  = d2Ux_;
    d2UyLimited_  = d2Uy_;
    d2UzLimited_  = d2Uz_;
    d2PLimited_   = d2P_;
    forAll(mesh_.cells(), cellI)
    {
        if (!trouble_[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<symmTensor> d2RhoList;
        std::vector<symmTensor> d2UxList;
        std::vector<symmTensor> d2UyList;
        std::vector<symmTensor> d2UzList;
        std::vector<symmTensor> d2PList;
        d2RhoList.push_back(d2Rho_[cellI]);
        d2UxList.push_back(d2Ux_[cellI]);
        d2UyList.push_back(d2Uy_[cellI]);
        d2UzList.push_back(d2Uz_[cellI]);
        d2PList.push_back(d2P_[cellI]);
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                d2RhoList.push_back(d2Rho_[j]);
                d2UxList.push_back(d2Ux_[j]);
                d2UyList.push_back(d2Uy_[j]);
                d2UzList.push_back(d2Uz_[j]);
                d2PList.push_back(d2P_[j]);
            }
        }
        d2RhoLimited_[cellI] = symmTensorWBAP(d2RhoList, np);
        d2UxLimited_[cellI]  = symmTensorWBAP(d2UxList, np);
        d2UyLimited_[cellI]  = symmTensorWBAP(d2UyList, np);
        d2UzLimited_[cellI]  = symmTensorWBAP(d2UzList, np);
        d2PLimited_[cellI]   = symmTensorWBAP(d2PList, np);
    }
    rhoGradLimited_ = rhoGrad_;
    UGradLimited_ = UGrad_;
    pGradLimited_ = pGrad_;
    forAll(mesh_.cells(), cellI)
    {
        if (!trouble_[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<vector> rhoGradList;
        std::vector<tensor> UGradList;
        std::vector<vector> pGradList;
        rhoGradList.push_back(rhoGrad_[cellI]);
        UGradList.push_back(UGrad_[cellI]);
        pGradList.push_back(pGrad_[cellI]);
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                const vector delta_ij = mesh_.C()[cellI] - mesh_.C()[j];
                const tensor temp
                (
                    delta_ij.x()*sqr(rLengthScale_[j].x()),
                    delta_ij.y()*rLengthScale_[j].x()*rLengthScale_[j].y(),
                    delta_ij.z()*rLengthScale_[j].x()*rLengthScale_[j].z(),
                    delta_ij.x()*rLengthScale_[j].y()*rLengthScale_[j].x(),
                    delta_ij.y()*sqr(rLengthScale_[j].y()),
                    delta_ij.z()*rLengthScale_[j].y()*rLengthScale_[j].z(),
                    delta_ij.x()*rLengthScale_[j].z()*rLengthScale_[j].x(),
                    delta_ij.y()*rLengthScale_[j].z()*rLengthScale_[j].y(),
                    delta_ij.z()*sqr(rLengthScale_[j].z())
                );
                const vector rhoGrad_j = rhoGrad_[j] + vector(d2RhoLimited_[j].xx()*temp.xx()+d2RhoLimited_[j].xy()*temp.xy()+d2RhoLimited_[j].xz()*temp.xz(),
                                                             d2RhoLimited_[j].xy()*temp.yx()+d2RhoLimited_[j].yy()*temp.yy()+d2RhoLimited_[j].yz()*temp.yz(),
                                                             d2RhoLimited_[j].xz()*temp.zx()+d2RhoLimited_[j].yz()*temp.zy()+d2RhoLimited_[j].zz()*temp.zz());
                const tensor UGrad_j = UGrad_[j] + tensor(d2UxLimited_[j].xx()*temp.xx()+d2UxLimited_[j].xy()*temp.xy()+d2UxLimited_[j].xz()*temp.xz(),
                                                         d2UxLimited_[j].xy()*temp.yx()+d2UxLimited_[j].yy()*temp.yy()+d2UxLimited_[j].yz()*temp.yz(),
                                                         d2UxLimited_[j].xz()*temp.zx()+d2UxLimited_[j].yz()*temp.zy()+d2UxLimited_[j].zz()*temp.zz(),
                                                         d2UyLimited_[j].xx()*temp.xx()+d2UyLimited_[j].xy()*temp.xy()+d2UyLimited_[j].xz()*temp.xz(),
                                                         d2UyLimited_[j].xy()*temp.yx()+d2UyLimited_[j].yy()*temp.yy()+d2UyLimited_[j].yz()*temp.yz(),
                                                         d2UyLimited_[j].xz()*temp.zx()+d2UyLimited_[j].yz()*temp.zy()+d2UyLimited_[j].zz()*temp.zz(),
                                                         d2UzLimited_[j].xx()*temp.xx()+d2UzLimited_[j].xy()*temp.xy()+d2UzLimited_[j].xz()*temp.xz(),
                                                         d2UzLimited_[j].xy()*temp.yx()+d2UzLimited_[j].yy()*temp.yy()+d2UzLimited_[j].yz()*temp.yz(),
                                                         d2UzLimited_[j].xz()*temp.zx()+d2UzLimited_[j].yz()*temp.zy()+d2UzLimited_[j].zz()*temp.zz());
                const vector pGrad_j = pGrad_[j] + vector(d2PLimited_[j].xx()*temp.xx()+d2PLimited_[j].xy()*temp.xy()+d2PLimited_[j].xz()*temp.xz(),
                                                         d2PLimited_[j].xy()*temp.yx()+d2PLimited_[j].yy()*temp.yy()+d2PLimited_[j].yz()*temp.yz(),
                                                         d2PLimited_[j].xz()*temp.zx()+d2PLimited_[j].yz()*temp.zy()+d2PLimited_[j].zz()*temp.zz());
                rhoGradList.push_back(rhoGrad_j);
                UGradList.push_back(UGrad_j);
                pGradList.push_back(pGrad_j);
            }
        }
        rhoGradLimited_[cellI] = vectorWBAP(rhoGradList, np);
        UGradLimited_[cellI] = tensorWBAP(UGradList, np);
        pGradLimited_[cellI] = vectorWBAP(pGradList, np);
    }
    // positive preserve
    const Switch positive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    if (positive)
    {
        forAll(mesh_.C(), cellI)
        {
            if (p0_[cellI]) continue;
            const symmTensor basisFuncs
            (
                0.5*(1-basisConst_[cellI].xx()), 1-basisConst_[cellI].xy(), 1-basisConst_[cellI].xz(),
                                           0.5*(1-basisConst_[cellI].yy()), 1-basisConst_[cellI].yz(),
                                                                       0.5*(1-basisConst_[cellI].yy())
            );
            const scalar deltaRho = cmptSum(cmptDivide(cmptMag(rhoGradLimited_[cellI]), rLengthScale_[cellI]))
                                  + cmptSum(cmptMultiply(cmptMag(d2RhoLimited_[cellI]), basisFuncs));
            const scalar deltaP   = cmptSum(cmptDivide(cmptMag(pGradLimited_[cellI]), rLengthScale_[cellI]))
                                  + cmptSum(cmptMultiply(cmptMag(d2PLimited_[cellI]), basisFuncs));
            scalar temp = min(rho_[cellI]-deltaRho, p_[cellI]-deltaP);
            if (temp < 0)
            {
                temp = min(rho_[cellI]/deltaRho, p_[cellI]/deltaP);
                rhoGradLimited_[cellI] *= temp;
                UGradLimited_[cellI] *= temp;
                pGradLimited_[cellI] *= temp;
                d2RhoLimited_[cellI] *= temp;
                d2UxLimited_[cellI] *= temp;
                d2UyLimited_[cellI] *= temp;
                d2UzLimited_[cellI] *= temp;
                d2PLimited_[cellI] *= temp;
            }
        }
    }
    rhoGradLimited_.correctBoundaryConditions();
    UGradLimited_.correctBoundaryConditions();
    pGradLimited_.correctBoundaryConditions();
    d2RhoLimited_.correctBoundaryConditions();
    d2UxLimited_.correctBoundaryConditions();
    d2UyLimited_.correctBoundaryConditions();
    d2UzLimited_.correctBoundaryConditions();
    d2PLimited_.correctBoundaryConditions();
}