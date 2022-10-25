void Foam::eulerPrimVar3rdSolver::limitCoefficients()
{
    // detect trouble cells
    const scalar IS = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("IS", 1.0);
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField delta = mesh_.delta();
    scalarField rho_max(rho_.primitiveField());
    scalarField rho_mid(mesh_.nCells());
    scalarField error(mesh_.nCells(), 0.0);
    forAll(mesh_.cells(), cellI)
    {
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
    forAll(mesh_.boundary(), patchI)
    {
        if (mesh_.boundary()[patchI].coupled())
        {
            const UList<label> &bfaceCells = mesh_.boundary()[patchI].faceCells();
            const vectorField delta = mesh_.boundary()[patchI].delta();
            const vectorField rDeltaXYZ_neigh = rDeltaXYZ_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField basisMean_neigh = basisMean_.boundaryField()[patchI].patchNeighbourField();
            const scalarField rho_neigh = rho_.boundaryField()[patchI].patchNeighbourField();
            const vectorField rhoGrad_neigh = rhoGrad_.boundaryField()[patchI].patchNeighbourField();
            const symmTensorField d2Rho_neigh = d2Rho_.boundaryField()[patchI].patchNeighbourField();
            forAll(bfaceCells, j)
            {
                const label i = bfaceCells[j];
                rho_max[i] = max(rho_max[i], rho_neigh[j]);
                error[i] += mag(rho_neigh[j]+scalarPoly(-delta[j], rhoGrad_neigh[j], d2Rho_neigh[j], rDeltaXYZ_neigh[j], basisMean_neigh[j])-rho_mid[i]);
            }
        }
    }
    label count = 0;
    isTrouble_ = std::vector<bool>(mesh_.nCells(), false);
    forAll(mesh_.cells(), cellI)
    {
        if (error[cellI]/(N_h_[cellI]*rho_max[cellI]) > IS)
        {
            isTrouble_[cellI] = true;
            count++;
        }
    }
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))/scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "Trouble cell [%] = " << percent << endl;

    // WBAP limiting
    const scalar np = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<scalar>("np", 10);
    d2RhoLimited_ = d2Rho_;
    d2UxLimited_  = d2Ux_;
    d2UyLimited_  = d2Uy_;
    d2UzLimited_  = d2Uz_;
    d2TLimited_   = d2T_;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble_[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<symmTensor> d2RhoList;
        std::vector<symmTensor> d2UxList;
        std::vector<symmTensor> d2UyList;
        std::vector<symmTensor> d2UzList;
        std::vector<symmTensor> d2TList;
        d2RhoList.push_back(d2Rho_[cellI]);
        d2UxList.push_back(d2Ux_[cellI]);
        d2UyList.push_back(d2Uy_[cellI]);
        d2UzList.push_back(d2Uz_[cellI]);
        d2TList.push_back(d2T_[cellI]);
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
                d2TList.push_back(d2T_[j]);
            }
        }
        d2RhoLimited_[cellI] = symmTensorWBAP(d2RhoList, np);
        d2UxLimited_[cellI]  = symmTensorWBAP(d2UxList, np);
        d2UyLimited_[cellI]  = symmTensorWBAP(d2UyList, np);
        d2UzLimited_[cellI]  = symmTensorWBAP(d2UzList, np);
        d2TLimited_[cellI]   = symmTensorWBAP(d2TList, np);
    }
    rhoGradLimited_ = rhoGrad_;
    UGradLimited_ = UGrad_;
    TGradLimited_ = TGrad_;
    forAll(mesh_.cells(), cellI)
    {
        if (!isTrouble_[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<vector> rhoGradList;
        std::vector<tensor> UGradList;
        std::vector<vector> TGradList;
        rhoGradList.push_back(rhoGrad_[cellI]);
        UGradList.push_back(UGrad_[cellI]);
        TGradList.push_back(TGrad_[cellI]);
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
                    delta_ij.x()*sqr(rDeltaXYZ_[j].x()),
                    delta_ij.y()*rDeltaXYZ_[j].x()*rDeltaXYZ_[j].y(),
                    delta_ij.z()*rDeltaXYZ_[j].x()*rDeltaXYZ_[j].z(),
                    delta_ij.x()*rDeltaXYZ_[j].y()*rDeltaXYZ_[j].x(),
                    delta_ij.y()*sqr(rDeltaXYZ_[j].y()),
                    delta_ij.z()*rDeltaXYZ_[j].y()*rDeltaXYZ_[j].z(),
                    delta_ij.x()*rDeltaXYZ_[j].z()*rDeltaXYZ_[j].x(),
                    delta_ij.y()*rDeltaXYZ_[j].z()*rDeltaXYZ_[j].y(),
                    delta_ij.z()*sqr(rDeltaXYZ_[j].z())
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
                const vector TGrad_j = TGrad_[j] + vector(d2TLimited_[j].xx()*temp.xx()+d2TLimited_[j].xy()*temp.xy()+d2TLimited_[j].xz()*temp.xz(),
                                                          d2TLimited_[j].xy()*temp.yx()+d2TLimited_[j].yy()*temp.yy()+d2TLimited_[j].yz()*temp.yz(),
                                                          d2TLimited_[j].xz()*temp.zx()+d2TLimited_[j].yz()*temp.zy()+d2TLimited_[j].zz()*temp.zz());
                rhoGradList.push_back(rhoGrad_j);
                UGradList.push_back(UGrad_j);
                TGradList.push_back(TGrad_j);
            }
        }
        rhoGradLimited_[cellI] = vectorWBAP(rhoGradList, np);
        UGradLimited_[cellI] = tensorWBAP(UGradList, np);
        TGradLimited_[cellI] = vectorWBAP(TGradList, np);
    }
    // positive preserve
    const Switch positive = mesh_.schemesDict().subDict("vrSchemes").lookupOrDefault<Switch>("positive", false);
    if (positive)
    {
        forAll(mesh_.C(), cellI)
        {
            if (isP0Cell_[cellI]) continue;
            const symmTensor basisFuncs
            (
                0.5*(1-basisMean_[cellI].xx()), 1-basisMean_[cellI].xy(), 1-basisMean_[cellI].xz(),
                                           0.5*(1-basisMean_[cellI].yy()), 1-basisMean_[cellI].yz(),
                                                                       0.5*(1-basisMean_[cellI].yy())
            );
            const scalar deltaRho = cmptSum(cmptDivide(cmptMag(rhoGradLimited_[cellI]), rDeltaXYZ_[cellI]))
                                  + cmptSum(cmptMultiply(cmptMag(d2RhoLimited_[cellI]), basisFuncs));
            const scalar deltaT   = cmptSum(cmptDivide(cmptMag(TGradLimited_[cellI]), rDeltaXYZ_[cellI]))
                                  + cmptSum(cmptMultiply(cmptMag(d2TLimited_[cellI]), basisFuncs));
            scalar temp = min(rho_[cellI]-deltaRho, T_[cellI]-deltaT);
            if (temp < 0)
            {
                temp = min(rho_[cellI]/deltaRho, T_[cellI]/deltaT);
                rhoGradLimited_[cellI] *= temp;
                UGradLimited_[cellI] *= temp;
                TGradLimited_[cellI] *= temp;
                d2RhoLimited_[cellI] *= temp;
                d2UxLimited_[cellI] *= temp;
                d2UyLimited_[cellI] *= temp;
                d2UzLimited_[cellI] *= temp;
                d2TLimited_[cellI] *= temp;
            }
        }
    }
    rhoGradLimited_.correctBoundaryConditions();
    UGradLimited_.correctBoundaryConditions();
    TGradLimited_.correctBoundaryConditions();
    d2RhoLimited_.correctBoundaryConditions();
    d2UxLimited_.correctBoundaryConditions();
    d2UyLimited_.correctBoundaryConditions();
    d2UzLimited_.correctBoundaryConditions();
    d2TLimited_.correctBoundaryConditions();
}