void Foam::vrScheme::detectP0(const volVectorField& rhoGrad)
{
    p0_ = std::vector<bool>(mesh_.nCells(), false);
    label count = 0;
    forAll(mesh_.C(), cellI)
    {
        if (cmptSum(cmptDivide(cmptMag(rhoGrad[cellI]), rLengthScale_[cellI])) < 1e-6)
        {
            p0_[cellI] = true;
            count++;
        }
    }
}

void Foam::vrScheme::limitWBAP
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    volVectorField& rhoGrad,
    volTensorField& UGrad,
    volVectorField& TGrad,
    volSymmTensorField& d2Rho,
    volSymmTensorField& d2Ux,
    volSymmTensorField& d2Uy,
    volSymmTensorField& d2Uz,
    volSymmTensorField& d2T,
    const scalar gamma
)
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    trouble_ = std::vector<bool>(mesh_.nCells(), false);
    const scalar temp = 1.0/((gamma-1.0)*gamma);
    label count = 0;
    forAll(mesh_.cells(), cellI)
    {
        scalar error = 0.0;
        const labelList& cellFaces = mesh_.cells()[cellI];
        const scalar rho_i = rho[cellI] + scalarPoly(vector::zero, rhoGrad[cellI], d2Rho[cellI], rLengthScale_[cellI], basisConst_[cellI]);
        const vector U_i = U[cellI] + vectorPoly(vector::zero, UGrad[cellI], d2Ux[cellI], d2Uy[cellI], d2Uz[cellI], rLengthScale_[cellI], basisConst_[cellI]);
        const scalar T_i = T[cellI] + scalarPoly(vector::zero, TGrad[cellI], d2T[cellI], rLengthScale_[cellI], basisConst_[cellI]);
        const scalar rhoE_i = rho_i*(T_i*temp+0.5*magSqr(U_i));
        const scalar rhoE_0 = rho[cellI]*(T[cellI]*temp+0.5*magSqr(U[cellI]));
        forAll(cellFaces,f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                const vector delta = mesh_.C()[cellI] - mesh_.C()[j];
                const scalar rho_j = rho[j] + scalarPoly(delta, rhoGrad[j], d2Rho[j], rLengthScale_[j], basisConst_[j]);
                const vector U_j = U[j] + vectorPoly(delta, UGrad[j], d2Ux[j], d2Uy[j], d2Uz[j], rLengthScale_[j], basisConst_[j]);
                const scalar T_j = T[j] + scalarPoly(delta, TGrad[j], d2T[j], rLengthScale_[j], basisConst_[j]);
                error += mag(rho_j*(T_j*temp+0.5*magSqr(U_j))-rhoE_i)*mesh_.magSf()[faceI];
            }
        }
        if (error/(mesh_.V()[cellI]*rhoE_0) > IS_)
        {
            trouble_[cellI] = true;
            count++;
        }
    }
    const scalar percent = 100*scalar(returnReduce(count, sumOp<label>()))
                            /scalar(returnReduce(mesh_.nCells(), sumOp<label>()));
    Info << "Trouble cell [%] = " << percent << endl;
    std::vector<symmTensor> d2RhoCorrected(count);
    std::vector<symmTensor> d2UxCorrected(count);
    std::vector<symmTensor> d2UyCorrected(count);
    std::vector<symmTensor> d2UzCorrected(count);
    std::vector<symmTensor> d2TCorrected(count);
    count = 0;

    // Limiting procedure
    // Secondary order
    forAll(mesh_.cells(), cellI)
    {
        if (!trouble_[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<symmTensor> d2RhoList;
        std::vector<symmTensor> d2UxList;
        std::vector<symmTensor> d2UyList;
        std::vector<symmTensor> d2UzList;
        std::vector<symmTensor> d2TList;
        d2RhoList.push_back(d2Rho[cellI]);
        d2UxList.push_back(d2Ux[cellI]);
        d2UyList.push_back(d2Uy[cellI]);
        d2UzList.push_back(d2Uz[cellI]);
        d2TList.push_back(d2T[cellI]);
        forAll(cellFaces, f)
        {
            const label faceI = cellFaces[f];
            if (mesh_.isInternalFace(faceI))
            {
                label j = owner[faceI];
                if (j == cellI) { j = neighbour[faceI]; }
                d2RhoList.push_back(d2Rho[j]);
                d2UxList.push_back(d2Ux[j]);
                d2UyList.push_back(d2Uy[j]);
                d2UzList.push_back(d2Uz[j]);
                d2TList.push_back(d2T[j]);
            }
        }
        d2RhoCorrected[count] = symmTensorWAP(d2RhoList, np_);
        d2UxCorrected[count] = symmTensorWAP(d2UxList, np_);
        d2UyCorrected[count] = symmTensorWAP(d2UyList, np_);
        d2UzCorrected[count] = symmTensorWAP(d2UzList, np_);
        d2TCorrected[count++] = symmTensorWAP(d2TList, np_);
    }
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!trouble_[cellI]) continue;
        d2Rho[cellI] = d2RhoCorrected[count];
        d2Ux[cellI] = d2UxCorrected[count];
        d2Uy[cellI] = d2UyCorrected[count];
        d2Uz[cellI] = d2UzCorrected[count];
        d2T[cellI] = d2TCorrected[count++];
    }
    std::vector<vector> rhoGradCorrected(count);
    std::vector<tensor> UGradCorrected(count);
    std::vector<vector> TGradCorrected(count);
    count = 0;

    forAll(mesh_.cells(), cellI)
    {
        if (!trouble_[cellI]) continue;
        const labelList& cellFaces = mesh_.cells()[cellI];
        std::vector<vector> rhoGradList;
        std::vector<tensor> UGradList;
        std::vector<vector> TGradList;
        rhoGradList.push_back(rhoGrad[cellI]);
        UGradList.push_back(UGrad[cellI]);
        TGradList.push_back(TGrad[cellI]);
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
                const vector rhoGrad_j = rhoGrad[j] + vector(d2Rho[j].xx()*temp.xx()+d2Rho[j].xy()*temp.xy()+d2Rho[j].xz()*temp.xz(),
                                                             d2Rho[j].xy()*temp.yx()+d2Rho[j].yy()*temp.yy()+d2Rho[j].yz()*temp.yz(),
                                                             d2Rho[j].xz()*temp.zx()+d2Rho[j].yz()*temp.zy()+d2Rho[j].zz()*temp.zz());
                const tensor UGrad_j = UGrad[j] + tensor(d2Ux[j].xx()*temp.xx()+d2Ux[j].xy()*temp.xy()+d2Ux[j].xz()*temp.xz(),
                                                         d2Ux[j].xy()*temp.yx()+d2Ux[j].yy()*temp.yy()+d2Ux[j].yz()*temp.yz(),
                                                         d2Ux[j].xz()*temp.zx()+d2Ux[j].yz()*temp.zy()+d2Ux[j].zz()*temp.zz(),
                                                         d2Uy[j].xx()*temp.xx()+d2Uy[j].xy()*temp.xy()+d2Uy[j].xz()*temp.xz(),
                                                         d2Uy[j].xy()*temp.yx()+d2Uy[j].yy()*temp.yy()+d2Uy[j].yz()*temp.yz(),
                                                         d2Uy[j].xz()*temp.zx()+d2Uy[j].yz()*temp.zy()+d2Uy[j].zz()*temp.zz(),
                                                         d2Uz[j].xx()*temp.xx()+d2Uz[j].xy()*temp.xy()+d2Uz[j].xz()*temp.xz(),
                                                         d2Uz[j].xy()*temp.yx()+d2Uz[j].yy()*temp.yy()+d2Uz[j].yz()*temp.yz(),
                                                         d2Uz[j].xz()*temp.zx()+d2Uz[j].yz()*temp.zy()+d2Uz[j].zz()*temp.zz());
                const vector TGrad_j = TGrad[j] + vector(d2T[j].xx()*temp.xx()+d2T[j].xy()*temp.xy()+d2T[j].xz()*temp.xz(),
                                                         d2T[j].xy()*temp.yx()+d2T[j].yy()*temp.yy()+d2T[j].yz()*temp.yz(),
                                                         d2T[j].xz()*temp.zx()+d2T[j].yz()*temp.zy()+d2T[j].zz()*temp.zz());
                rhoGradList.push_back(rhoGrad_j);
                UGradList.push_back(UGrad_j);
                TGradList.push_back(TGrad_j);
            }
        }
        rhoGradCorrected[count] = vectorWAP(rhoGradList, np_);
        UGradCorrected[count] = tensorWAP(UGradList, np_);
        TGradCorrected[count++] = vectorWAP(TGradList, np_);
    }
    count = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (!trouble_[cellI]) continue;
        rhoGrad[cellI] = rhoGradCorrected[count];
        UGrad[cellI] = UGradCorrected[count];
        TGrad[cellI] = TGradCorrected[count++];
    }
}

void Foam::vrScheme::positiveCorrect
(
    const volScalarField& rho,
    const volVectorField& U,
    const volScalarField& T,
    volVectorField& rhoGrad,
    volTensorField& UGrad,
    volVectorField& TGrad,
    volSymmTensorField& d2Rho,
    volSymmTensorField& d2Ux,
    volSymmTensorField& d2Uy,
    volSymmTensorField& d2Uz,
    volSymmTensorField& d2T
)
{
    forAll(mesh_.C(), cellI)
    {
        if (p0_[cellI] || !trouble_[cellI]) continue;
        const symmTensor basisFuncs
        (
            0.5*(1-basisConst_[cellI].xx()), 1-basisConst_[cellI].xy(), 1-basisConst_[cellI].xz(),
                                       0.5*(1-basisConst_[cellI].yy()), 1-basisConst_[cellI].yz(),
                                                                   0.5*(1-basisConst_[cellI].yy())
        );
        const scalar deltaRho = cmptSum(cmptDivide(cmptMag(rhoGrad[cellI]), rLengthScale_[cellI]))
                              + cmptSum(cmptMultiply(cmptMag(d2Rho[cellI]), basisFuncs));
        const scalar deltaT   = cmptSum(cmptDivide(cmptMag(TGrad[cellI]), rLengthScale_[cellI]))
                              + cmptSum(cmptMultiply(cmptMag(d2T[cellI]), basisFuncs));
        scalar temp = min(rho[cellI]-deltaRho, T[cellI]-deltaT);
        if (temp < 0)
        {
            temp = min(rho[cellI]/deltaRho, T[cellI]/deltaT);
            rhoGrad[cellI] *= temp;
            UGrad[cellI] *= temp;
            TGrad[cellI] *= temp;
            d2Rho[cellI] *= temp;
            d2Ux[cellI] *= temp;
            d2Uy[cellI] *= temp;
            d2Uz[cellI] *= temp;
            d2T[cellI] *= temp;
        }
    }
}