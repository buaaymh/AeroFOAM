symmTensor Foam::basisPoly
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    return symmTensor
    (
        0.5*(sqr(B.x())-basisMean.xx()), B.x()*B.y()-basisMean.xy(),  B.x()*B.z()-basisMean.xz(),
                                     0.5*(sqr(B.y())-basisMean.yy()), B.y()*B.z()-basisMean.yz(),
                                                                  0.5*(sqr(B.z())-basisMean.zz())
    );
}

symmTensor Foam::basisPolyDn1
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const vector& normal
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    const vector rLn = cmptMultiply(rDeltaXYZ, normal);
    return symmTensor
    (
        B.x()*rLn.x(), B.y()*rLn.x()+B.x()*rLn.y(), B.z()*rLn.x()+B.x()*rLn.z(),
                       B.y()*rLn.y(),               B.z()*rLn.y()+B.y()*rLn.z(),
                                                    B.z()*rLn.z()
    );
}

symmTensor Foam::basisPolyDn2
(
    const vector& rDeltaXYZ,
    const vector& normal
)
{
    const vector Bn2 = cmptMultiply(rDeltaXYZ, normal);
    return symmTensor
    (
        sqr(Bn2.x()), 2*Bn2.x()*Bn2.y(), 2*Bn2.x()*Bn2.z(),
                           sqr(Bn2.y()), 2*Bn2.y()*Bn2.z(),
                                               sqr(Bn2.z())
    );
}

Foam::Col6X1 Foam::basisPolyEigen
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    return Col6X1
    (
        0.5*(sqr(B.x())-basisMean.xx()), B.x()*B.y()-basisMean.xy(),  B.x()*B.z()-basisMean.xz(),
                                     0.5*(sqr(B.y())-basisMean.yy()), B.y()*B.z()-basisMean.yz(),
                                                                  0.5*(sqr(B.z())-basisMean.zz())
    );
}

Foam::Col6X1 Foam::basisPolyDn1Eigen
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const vector& normal
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    const vector rLn = cmptMultiply(rDeltaXYZ, normal);
    return Col6X1
    (
        B.x()*rLn.x(), B.y()*rLn.x()+B.x()*rLn.y(), B.z()*rLn.x()+B.x()*rLn.z(),
                       B.y()*rLn.y(),               B.z()*rLn.y()+B.y()*rLn.z(),
                                                    B.z()*rLn.z()
    );
}

scalar Foam::scalarPoly
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
)
{
    symmTensor basisPoly = Foam::basisPoly(delta, rDeltaXYZ, basisMean);
    return (grad&delta) + cmptSum(cmptMultiply(d2Val, basisPoly));
}

vector Foam::vectorPoly
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
)
{
    symmTensor basisPoly = Foam::basisPoly(delta, rDeltaXYZ, basisMean);
    return (grad&delta) + vector
    (
        cmptSum(cmptMultiply(d2Ux, basisPoly)),
        cmptSum(cmptMultiply(d2Uy, basisPoly)),
        cmptSum(cmptMultiply(d2Uz, basisPoly))
    );
}

vector Foam::scalarPolyGrad
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rDeltaXYZ
)
{
    const tensor gradTensor = rDeltaXYZ * vector(cmptMultiply(delta, rDeltaXYZ));
    return grad + vector(gradTensor.x()&vector(d2Val.xx(), d2Val.xy(), d2Val.xz()),
                         gradTensor.y()&vector(d2Val.xy(), d2Val.yy(), d2Val.yz()),
                         gradTensor.z()&vector(d2Val.xz(), d2Val.yz(), d2Val.zz()));
}

tensor Foam::vectorPolyGrad
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rDeltaXYZ
)
{
    const tensor gradTensor = rDeltaXYZ * vector(cmptMultiply(delta, rDeltaXYZ));
    return grad + tensor
    (
        gradTensor.x()&vector(d2Ux.xx(), d2Ux.xy(), d2Ux.xz()),
        gradTensor.y()&vector(d2Ux.xy(), d2Ux.yy(), d2Ux.yz()),
        gradTensor.z()&vector(d2Ux.xz(), d2Ux.yz(), d2Ux.zz()),
        gradTensor.x()&vector(d2Uy.xx(), d2Uy.xy(), d2Uy.xz()),
        gradTensor.y()&vector(d2Uy.xy(), d2Uy.yy(), d2Uy.yz()),
        gradTensor.z()&vector(d2Uy.xz(), d2Uy.yz(), d2Uy.zz()),
        gradTensor.x()&vector(d2Uz.xx(), d2Uz.xy(), d2Uz.xz()),
        gradTensor.y()&vector(d2Uz.xy(), d2Uz.yy(), d2Uz.yz()),
        gradTensor.z()&vector(d2Uz.xz(), d2Uz.yz(), d2Uz.zz())
    );
}

Foam::Mat5X3 Foam::evaluateGrad
(
    const vector& delta,
    const Mat5X3& d1Var,
    const Mat5X6& d2Var,
    const vector& rDeltaXYZ
)
{
    const tensor gradTensor = rDeltaXYZ * vector(cmptMultiply(delta, rDeltaXYZ));
    Mat5X3 d2VarGrad;
    d2VarGrad << (gradTensor.x()&vector(d2Var(0,0), d2Var(0,1), d2Var(0,2))),
                 (gradTensor.y()&vector(d2Var(0,1), d2Var(0,3), d2Var(0,4))),
                 (gradTensor.z()&vector(d2Var(0,2), d2Var(0,4), d2Var(0,5))),

                 (gradTensor.x()&vector(d2Var(1,0), d2Var(1,1), d2Var(1,2))),
                 (gradTensor.y()&vector(d2Var(1,1), d2Var(1,3), d2Var(1,4))),
                 (gradTensor.z()&vector(d2Var(1,2), d2Var(1,4), d2Var(1,5))),

                 (gradTensor.x()&vector(d2Var(2,0), d2Var(2,1), d2Var(2,2))),
                 (gradTensor.y()&vector(d2Var(2,1), d2Var(2,3), d2Var(2,4))),
                 (gradTensor.z()&vector(d2Var(2,2), d2Var(2,4), d2Var(2,5))),

                 (gradTensor.x()&vector(d2Var(3,0), d2Var(3,1), d2Var(3,2))),
                 (gradTensor.y()&vector(d2Var(3,1), d2Var(3,3), d2Var(3,4))),
                 (gradTensor.z()&vector(d2Var(3,2), d2Var(3,4), d2Var(3,5))),
                 
                 (gradTensor.x()&vector(d2Var(4,0), d2Var(4,1), d2Var(4,2))),
                 (gradTensor.y()&vector(d2Var(4,1), d2Var(4,3), d2Var(4,4))),
                 (gradTensor.z()&vector(d2Var(4,2), d2Var(4,4), d2Var(4,5)));
    return d2VarGrad + d1Var;
}