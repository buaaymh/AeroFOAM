symmTensor Foam::basisPoly
(
    const vector& delta,
    const vector& rLengthScale,
    const symmTensor& basisConst
)
{
    const vector B = cmptMultiply(delta, rLengthScale);
    return symmTensor
    (
        0.5*(sqr(B.x())-basisConst.xx()), B.x()*B.y()-basisConst.xy(), B.x()*B.z()-basisConst.xz(),
                                     0.5*(sqr(B.y())-basisConst.yy()), B.y()*B.z()-basisConst.yz(),
                                                                   0.5*(sqr(B.z())-basisConst.zz())
    );
}

symmTensor Foam::basisPolyDn1
(
    const vector& delta,
    const vector& rLengthScale,
    const vector& normal
)
{
    const vector B = cmptMultiply(delta, rLengthScale);
    const vector rLn = cmptMultiply(rLengthScale, normal);
    return symmTensor
    (
        B.x()*rLn.x(), B.y()*rLn.x()+B.x()*rLn.y(), B.z()*rLn.x()+B.x()*rLn.z(),
                       B.y()*rLn.y(),               B.z()*rLn.y()+B.y()*rLn.z(),
                                                    B.z()*rLn.z()
    );
}

symmTensor Foam::basisPolyDn2
(
    const vector& rLengthScale,
    const vector& normal
)
{
    const vector Bn2 = cmptMultiply(rLengthScale, normal);
    return symmTensor
    (
        sqr(Bn2.x()), 2*Bn2.x()*Bn2.y(), 2*Bn2.x()*Bn2.z(),
                           sqr(Bn2.y()), 2*Bn2.y()*Bn2.z(),
                                               sqr(Bn2.z())
    );
}

Foam::vrScheme::Column Foam::basisPolyEigen
(
    const vector& delta,
    const vector& rLengthScale,
    const symmTensor& basisConst
)
{
    const vector B = cmptMultiply(delta, rLengthScale);
    return Foam::vrScheme::Column
    (
        0.5*(sqr(B.x())-basisConst.xx()), B.x()*B.y()-basisConst.xy(), B.x()*B.z()-basisConst.xz(),
                                     0.5*(sqr(B.y())-basisConst.yy()), B.y()*B.z()-basisConst.yz(),
                                                                   0.5*(sqr(B.z())-basisConst.zz())
    );
}

Foam::vrScheme::Column Foam::basisPolyDn1Eigen
(
    const vector& delta,
    const vector& rLengthScale,
    const vector& normal
)
{
    const vector B = cmptMultiply(delta, rLengthScale);
    const vector rLn = cmptMultiply(rLengthScale, normal);
    return Foam::vrScheme::Column
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
    const vector& rLengthScale,
    const symmTensor& basisConst
)
{
    symmTensor basisPoly = Foam::basisPoly(delta, rLengthScale, basisConst);
    return (grad&delta) + d2Val.xx()*basisPoly.xx() + d2Val.xy()*basisPoly.xy() + d2Val.xz()*basisPoly.xz()
                                                    + d2Val.yy()*basisPoly.yy() + d2Val.yz()*basisPoly.yz()
                                                                                + d2Val.zz()*basisPoly.zz();
}

vector Foam::vectorPoly
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rLengthScale,
    const symmTensor& basisConst
)
{
    symmTensor basisPoly = Foam::basisPoly(delta, rLengthScale, basisConst);
    return (grad&delta) + vector
    (
        d2Ux.xx()*basisPoly.xx() + d2Ux.xy()*basisPoly.xy() + d2Ux.xz()*basisPoly.xz() +
        d2Ux.yy()*basisPoly.yy() + d2Ux.yz()*basisPoly.yz() + d2Ux.zz()*basisPoly.zz(),
        d2Uy.xx()*basisPoly.xx() + d2Uy.xy()*basisPoly.xy() + d2Uy.xz()*basisPoly.xz() +
        d2Uy.yy()*basisPoly.yy() + d2Uy.yz()*basisPoly.yz() + d2Uy.zz()*basisPoly.zz(),
        d2Uz.xx()*basisPoly.xx() + d2Uz.xy()*basisPoly.xy() + d2Uz.xz()*basisPoly.xz() +
        d2Uz.yy()*basisPoly.yy() + d2Uz.yz()*basisPoly.yz() + d2Uz.zz()*basisPoly.zz()
    );
}

vector Foam::scalarPolyGrad
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rLengthScale
)
{
    const tensor gradTensor = rLengthScale * vector(cmptMultiply(delta, rLengthScale));
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
    const vector& rLengthScale
)
{
    const tensor gradTensor = rLengthScale * vector(cmptMultiply(delta, rLengthScale));
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