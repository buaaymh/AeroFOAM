scalar Foam::scalarPoly
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
)
{
    symmTensor basis = basisPoly<symmTensor>(delta, rDeltaXYZ, basisMean);
    return (grad&delta) + cmptSum(cmptMultiply(d2Val, basis));
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
    symmTensor basis = basisPoly<symmTensor>(delta, rDeltaXYZ, basisMean);
    return (grad&delta) + vector
    (
        cmptSum(cmptMultiply(d2Ux, basis)),
        cmptSum(cmptMultiply(d2Uy, basis)),
        cmptSum(cmptMultiply(d2Uz, basis))
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

void Foam::vr3rdSolver::evaluateVars
(
    scalar& rho,
    vector& U,
    scalar& T,
    const scalar& rho_0,
    const vector& U_0,
    const scalar& T_0,
    const vector& rhoGrad,
    const tensor& UGrad,
    const vector& TGrad,
    const symmTensor& d2Rho,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const symmTensor& d2T,
    const symmTensor& basisMean,
    const vector& rDeltaXYZ,
    const vector& delta
) const
{
    symmTensor basis = basisPoly<symmTensor>(delta, rDeltaXYZ, basisMean);
    rho = rho_0 + (rhoGrad&delta) + cmptSum(cmptMultiply(d2Rho, basis));
    U   = U_0   + (UGrad  &delta) + vector(cmptSum(cmptMultiply(d2Ux, basis)),
                                           cmptSum(cmptMultiply(d2Uy, basis)),
                                           cmptSum(cmptMultiply(d2Uz, basis)));
    T   = T_0   + (TGrad  &delta) + cmptSum(cmptMultiply(d2T, basis));
}