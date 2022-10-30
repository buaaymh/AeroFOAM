/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    vr3rdSolver

Description
    Variational reconstruction

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef vr3rdSolver_H
#define vr3rdSolver_H

#include "solver.H"
#include <array>
#include <map>
#include <vector>
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

using Col3X1 = Eigen::Matrix<scalar, 3, 1>;
using Col6X1 = Eigen::Matrix<scalar, 6, 1>;
using Col5X1 = Eigen::Matrix<scalar, 5, 1>;
using Mat5X3 = Eigen::Matrix<scalar, 5, 3>;
using Mat5X6 = Eigen::Matrix<scalar, 5, 6>;
using Mat6X5 = Eigen::Matrix<scalar, 6, 5>;
using Mat6X6 = Eigen::Matrix<scalar, 6, 6>;

/*---------------------------------------------------------------------------*\
                           Class vr3rdSolver Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vr3rdSolver : public solver
{
    // Private Member Functions
    //- Disallow default bitwise copy construct
    vr3rdSolver(const vr3rdSolver&);
    //- Disallow default bitwise assignment
    void operator=(const vr3rdSolver&);

public:

    // Constructors
    vr3rdSolver
    (
        const fluidProperties& fluidProps,
        volScalarField& rho,
        volVectorField& U,
        volScalarField& p
    );

    // Destructor - default
    virtual ~vr3rdSolver() {};

    // Member Functions
    virtual void evaluateFlowRes
    (
        scalarField& resRho,
        vectorField& resRhoU,
        scalarField& resRhoE
    ) = 0;

protected:
    // Private Data
    std::vector<Mat5X3> d1Var_;
    std::vector<Mat5X6> d2Var_;
    //- Mesh data
    volVectorField rDeltaXYZ_;
    volSymmTensorField basisMean_;
    scalarField N_h_;
    std::vector<bool> isP0Cell_;
    //- Linear system
    std::vector<Mat6X6> rA_;
    std::vector<Mat6X6> B_;
    std::map<label, std::vector<Mat6X6>> coupledB_;
    //- Quadrature points
    std::vector<std::vector<vector>> quad_;

    // Member Functions
    void initVrLinearSystem();

    void evaluateVars
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
    ) const;

    // * * * * * * * * * * * * Gauss Points * * * * * * * ** * * * * * //

    void gaussQuad4
    (
        label faceI,
        std::vector<vector>& quadPoints
    );
    
    void gaussQuad9
    (
        label faceI,
        std::vector<vector>& quadPoints
    );

    void gaussHex8
    (
        label cellI,
        std::vector<vector>& quadPoints
    );
    
};

/*---------------------------------------------------------------------------*\
                           Methods Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * ** * * * Taylor Basis Functions * * * * * ** * * * * //

template<typename T>
T basisPoly
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    return T
    (
        0.5*(sqr(B.x())-basisMean.xx()), B.x()*B.y()-basisMean.xy(),  B.x()*B.z()-basisMean.xz(),
                                     0.5*(sqr(B.y())-basisMean.yy()), B.y()*B.z()-basisMean.yz(),
                                                                  0.5*(sqr(B.z())-basisMean.zz())
    );
}

template<typename T>
T basisPolyDn1
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const vector& normal
)
{
    const vector B = cmptMultiply(delta, rDeltaXYZ);
    const vector rLn = cmptMultiply(rDeltaXYZ, normal);
    return T
    (
        B.x()*rLn.x(), B.y()*rLn.x()+B.x()*rLn.y(), B.z()*rLn.x()+B.x()*rLn.z(),
                       B.y()*rLn.y(),               B.z()*rLn.y()+B.y()*rLn.z(),
                                                    B.z()*rLn.z()
    );
}

template<typename T>
T basisPolyDn2
(
    const vector& rDeltaXYZ,
    const vector& normal
)
{
    const vector Bn2 = cmptMultiply(rDeltaXYZ, normal);
    return T
    (
        sqr(Bn2.x()), 2*Bn2.x()*Bn2.y(), 2*Bn2.x()*Bn2.z(),
                           sqr(Bn2.y()), 2*Bn2.y()*Bn2.z(),
                                               sqr(Bn2.z())
    );
}

scalar scalarPoly
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
);

vector vectorPoly
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
);

vector scalarPolyGrad
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rDeltaXYZ
);

tensor vectorPolyGrad
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rDeltaXYZ
);

Mat5X3 evaluateGrad
(
    const vector& delta,
    const Mat5X3& d1Var,
    const Mat5X6& d2Var,
    const vector& rDeltaXYZ
);

// * * * * * * * * * * * * WBAP Functions * * * * * * * * * * * * //
template<int rows, int cols>
Eigen::Matrix<scalar, rows, cols> WBAP_L2_Limiter
(
    const std::vector<Eigen::Matrix<scalar, rows, cols>>& list,
    scalar value
)
{
    const label n = list.size();
    Eigen::Matrix<scalar, rows, cols> result = list[0];
    for (label i = 0; i != rows; i++)
    {
        for (label j = 0; j != cols; j++)
        {
            if (mag(list[0](i, j)) < 1e-10) continue;
            scalar temp3 = value;
            scalar temp4 = value;
            for (label k = 1; k != n; k++)
            {
                const scalar theta = list[k](i, j)/list[0](i, j);
                if (theta <= 1e-6)
                {
                    temp3 = 0;
                    break;
                }
                temp3 += 1.0 / Foam::pow3(theta);
                temp4 += 1.0 / Foam::pow4(theta);
            }
            result(i, j) *= temp3 / temp4;
        }
    }
    return result;
}

/*---------------------------------------------------------------------------*\
                           Gauss Namespace  Declaration
\*---------------------------------------------------------------------------*/

namespace Quad4
{
    static constexpr std::array<scalar, 4> x{-0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, +0.577350269189626};
    static constexpr std::array<scalar, 4> y{+0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, -0.577350269189626};
    static constexpr std::array<scalar, 4> w{0.25, 0.25, 0.25, 0.25};
}

namespace Quad9
{
    static constexpr std::array<scalar, 9> x{-0.774596669241483,                0.0, +0.774596669241483,
                                             -0.774596669241483,                0.0, +0.774596669241483,
                                             -0.774596669241483,                0.0, +0.774596669241483};
    static constexpr std::array<scalar, 9> y{+0.774596669241483, +0.774596669241483, +0.774596669241483,
                                                            0.0,                0.0,                0.0,
                                             -0.774596669241483, -0.774596669241483, -0.774596669241483};
    static constexpr std::array<scalar, 9> w{0.30864197530864197531, 0.49382716049382716049, 0.30864197530864197531,
                                             0.49382716049382716049, 0.79012345679012345679, 0.49382716049382716049,
                                             0.30864197530864197531, 0.49382716049382716049, 0.30864197530864197531};
}

namespace Hex8
{
    static constexpr std::array<scalar, 8> x{-0.577350269189626, +0.577350269189626, -0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, +0.577350269189626, -0.577350269189626, +0.577350269189626};
    static constexpr std::array<scalar, 8> y{+0.577350269189626, +0.577350269189626, -0.577350269189626, -0.577350269189626,
                                             +0.577350269189626, +0.577350269189626, -0.577350269189626, -0.577350269189626};
    static constexpr std::array<scalar, 8> z{+0.577350269189626, +0.577350269189626, +0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, -0.577350269189626, -0.577350269189626, -0.577350269189626};
    static constexpr std::array<scalar, 8> w{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //