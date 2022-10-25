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

    using Col6X1 = Eigen::Matrix<scalar, 6, 1>;
    using Col5X1 = Eigen::Matrix<scalar, 5, 1>;
    using Mat6X5 = Eigen::Matrix<scalar, 6, 5>;
    using Mat5X6 = Eigen::Matrix<scalar, 5, 6>;
    using Mat6X6 = Eigen::Matrix<scalar, 6, 6>;

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
    //- Mesh data
    volVectorField rDeltaXYZ_;
    volSymmTensorField basisMean_;
    scalarField N_h_;
    std::vector<bool> isP0Cell_;
    std::vector<bool> isTrouble_;
    //- Linear system
    std::vector<Mat6X6> rA_;
    std::vector<Mat6X6> B_;
    std::map<label, std::vector<Mat6X6>> coupledB_;
    //- Quadrature points
    std::vector<std::vector<vector>> quad_;

    // Member Functions
    void initVrLinearSystem();

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

// * * * * * * * * * * * *Taylor Basis Functions * * * * * * * * * * * * * * //

symmTensor basisPoly
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
);

symmTensor basisPolyDn1
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const vector& normal
);

symmTensor basisPolyDn2
(
    const vector& rDeltaXYZ,
    const vector& normal
);

vr3rdSolver::Col6X1 basisPolyEigen
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const symmTensor& basisMean
);

vr3rdSolver::Col6X1 basisPolyDn1Eigen
(
    const vector& delta,
    const vector& rDeltaXYZ,
    const vector& normal
);

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

// * * * * * * * * * * * * * * * WBAP Functions * * * * * * * * * * * * * * //

vector vectorWBAP
(
    const std::vector<vector>& list,
    scalar value
);

symmTensor symmTensorWBAP
(
    const std::vector<symmTensor>& list,
    scalar value
);

tensor tensorWBAP
(
    const std::vector<tensor>& list,
    scalar value
);

/*---------------------------------------------------------------------------*\
                           Gauss Namespace  Declaration
\*---------------------------------------------------------------------------*/

namespace Quad4
{
    static constexpr std::array<scalar, 4> x{-0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, +0.577350269189626};
    static constexpr std::array<scalar, 4> y{+0.577350269189626, +0.577350269189626,
                                             -0.577350269189626, -0.577350269189626};
    static constexpr std::array<scalar, 4> w{1.0, 1.0, 1.0, 1.0};
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