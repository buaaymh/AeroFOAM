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
    vrScheme

Description
    Variational reconstruction

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef vrScheme_H
#define vrScheme_H

#include "fvCFD.H"
#include "solver.H"
#include <array>
#include <map>
#include <vector>
#include <Eigen/Dense>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class vrScheme Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vrScheme
{

public:
    using Matrix = Eigen::Matrix<scalar, 6, 6>;
    using Column = Eigen::Matrix<scalar, 6, 1>;

private:
    // Private Data
    //- Mesh data
    const fvMesh& mesh_;
    volVectorField rLengthScale_;
    volSymmTensorField basisConst_;
    std::vector<bool> p0_;
    std::vector<bool> trouble_;
    //- Linear system
    std::vector<Matrix> rA_;
    std::vector<Matrix> B_;
    std::vector<Column> bRho_;
    std::vector<Column> bUx_;
    std::vector<Column> bUy_;
    std::vector<Column> bUz_;
    std::vector<Column> bT_;
    //- Quadrature points
    std::vector<std::vector<vector>> quad_;
    scalarField delta_;
    scalarField limit_;
    scalar np_;
    scalar IS_;
    vector vrWeight_;
    Switch adaptive_;
    Switch positive_;
    // Private Member Functions
    //- Disallow default bitwise copy construct
    vrScheme(const vrScheme&);
    //- Disallow default bitwise assignment
    void operator=(const vrScheme&);

    void detectP0(const volVectorField& rhoGrad);

    void limitWBAP
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
    );

    void positiveCorrect
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
    );

public:

    // Constructors
    vrScheme(const fvMesh& mesh);

    // Destructor - default
    virtual ~vrScheme() {};

    // Member Functions

    void rLengthScaleInit();

    void basisConstInit();

    void matInit();

    void updateCoefficients
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
    );

    const volVectorField& rLengthScale() const { return rLengthScale_; }

    const volSymmTensorField& basisConst() const { return basisConst_; }

    const std::vector<vector>& quadrature(label faceI) const { return quad_[faceI]; }

    const scalarField& delta() const { return delta_; }

    const std::vector<bool>& p0() const { return p0_; }

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


/*---------------------------------------------------------------------------*\
                           Methods Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * *Taylor Basis Functions * * * * * * * * * * * * * * //

symmTensor basisPoly
(
    const vector& delta,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

symmTensor basisPolyDn1
(
    const vector& delta,
    const vector& rLengthScale,
    const vector& normal
);

symmTensor basisPolyDn2
(
    const vector& rLengthScale,
    const vector& normal
);

vrScheme::Column basisPolyEigen
(
    const vector& delta,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

vrScheme::Column basisPolyDn1Eigen
(
    const vector& delta,
    const vector& rLengthScale,
    const vector& normal
);

scalar scalarPoly
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

vector vectorPoly
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rLengthScale,
    const symmTensor& basisConst
);

vector scalarPolyGrad
(
    const vector& delta,
    const vector& grad,
    const symmTensor& d2Val,
    const vector& rLengthScale
);

tensor vectorPolyGrad
(
    const vector& delta,
    const tensor& grad,
    const symmTensor& d2Ux,
    const symmTensor& d2Uy,
    const symmTensor& d2Uz,
    const vector& rLengthScale
);

// * * * * * * * * * * * * * * * WAP Functions * * * * * * * * * * * * * * //

vector vectorWAP
(
    const std::vector<vector>& list,
    scalar value
);

symmTensor symmTensorWAP
(
    const std::vector<symmTensor>& list,
    scalar value
);

tensor tensorWAP
(
    const std::vector<tensor>& list,
    scalar value
);

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //