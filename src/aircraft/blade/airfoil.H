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
    Face

Description
    Face for 3rd order solver.

Author
    Minghao Yang

\*---------------------------------------------------------------------------*/

#ifndef airfoil_H
#define airfoil_H

#include "definitions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Airfoil Class Declaration
\*---------------------------------------------------------------------------*/

class NACA0012
{
public:
    static constexpr auto getLiftCoefficients()
    {
        return std::array<scalar, 23>
        {
         /* ---------  ---------  ---------  ---------  ---------  --------- */
            0.0000000, 0.1096268, 0.2185442, 0.3260494, 0.4312037, 0.5387641,
            0.665721,  0.7994149, 0.9245446, 1.007842,  1.094685,  1.183979,
            1.272707,  1.350418,  1.419637,  1.462286,  1.472421,  1.481117,
            1.444277,  1.33457,   1.204589,  1.124734,  1.067576
         /* ---------  ---------  ---------  ---------  ---------  ---------  */
        };
    }
    static constexpr auto getDragCoefficients()
    {
        return std::array<scalar, 23>
        {
         /* -----------  -----------  -----------  -----------  -----------  ----------- */
            0.005224612, 0.005296723, 0.005549942, 0.006014264, 0.006738895, 0.007691279,
            0.008764942, 0.00986979,  0.01095572,  0.01214015,  0.01341638,  0.01495878,
            0.01665783,  0.01906542,  0.02139354,  0.0249273,   0.03219501,  0.04226846,
            0.06086352,  0.09298569,  0.1316034,   0.1644066,   0.1961895
         /* ---------  ---------  ---------  ---------  ---------  ---------  --------- */
        };
    }
};

class NACA0014
{
public:
    static constexpr auto getLiftCoefficients()
    {
        return std::array<scalar, 29>
        {
         /* ---------  ---------  ---------  ---------  ---------  ---------  --------- */
            0.0000000, 0.1137715, 0.2267232, 0.3391929, 0.4506533, 0.5607037, 0.6689544,
            0.7776738, 0.8828297, 0.9946617, 1.121721,  1.244619,  1.321202,  1.410113,
            1.490664,  1.545651,  1.601115,  1.649106,  1.685385,  1.702667,  1.687639,
            1.623978,  1.505279,  1.348773,  1.225356,  1.133107,  1.08337,   1.061561,
            1.044248
         /* ---------  ---------  ---------  ---------  ---------  ---------  */
        };
    }
    static constexpr auto getDragCoefficients()
    {
        return std::array<scalar, 29>
        {
         /* -----------  -----------  -----------  -----------  -----------  ----------- */
            0.005383429, 0.00545056,  0.005602818, 0.005876952, 0.006284965, 0.006779944,
            0.00733932,  0.008054875, 0.008779963, 0.009670529, 0.01080122,  0.0119718,
            0.01318092,  0.01432266,  0.01603245,  0.01848413,  0.02177412,  0.02520666,
            0.03085292,  0.03935539,  0.05277494,  0.07395128,  0.1045583,   0.1439745,
            0.1820613,   0.2200939,   0.2523322,   0.2821259,   0.3122979
         /* ---------  ---------  ---------  ---------  ---------  ---------  --------- */
        };
    }
};

class NACA0015
{
public:
    static constexpr auto getLiftCoefficients()
    {
        return std::array<scalar, 29>
        {
         /* ---------  ---------  ---------  ---------  ---------  --------- */
            0.000000,  0.1204435, 0.2407259, 0.360358,  0.4794342, 0.5981435,
            0.7164582, 0.8340779, 0.9503337, 1.06396,   1.171626,  1.292637,
            1.41336,   1.487677,  1.556699,  1.622683,  1.668047,  1.681242,
            1.666333,  1.598837,  1.482079,  1.347875,  1.240259,  1.171396,
            1.130404,  1.115472,  1.117216,  1.124022,  1.111072
         /* ---------  ---------  ---------  ---------  ---------  ---------  */
        };
    }
    static constexpr auto getDragCoefficients()
    {
        return std::array<scalar, 29>
        {
         /* -----------  -----------  -----------  -----------  -----------  ----------- */
            0.005557139, 0.005612116, 0.005773551, 0.006028246, 0.006386258, 0.006814381,
            0.007310459, 0.007893265, 0.008605149, 0.009423933, 0.01042989,  0.0117386,
            0.01325721,  0.01481191,  0.01736391,  0.02076258,  0.02646984,  0.03626848,
            0.05040023,  0.07144987,  0.09975926,  0.1327155,   0.1652533,   0.1947299,
            0.22183,     0.2448496,   0.2646786,   0.2839075,   0.312397
         /* ---------  ---------  ---------  ---------  ---------  ---------  --------- */
        };
    }
};

class SC1095
{
public:
    static constexpr auto getLiftCoefficients()
    {
        return std::array<scalar, 181>
        {
            0.0, 0.25667, 0.5133, 0.77, 0.755,
            0.74, 0.725, 0.710, 0.695, 0.68,
            0.665, 0.65, 0.725, 0.8, 0.875,
            0.95, 0.9175, 0.886, 0.8525, 0.82,
            0.7875, 0.755, 0.7225, 0.69, 0.6575,
            0.625, 0.5925, 0.56, 0.5275, 0.495,
            0.4625, 0.43, 0.3975, 0.365, 0.3325,
            0.3, 0.2675, 0.235, 0.2025, 0.17,
            0.1375, 0.105, 0.0725, 0.04, 0.0075,
            -0.025, -0.0575, -0.09, -0.1225, -0.155,
            -0.1875, -0.22, -0.2525, -0.285, -0.3175,
            -0.35, -0.3825, -0.415, -0.4475, -0.48,
            -0.5125, -0.545, -0.5775, -0.61, -0.6425,
            -0.675, -0.7075, -0.74, -0.7725, -0.805,
            -0.8357, -0.87, -0.9025, -0.935, -0.9675,
            -1.0, -0.996, -0.992, -0.988, -0.984,
            -0.98, -0.976, -0.972, -1.07, -0.724,
            -0.37, -0.19, -0.39, -0.46, -0.19,
            0.03, 0.243, 0.46, 0.67, 0.89,
            1.1, 1.25, 1.1, 0.98, 0.9828,
            0.9856, 0.9884, 0.9912, 0.994, 0.997,
            1.0, 0.9675, 0.935, 0.9025, 0.87,
            0.8375, 0.805, 0.7725, 0.74, 0.7075,
            0.675, 0.6425, 0.61, 0.5775, 0.548,
            0.5125, 0.48, 0.4475, 0.415, 0.3825,
            0.35, 0.3175, 0.285, 0.2525, 0.22,
            0.1875, 0.155, 0.1225, 0.09, 0.0575,
            0.025, -0.0075, -0.04, -0.0725, -0.105,
            -0.1375, -0.17, -0.2025, -0.236, -0.2675,
            -0.3, -0.3325, -0.365, -0.3975, -0.43,
            -0.4625, -0.495, -0.5275, -0.56, -0.5925,
            -0.625, -0.6575, -0.69, -0.7225, -0.755,
            -0.7875, -0.82, -0.8525, -0.885, -0.9175,
            -0.95, -0.875, -0.8, -0.725, -0.65,
            -0.665, -0.68, -0.695, -0.71, -0.725,
            -0.74, -0.755, -0.77, -0.5133, -0.25667,
            0.0
        };
    }
    static constexpr auto getDragCoefficients()
    {
        return std::array<scalar, 181>
        {
            0.0, 0.367, 0.733, 1.1, 1.065,
            1.03, 0.995, 0.96, 0.92, 0.88,
            0.84, 0.8, 0.76, 0.72, 0.68,
            0.64, 0.6875, 0.735, 0.7825, 0.83,
            0.87875, 0.9275, 0.97625, 1.025, 1.07375,
            1.1225, 1.17125, 1.22, 1.26875, 1.3175,
            1.36625, 1.415, 1.46375, 1.5125, 1.56125,
            1.61, 1.66, 1.71, 1.76, 1.81,
            1.8575, 1.905, 1.9525, 2.0, 2.02625,
            2.0525, 2.07875, 2.105, 2.09875, 2.0925,
            2.08625, 2.08, 2.055, 2.03, 2.005,
            1.98, 1.92875, 1.8775, 1.82625, 1.775,
            1.70125, 1.6275, 1.55375, 1.48, 1.4075,
            1.335, 1.2625, 1.19, 1.1175, 1.045,
            0.9725, 0.9, 0.8325, 0.765, 0.6975,
            0.63, 0.562, 0.488, 0.417, 0.34,
            0.267, 0.195, 0.12, 0.045, 0.018,
            0.012, 0.008, 0.00775, 0.0075, 0.0075,
            0.0075, 0.008, 0.0085, 0.009, 0.011,
            0.017, 0.026, 0.145, 0.23, 0.293,
            0.345, 0.4, 0.455, 0.507, 0.56,
            0.63, 0.6975, 0.765, 0.8325, 0.9,
            0.9725, 1.045, 1.1175, 1.19, 1.2625,
            1.335, 1.4075, 1.48, 1.55375, 1.6275,
            1.70125, 1.775, 1.82625, 1.8775, 1.92875,
            1.98, 2.005, 2.03, 2.055, 2.08,
            2.08625, 2.0925, 2.09875, 2.105, 2.07875,
            2.0525, 2.02625, 2.0, 1.9525, 1.905,
            1.8575, 1.81, 1.76, 1.71, 1.66,
            1.61, 1.56125, 1.5125, 1.46375, 1.415,
            1.36625, 1.3175, 1.26875, 1.22, 1.17125,
            1.1225, 1.07375, 1.025, 0.97625, 0.9275,
            0.87875, 0.83, 0.7825, 0.735, 0.6875,
            0.64, 0.68, 0.72, 0.76, 0.8,
            0.84, 0.88, 0.92, 0.96, 0.996,
            1.03, 1.065, 1.1, 0.733, 0.367,
            0.0
        };
    }
};

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //