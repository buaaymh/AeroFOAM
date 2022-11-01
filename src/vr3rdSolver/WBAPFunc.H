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

Foam::Mat5X3 Foam::d1VarWBAP
(
    const std::vector<Foam::Mat5X3>& list,
    scalar value
)
{
    const label n = list.size();
    Foam::Mat5X3 result = list[0];
    for (label i = 0; i != 5; i++)
    {
        for (label j = 0; j != 3; j++)
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

Foam::Mat5X6 Foam::d2VarWBAP
(
    const std::vector<Foam::Mat5X6>& list,
    scalar value
)
{
    const label n = list.size();
    Foam::Mat5X6 result = list[0];
    for (label i = 0; i != 5; i++)
    {
        for (label j = 0; j != 6; j++)
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
