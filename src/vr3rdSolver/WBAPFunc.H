template<int cols>
Eigen::Matrix<scalar, 5, cols> WBAP_L2_Limiter
(
    const std::vector<Eigen::Matrix<scalar, 5, cols>>& list,
    const scalar& value
)
{
    const label n = list.size();
    auto result = list[0];
    for (label i = 0; i != 5; i++)
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