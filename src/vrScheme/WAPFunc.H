vector Foam::vectorWAP
(
    const std::vector<vector>& list,
    scalar value
)
{   
    const label n = list.size();
    vector result = list[0];
    for (label i = 0; i != 3; ++i)
    {
        if (mag(list[0][i]) < 1e-6) continue;
        scalar temp3 = value;
        scalar temp4 = value;
        for (label j = 1; j != n; ++j)
        {
            const scalar theta = list[j][i] / list[0][i];
            if (theta <= 1e-6)
            {
                temp3 = 0;
                break;
            }
            temp3 += 1.0 / Foam::pow3(theta);
            temp4 += 1.0 / Foam::pow4(theta);
        }
        result[i] *= temp3 / temp4;
    }
    return result;
}

symmTensor Foam::symmTensorWAP
(
    const std::vector<symmTensor>& list,
    scalar value
)
{
    const label n = list.size();
    symmTensor result = list[0];
    for (label i = 0; i != 6; ++i)
    {
        if (mag(list[0][i]) < 1e-6) continue;
        scalar temp3 = value;
        scalar temp4 = value;
        for (label j = 1; j != n; ++j)
        {
            const scalar theta = list[j][i] / list[0][i];
            if (theta <= 1e-6)
            {
                temp3 = 0;
                break;
            }
            temp3 += 1.0 / Foam::pow3(theta);
            temp4 += 1.0 / Foam::pow4(theta);
        }
        result[i] *= temp3 / temp4;
    }
    return result;
}

tensor Foam::tensorWAP
(
    const std::vector<tensor>& list,
    scalar value
)
{
    const label n = list.size();
    tensor result = list[0];
    for (label i = 0; i != 9; ++i)
    {
        if (mag(list[0][i]) < 1e-6) continue;
        scalar temp3 = value;
        scalar temp4 = value;
        for (label j = 1; j != n; ++j)
        {
            const scalar theta = list[j][i] / list[0][i];
            if (theta <= 1e-6)
            {
                temp3 = 0;
                break;
            }
            temp3 += 1.0 / Foam::pow3(theta);
            temp4 += 1.0 / Foam::pow4(theta);
        }
        result[i] *= temp3 / temp4;
    }
    return result;
}