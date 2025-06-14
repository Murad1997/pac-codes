
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include "gaussian_approximation.h"

double phi(double x) {
    double y;
    if (x < 10) {

        y = std::exp((-0.4527 * std::pow(x, 0.86)) + 0.0218);

    }
    else
    {
        y = std::sqrt(PI / x) * (1 - 10 / (7 * x)) * std::exp(-x / 4);

    }
    return y;
}

double phi_residual(double x, double value) {

    return phi(x) - value;

}

double bisection(double value, double a, double b) {
    double c = a;
    while ((b - a) >= 0.01)
    {
        c = (a + b) / 2;
        double half_residual = phi_residual(c, value);
        if (half_residual == 0.0)
            break;
        if (half_residual * phi_residual(a, value) < 0.0)
        {
            b = c;
        }
        else
        {
            a = c;
        }
    }
    return c;

}
double phi_inverse(double value) {
    return bisection(value, 0, 10000);

}
std::vector<double> gaussian_approximation(std::vector<double> init_mean) {


    int n = std::log2(N);
    std::vector<double > mean_llr(N, 0);
    std::vector<std::vector<double>> arr(n + 1, std::vector<double>(N, 0));

    for (int i = 0; i < N; ++i)
    {
        arr[0][i] = init_mean[i];
    }

    for (int j = 1; j <= n; ++j)
    {
        int u = std::pow(2, j);
        for (int i = 0; i < u; ++i)
        {


            if (i % 2 == 0)
            {

                int inds = (int)((i + 1) / 2);
                double T = arr[j - 1][inds];
                arr[j][i] = phi_inverse(1 - (1 - phi(T)) * (1 - phi(T)));
            }
            if (i % 2 == 1)
            {
                int inds = (int)((i) / 2);
                double T = arr[j - 1][inds];
                arr[j][i] = 2 * T;

            }
        }
    }
    for (int i = 0; i < N; ++i) {
        mean_llr[i] = arr[n][i];
    }


    return mean_llr;

}
