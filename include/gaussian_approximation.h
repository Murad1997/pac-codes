/*
 * File: gaussian_approximation.h
 * Created: 2025-06-14
 * Author: murad
 *
 * Description:
 *
 */

#ifndef GAUSSIAN_APPROXIMATION_H_
#define GAUSSIAN_APPROXIMATION_H_

#include "general_headers.h"
#include <iostream>

double phi(double x);
double phi_residual(double x, double value);
double bisection(double value, double a, double b);
double phi_inverse(double value);
std::vector<double> gaussian_approximation(std::vector<double> init_mean);

#endif /* POLAR_HEADERS_H_ */
