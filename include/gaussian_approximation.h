/*
 * Polar_Headers.h
 *
 *  Created on: Jun 19, 2020
 *      Author: murad
 */

#ifndef GAUSSIAN_APPROXIMATION_H_
#define GAUSSIAN_APPROXIMATION_H_

#include "general_headers.h"

double phi(double x);
double phi_residual(double x, double value);
double bisection(double value, double a, double b);
double phi_inverse(double value);
std::vector<double> gaussian_approximation(std::vector<double> init_mean);





#endif /* POLAR_HEADERS_H_ */
