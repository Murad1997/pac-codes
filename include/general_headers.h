/*
 * File: general_headers.h
 * Created: 2025-06-14
 * Author: murad
 *
 * Description:
 *
 */

#ifndef GENERAL_HEADERS_H_
#define GENERAL_HEADERS_H_

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <stack>
#include <vector>
#include "conv_encoder.h"

#define MAX(x, y) (x - y > 0 ? x : y)
#define MIN(x, y) (x - y < 0 ? x : y)
#define ABS(x) (x < 0 ? -x : x)
#define SIGN(x) (x < 0 ? -1 : 1)

#define PI 3.141592653589793238
#define N 64
#define K 32

#endif /* GENERAL_HEADERS_H_ */
