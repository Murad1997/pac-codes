/*
 * File: polar_encoder.h
 * Created: 2025-06-14
 * Author: murad
 *
 * Description:
 *
 */

#ifndef POLAR_ENCODER_H_
#define POLAR_ENCODER_H_
#include "general_headers.h"

std::vector<uint8_t> polarEncoder(std::vector<uint8_t>);
std::vector<uint8_t> rate_profiler(std::vector<uint8_t>, std::vector<int>);

int In_formation_enc(int phi, int lambda_index, int beta);

std::vector<uint8_t> polarEncoder_new(std::vector<uint8_t>);

void rpolarenc(std::vector<std::vector<uint8_t>> &P, int l, int phi);

#endif /* POLAR_ENCODER_H_ */
