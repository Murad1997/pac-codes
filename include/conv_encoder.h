/*
 * conv_encoder.h
 *
 *  Created on: Jun 19, 2020
 *      Author: murad
 */

#ifndef CONV_ENCODER_H_
#define CONV_ENCODER_H_
#include "general_headers.h"

std::pair<uint8_t, std::vector<uint8_t>> conv1bitEnc(uint8_t v, std::vector<uint8_t> currState, std::vector<uint8_t> g);
std::vector<uint8_t> convEncoder(std::vector<uint8_t>, std::vector<uint8_t>);




#endif /* CONV_ENCODER_H_ */
