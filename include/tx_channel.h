/*
 * File: tx_channel.h
 * Created: 2025-06-14
 * Author: murad
 *
 * Description:
 *
 */

#ifndef TX_CHANNEL_H_
#define TX_CHANNEL_H_

#include "general_headers.h"

std::vector<uint8_t> gen_msg();
std::vector<double> modulate(std::vector<uint8_t> coded_bits);
std::vector<double> transmitter(std::vector<double> modulated, double sigma_);

#endif /* TX_CHANNEL_H_ */
