/*
 * File: polar_decoders.h
 * Created: 2025-06-14
 * Author: murad
 *
 * Description:
 *
 */

#ifndef POLAR_DECODERS_H_
#define POLAR_DECODERS_H_

#define CAP 1300000

#include "conv_encoder.h"
#include "general_headers.h"

int In_formation(int phis, int lambda_index, int beta);
void updatellrs(std::vector<std::vector<double>> &llr,
                std::vector<std::vector<uint8_t>> &u_matrix, int l, int phi);
void updatepartialsum(std::vector<std::vector<uint8_t>> &u_matrix, int l,
                      int phi);

bool frozen_check(int ind, std::vector<int> B);
std::vector<int> reverse_order();

std::vector<uint8_t> sc_decoding(std::vector<double> rx_signal,
                                 std::vector<int> B);
// std::vector<uint8_t> scl_decoding(std::vector<double> rx_signal,
// std::vector<int> B);

std::pair<int, std::vector<uint8_t>> pac_decoding(std::vector<double> rx_signal,
                                                  std::vector<int> Bvec_us,
                                                  std::vector<double> pe,
                                                  std::vector<uint8_t> g);

double update_Th(double T, double delta, double tau);
std::vector<double> BACK(std::vector<double> betas, int j, double T,
                         std::vector<uint8_t> kappa, double delta);
double actual_metric_cal(double lru, uint8_t ucap);
double cal_PM(int i, std::vector<std::vector<double>> llr,
              std::vector<std::vector<uint8_t>> u_matrix);
double cal_pe(int i, std::vector<double> pe);

#endif /* POLAR_DECODERS_H_ */
