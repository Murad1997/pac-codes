/*
 * File: scl_decdr.h
 * Created: 2025-06-14
 * Author: murad
 *
 * Description:
 *
 */

#ifndef SCL_DECDR_H_
#define SCL_DECDR_H_
#include "conv_encoder.h"
#include "general_headers.h"

class SCL_DECODER {
public:
  SCL_DECODER();

  std::vector<uint8_t> scl_deocder(std::vector<double> rx_signal,
                                   std::vector<int> B);
  // void __MEM_FREE__();
  void reset_data();
  static SCL_DECODER *Instance(uint32_t list_size, std::vector<uint8_t> reg);

private:
  uint32_t L;
  uint32_t m;
  std::vector<uint8_t> activePath;
  std::vector<uint8_t> g;
  std::stack<uint32_t> inactivePathIndices;
  // std::vector<std::vector<double>> PM;
  // std::vector<std::vector<std::vector<double>>> llr;
  // std::vector<std::vector<std::vector<uint8_t>>> u_matrix;
  double ***llr;
  uint8_t ***u_matrix;
  double **PM;
  uint8_t **v_bits;

  static SCL_DECODER *_instance;

  std::vector<std::vector<uint8_t>> curstate;
  void killPath(uint32_t l);

  void duplicatePath(uint32_t l, uint32_t phi, double z_bit, double o_bit,
                     std::pair<uint8_t, std::vector<uint8_t>> mkz,
                     std::pair<uint8_t, std::vector<uint8_t>> mko);
  double pathmetric(double mu, double lamb, uint8_t u);

  void continueFrz_bit(uint32_t i);
  void conitnueUnFrz_bit(uint32_t i);

  int In_formation_(int phis, int lambda_index, int beta);

  void updatellrs_(double ***llr, uint8_t ***u_matrix, int l, int phi,
                   uint16_t lst);

  void updatepartialsum_(uint8_t ***u_matrix, int l, int phi, uint16_t lst);

  bool frozen_check_(int ind, std::vector<int> B);

  void allocate_mem();
};
#endif
