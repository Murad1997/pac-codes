

#include <armadillo>
#include "gaussian_approximation.h" // should be imported after the 
                                    // armadillo

#include "polar_decoders.h"
#include "polar_encoder.h"
#include "scl_decdr.h"
#include "tx_channel.h"

std::vector<int> RM_Design_rule();
void converting_to_imp_(unsigned int k, std::vector<uint8_t> *imp);
// START MAIN Here...
int main() {
  bool pac = false;
  std::cout << arma::arma_version::major << std::endl;

  std::vector<int> rm_ = RM_Design_rule();

  arma::vec v_rm = arma::conv_to<arma::vec>::from(rm_);
  arma::uvec indx_ = arma::sort_index(v_rm, "ascend");
  std::vector<int> B_unsorted(K, 0);

  // indx_={0,2,4,6,8,10,12,14};
  // indx_ = {0,1,2,3,4,5,6,7};

  // indx_ = {0,2,5,6,8,10,13,14};

  unsigned int rg_no = 66;
  for (int i = 0; i < K; ++i) {
    B_unsorted.at(i) = indx_(i);
    std::cout << B_unsorted.at(i) << ",";
  }
  std::cout << std::endl;
  // if(!pac){
  //     B_unsorted = {43,  84,  98,  71,  88,  45, 100,  51,  46, 104,  75, 53,
  //     112, 54,  77,  57,  83,  78,  58,  31,  85,  60,  99,  86,  89, 101,
  //     90,  47, 102,  92, 105, 106, 113, 108,  55, 114, 116,  79, 120,
  //     59,  61,  62,  87,  91, 103,  93,  94, 107, 109, 110, 115, 117,
  //    118, 121, 122, 124,  63,  95, 111, 119, 123, 125, 127, 126};

  // B_unsorted = { 97,  84,  98,  88,  45,  71, 100,  46,  51, 104, 112,  75,
  // 53,
  //     54,  77,  57,  78,  58,  83,  31,  60,  85,  86,  99,  89,  90,
  //    101,  92, 102,  47, 105, 106, 113, 108, 114, 116,  55, 120,  79,
  //     59,  61,  62,  87,  91,  93,  94, 103, 107, 109, 110, 115, 117,
  //    125, 124, 122, 121, 123,  63, 118, 111,  95, 126, 119, 127};

  // }

  int L = 32;
  float Rate = (float)K / N;
  double start_snr = 0.0;
  double end_snr = 5.2;
  double step_snr = 0.25;
  const int NBLK = 1e9;
  std::vector<double> snr_array;
  while (start_snr <= end_snr) {
    snr_array.push_back(start_snr);
    start_snr += step_snr;
  }
  for (int i = 0; i < snr_array.size(); ++i) {
    std::cout << snr_array[i] << ",";
  }
  std::cout << "\n";

  int n = std::log2(N);

  std::vector<uint8_t> indxs{0, 33, 43, 59}; //,6}; //4,9};//,13};//,33,43,59};
  std::vector<uint8_t> g(N, 0);              //{1,0,1,1,0,1,1};
  // std::vector<uint8_t> g{1,0,1,1,0,1,1};

  for (uint8_t i = 0; i < indxs.size(); ++i) {
    g.at(indxs.at(i)) = 1;
  }

  // for(uint8_t i=0; i<N; ++i){
  //     std::cout << (int)g.at(i) << ",";
  // }
  // return 0;
  // converting_to_imp_(rg_no,&g);

  double p_snr;
  int blk_errs;
  int bit_errs;
  double p_snr_lr;
  double sigma_;
  std::vector<uint8_t> msg;
  std::vector<uint8_t> v;
  std::vector<uint8_t> f_coded;
  std::vector<double> modulated;
  std::vector<double> rx_signal;
  std::vector<uint8_t> msg_cap;
  int errs;

  SCL_DECODER *scl = SCL_DECODER::Instance(L, g);
  for (int snr_indx = 0; snr_indx < snr_array.size(); ++snr_indx) {
    p_snr = snr_array[snr_indx];

    blk_errs = 0;
    bit_errs = 0;

    p_snr_lr = std::pow(10, p_snr / 10);
    sigma_ = std::sqrt(1 / (2 * Rate * p_snr_lr));
    // STARTING Simulation here ...
    for (int blk = 0; blk < NBLK; ++blk) {
      msg = gen_msg();
      v = rate_profiler(msg, B_unsorted);

      v = convEncoder(v, g);

      f_coded = polarEncoder_new(v);
      modulated = modulate(f_coded);
      rx_signal = transmitter(modulated, sigma_);
      for (int i = 0; i < rx_signal.size(); ++i) {
        rx_signal[i] = 2 * rx_signal[i] / (sigma_ * sigma_);
      }

      msg_cap = scl->scl_deocder(rx_signal, B_unsorted);
      //                 for(int i=0; i<N; ++i)
      // std::cout << (int)g.at(i) << ",";
      // std::cout << std::endl;
      // scl.__MEM_FREE__();

      errs = 0;
      // std::cout << "\nMsg_cap: " ;
      for (int i = 0; i < K; ++i) {
        // std::cout << (int)msg_cap[i] << ",";
        if (msg[i] != msg_cap[i]) {
          errs++;
        }
      }
      // std::cout << std::endl;
      if (errs > 0) {
        bit_errs += errs;
        blk_errs++;
      }

      if (blk_errs >= 110) {
        std::cout << "BLK: " << blk + 1 << ",SNR: " << p_snr
                  << ",FER: " << blk_errs << ",BER: " << bit_errs << "\n";
        std::ofstream outfile;
        outfile.open("sim_result.txt", std::fstream::app);

        outfile << "SCL: L: " << L << ",BLK:," << blk + 1 << ",SNR:," << p_snr
                << "(dB),FER:," << blk_errs << ",BER:," << bit_errs << ",N:,"
                << N << ",K:," << K << ",Rate: " << Rate << "\n";

        outfile.close();

        break;
      }
      if ((blk + 1) % 20 == 0) {
        std::cout << "SCL: L: " << L << ", BLK: " << blk + 1
                  << ",SNR: " << p_snr << "(dB),FER: " << blk_errs
                  << ",BER: " << bit_errs << ",N: " << N << ",K: " << K
                  << ",Rate: " << Rate << "\n";
      }
    }

  } // snrs _loop ends here

  return 0;
} // End Main here....

std::vector<int> RM_Design_rule() {
  int n = std::log2(N);
  std::vector<int> ones_(N, 0);
  for (int i = 0; i < N; ++i) {
    int k = 0;
    int j = i;
    int r = 0;
    while (k < n) {
      r = r + (j & 1);
      j = j >> 1;
      k++;
    }
    ones_[i] = -r;
  }
  return ones_;
}
void converting_to_imp_(unsigned int k, std::vector<uint8_t> *imp) {

  uint8_t indx = 0;
  imp->at(indx) = 1;
  while (indx < N - 1) {
    if (k & 1) {
      imp->at(N - 1 - indx) = 1;
    } else {
      imp->at(N - 1 - indx) = 0;
    }
    indx++;
    k >>= 1;
    // std::cout << "Indx: " << (N-1-indx) << std::endl;
  }
}
