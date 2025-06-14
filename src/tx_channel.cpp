
#include "tx_channel.h"

// message generation ...
std::vector<uint8_t> gen_msg() {
  std::random_device rd{};
  std::mt19937 generator{rd()};
  std::bernoulli_distribution distribution(0.5);
  std::vector<uint8_t> msg;
  for (int i = 0; i < K; ++i) {
    msg.push_back(distribution(generator));
  }

  return msg;
}

// BPSK modulation function....
std::vector<double> modulate(std::vector<uint8_t> coded_bits) {
  std::vector<double> modulated(coded_bits.size(), 0);
  for (int i = 0; i < coded_bits.size(); ++i)
    modulated[i] = (1 - 2 * (double)coded_bits[i]);

  return modulated;
}

// Transmitter function....
std::vector<double> transmitter(std::vector<double> modulated, double sigma_) {
  std::random_device rd{};
  std::mt19937 generator{rd()};
  std::normal_distribution<double> distribution(0, sigma_);
  double u;
  for (int i = 0; i < modulated.size(); ++i) {
    modulated[i] = modulated[i] + distribution(generator);
  }
  return modulated;
}
