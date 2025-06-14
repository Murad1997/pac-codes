
#include "scl_decdr.h"

SCL_DECODER::SCL_DECODER() {}

SCL_DECODER *SCL_DECODER::_instance = 0;
SCL_DECODER *SCL_DECODER::Instance(uint32_t list_size,
                                   std::vector<uint8_t> reg) {
  if (_instance == 0) {
    _instance = new SCL_DECODER();
    _instance->L = list_size;
    _instance->m = (uint32_t)std::log2(N);
    _instance->g.assign(reg.begin(), reg.end());
    _instance->allocate_mem();
  }
  return _instance;
}

void SCL_DECODER::reset_data() {
  while (!inactivePathIndices.empty()) {
    inactivePathIndices.pop();
  }

  for (int l = 0; l < L; ++l) {
    for (int i = 0; i < N; ++i) {
      for (int m_ = 0; m_ < m; ++m_) {

        _instance->llr[m_][i][l] = 0;
        _instance->u_matrix[m_][i][l] = 0;
      }

      _instance->PM[l][i] = 0;
      _instance->v_bits[l][i] = 0;
    }
    for (int t = 0; t < g.size() - 1; ++t) {
      _instance->curstate.at(l).at(t) = 0;
    }

    _instance->activePath.at(l) = 0;
    _instance->inactivePathIndices.push(l);
  }
}

void SCL_DECODER::allocate_mem() {

  _instance->llr = new double **[m + 1];
  _instance->u_matrix = new uint8_t **[m + 1];
  _instance->PM = new double *[L];
  _instance->v_bits = new uint8_t *[L];

  for (uint16_t i = 0; i < L; ++i) {
    _instance->inactivePathIndices.push(i);
    _instance->activePath.push_back(0);
    _instance->curstate.push_back(std::vector<uint8_t>((int)(g.size() - 1), 0));
    _instance->PM[i] = new double[N];
    _instance->v_bits[i] = new uint8_t[N];
  }

  for (uint8_t j = 0; j <= m; ++j) {
    _instance->llr[j] = new double *[N];
    _instance->u_matrix[j] = new uint8_t *[N];
    for (uint16_t k = 0; k < N; ++k) {
      _instance->llr[j][k] = new double[L];
      _instance->u_matrix[j][k] = new uint8_t[L];
    }
  }
}

// void SCL_DECODER::__MEM_FREE__(){
//     for(uint8_t m_=0; m_<m; ++m_){
//         for(uint16_t i=0; i<N; ++i){
//             delete[] llr[m_][i];
//             delete[] u_matrix[m_][i];
//         }
//         delete[] llr[m_];
//         delete[] u_matrix[m_];
//     }
//     for(uint16_t l =0; l<L; ++l){
//         delete[] PM[l];
//         delete[] v_bits[l];
//     }
//     llr = NULL;
//     u_matrix=NULL;
//     PM = NULL;
//     v_bits = NULL;
// }

void SCL_DECODER::killPath(uint32_t l) {
  _instance->activePath[l] = 0;
  _instance->inactivePathIndices.push(l);
}
void SCL_DECODER::duplicatePath(uint32_t l, uint32_t phi, double z_bit,
                                double o_bit,
                                std::pair<uint8_t, std::vector<uint8_t>> mkz,
                                std::pair<uint8_t, std::vector<uint8_t>> mko) {
  uint32_t lp = _instance->inactivePathIndices.top();
  _instance->inactivePathIndices.pop();
  _instance->activePath[lp] = 1;
  for (uint16_t c = 0; c < N; c++) {
    for (uint16_t r = 0; r <= m; r++) {
      _instance->llr[r][c][lp] = _instance->llr[r][c][l];
      _instance->u_matrix[r][c][lp] = _instance->u_matrix[r][c][l];
    }
    v_bits[lp][c] = v_bits[l][c];
    PM[lp][c] = PM[l][c];
  }

  v_bits[l][phi] = 0;
  curstate.at(l) = mkz.second;
  u_matrix[m][phi][l] = mkz.first;
  PM[l][phi] = z_bit;

  v_bits[lp][phi] = 1;
  curstate.at(lp) = mko.second;
  u_matrix[m][phi][lp] = mko.first;
  PM[lp][phi] = o_bit;
}
double SCL_DECODER::pathmetric(double mu, double lamb, uint8_t u) {
  if (u == (uint8_t)(0.5 * (1 - SIGN(lamb)))) {
    return mu;
  } else {
    return mu + ABS(lamb);
  }
}
void SCL_DECODER::conitnueUnFrz_bit(uint32_t i) {
  SCL_DECODER *scl = SCL_DECODER::Instance(L, g);

  std::vector<double> probForks(2 * L, 0);
  std::vector<bool> contForks(2 * L, 0);
  std::vector<double> probabilities;
  std::vector<std::pair<uint8_t, std::vector<uint8_t>>> mkz;
  std::vector<std::pair<uint8_t, std::vector<uint8_t>>> mko;
  std::vector<uint16_t> pair_list_index(L, 0);

  int ci = 0;
  for (int l = 0; l < L; ++l) {
    if (activePath[l]) {

      mkz.push_back(conv1bitEnc(0, curstate.at(l), g));
      mko.push_back(conv1bitEnc(1, curstate.at(l), g));

      if (i == 0) {
        probForks[2 * l] = pathmetric(0, llr[m][i][l], mkz.at(ci).first);
        probForks[2 * l + 1] = pathmetric(0, llr[m][i][l], mko.at(ci).first);
      } else {
        probForks[2 * l] =
            pathmetric(PM[l][i - 1], llr[m][i][l], mkz.at(ci).first);

        probForks[2 * l + 1] =
            pathmetric(PM[l][i - 1], llr[m][i][l], mko.at(ci).first);
      }

      probabilities.push_back(-probForks.at(2 * l));
      probabilities.push_back(-probForks.at(2 * l + 1));
      pair_list_index.at(l) = ci;

      ci++;
    } else {
      probForks[2 * l] = 1e20;
      probForks[2 * l + 1] = 1e20;

      // probabilities.push_back(probForks.at(2*l));
      // probabilities.push_back(probForks.at(2*l + 1));
    }
  }

  int rho = MIN(2 * ci, (int)L);

  std::sort(probabilities.begin(), probabilities.end(),
            std::greater<double>()); // in descending order

  double thr = -probabilities.at(rho - 1); // changed rho-1

  int num_paths_continued = 0;
  for (int l = 0; l < 2 * L; ++l) {
    if (probForks[l] < thr) {
      contForks[l] = 1;
      num_paths_continued++;
    }
    if (num_paths_continued == rho) {
      break;
    }
  }

  if (num_paths_continued < rho) {
    for (int l = 0; l < 2 * L; ++l) {
      if (probForks[l] == thr) {
        contForks[l] = 1;
        num_paths_continued++;
      }
      if (num_paths_continued == rho) {
        break;
      }
    }
  }

  for (int l = 0; l < L; ++l) {
    if (activePath[l] == 0) {
      continue;
    }
    if (contForks[2 * l] == 0 && contForks[2 * l + 1] == 0) {
      killPath(l);
    }
  }

  for (int l = 0; l < L; ++l) {
    if (contForks[2 * l] == 0 && contForks[2 * l + 1] == 0) {
      continue;
    }

    if (contForks[2 * l] == 1 && contForks[2 * l + 1] == 1) {

      duplicatePath(l, i, probForks[2 * l], probForks[2 * l + 1],
                    mkz.at(pair_list_index.at(l)),
                    mko.at(pair_list_index.at(l)));

    } else {
      if (contForks[2 * l] == 1) {
        v_bits[l][i] = 0;

        curstate.at(l) = mkz.at(pair_list_index.at(l)).second;
        u_matrix[m][i][l] = mkz.at(pair_list_index.at(l)).first;
        PM[l][i] = probForks[2 * l];
      } else {
        v_bits[l][i] = 1;
        curstate.at(l) = mko.at(pair_list_index.at(l)).second;
        u_matrix[m][i][l] = mko.at(pair_list_index.at(l)).first;
        PM[l][i] = probForks[2 * l + 1];
      }
    }
  }
}

void SCL_DECODER::continueFrz_bit(uint32_t i) {
  for (uint16_t l = 0; l < L; ++l) {
    if (!activePath.at(l)) {
      continue;
    }
    // u_matrix[n][In_formation(phi, n, 0)] = 0;

    std::pair<uint8_t, std::vector<uint8_t>> mkp;
    mkp = conv1bitEnc(0, curstate.at(l), g);
    v_bits[l][i] = 0;
    curstate.at(l) = mkp.second;
    uint8_t v = mkp.first;
    u_matrix[m][In_formation_(i, m, 0)][l] = v;
    if (i == 0) {
      PM[l][i] = pathmetric(0, llr[m][i][l], v);
    } else {
      PM[l][i] = pathmetric(PM[l][i - 1], llr[m][i][l], v);
    }
  }
}

std::vector<uint8_t> SCL_DECODER::scl_deocder(std::vector<double> rx_signal,
                                              std::vector<int> B) {
  /// this->allocate_mem();

  reset_data();

  for (uint16_t l = 0; l < L; ++l) {
    for (uint16_t j = 0; j < N; ++j) {
      llr[0][j][l] = rx_signal.at(j);
    }
  }

  uint32_t l_c = inactivePathIndices.top();
  inactivePathIndices.pop();

  activePath.at(l_c) = 1;

  for (uint16_t i = 0; i < N; ++i) {

    for (uint16_t l = 0; l < L; ++l) {
      if (activePath.at(l)) {
        updatellrs_(llr, u_matrix, m, i, l);
      }
    }

    if (frozen_check_(i, B)) {
      // std::cout << "i: " << (int)i << std::endl;
      continueFrz_bit(i);
      // std::cout << "i: " << (int)i << " after " << std::endl;
    } else {
      conitnueUnFrz_bit(i);
    }
    if (i % 2 == 1) {
      for (uint16_t l = 0; l < L; ++l) {
        if (activePath.at(l)) {
          updatepartialsum_(u_matrix, m, i, l);
        }
      }
    }
  }

  double min_pm = 1e20;
  uint16_t min_indx = 0;
  // std::cout << "\n Path metrics: \n";
  for (uint16_t l = 0; l < L; ++l) {
    // for(int upt = 0; upt<N; ++upt){
    //     std::cout << PM[l][upt] << ",";
    // }
    // std::cout << std::endl;
    // std::cout << "MinPm: " << min_pm << std::endl;
    if (min_pm > PM[l][N - 1]) {
      min_indx = l;
      min_pm = PM[l][N - 1];
    }
  }
  // std::cout << std::endl;
  std::vector<uint8_t> msg_cap(K, 0);

  // std::cout << "MIN_indx: " << (int) min_indx << std::endl;
  for (uint16_t i = 0; i < K; ++i) {
    // msg_cap[i] = u_matrix[m][B[i]][min_indx];
    msg_cap[i] = v_bits[min_indx][B[i]];
  }

  // this->__MEM_FREE__();

  return msg_cap;
}

int SCL_DECODER::In_formation_(int phis, int lambda_index, int beta) {

  return phis + (std::pow(2, lambda_index) * beta);
}
void SCL_DECODER::updatellrs_(double ***llr, uint8_t ***u_matrix, int l,
                              int phi, uint16_t lst) {
  if (l == 0) {

    return;
  }
  int si = phi / 2;
  int l_ = l - 1;
  if (phi % 2 == 0) {
    updatellrs_(llr, u_matrix, l_, si, lst);
  }
  for (int beta = 0; beta < std::pow(2, std::log2(N) - l); ++beta) {

    double a = llr[l_][In_formation_(si, l_, 2 * beta)][lst];
    double b = llr[l_][In_formation_(si, l_, 2 * beta + 1)][lst];
    if (phi % 2 == 0) {
      // first case here...
      llr[l][In_formation_(phi, l, beta)][lst] =
          (double)SIGN(a) * (double)SIGN(b) * (double)MIN(ABS(a), ABS(b));
    } else {
      // second case here
      int u_p = (int)u_matrix[l][In_formation_(phi - 1, l, beta)][lst];
      llr[l][In_formation_(phi, l, beta)][lst] =
          ((1 - 2 * (double)u_p)) * a + b;
    }
  }
}
void SCL_DECODER::updatepartialsum_(uint8_t ***u_matrix, int l, int phi,
                                    uint16_t lst) {
  int si = phi / 2;
  int l_ = l - 1;

  for (int beta = 0; beta < std::pow(2, m - l); ++beta) {

    u_matrix[l_][In_formation_(si, l_, 2 * beta)][lst] =
        u_matrix[l][In_formation_(phi - 1, l, beta)][lst] ^
        u_matrix[l][In_formation_(phi, l, beta)][lst];
    u_matrix[l_][In_formation_(si, l_, 2 * beta + 1)][lst] =
        u_matrix[l][In_formation_(phi, l, beta)][lst];
  }

  if (si % 2 == 1) {
    updatepartialsum_(u_matrix, l_, si, lst);
  }
}

bool SCL_DECODER::frozen_check_(int ind, std::vector<int> B) {
  bool flag = true;

  for (int i = 0; i < B.size(); ++i) {

    if (ind == B.at(i)) {
      flag = false;
      break;
    }
  }

  return flag;
}
