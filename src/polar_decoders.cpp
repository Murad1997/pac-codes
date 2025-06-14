
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

#include "polar_decoders.h"






std::pair<int ,std::vector<uint8_t>> pac_decoding(std::vector<double>rx_signal, std::vector<int> Bvec_us, std::vector<double> pe,std::vector<uint8_t> g) {
    std::vector<int>Bvec;
    Bvec.assign(Bvec_us.begin(),Bvec_us.end());
    std::sort(Bvec.begin(), Bvec.end());

    std::vector<uint8_t> msg_cap(K, 0);
    int n = std::log2(N);
    std::vector<std::vector<double>> llr;
    llr.resize(n + 1, std::vector<double>(N, 0));
    std::vector<std::vector<uint8_t>> u_matrix;
    u_matrix.resize(n + 1, std::vector<uint8_t>(N, 0));
//    std::vector<uint8_t> msg_cap(K, 0);
    for (int i = 0; i < rx_signal.size(); ++i) {
        llr[0][i] = rx_signal[i];
    }
    int i = 0;
    int j = -1;
    double T = 0;
    double B = -1;
    double delta = 1;
    std::vector<double> betas(K, 0);
    std::vector<uint8_t> kappa(K, 0);
    std::vector<double> mus(N, 0);
    std::vector<std::vector<uint8_t>> currState;
    currState.resize(K + 1, std::vector<uint8_t>(g.size() - 1, 0));
    std::vector<uint8_t> cState(g.size() - 1, 0);
    std::vector<uint8_t> vis(N, 0);
    double mu = 0;
    double expected_metric = cal_pe(N, pe);

    int v_nodes = 0;
    while (i != N) {
    	v_nodes++;


        updatellrs(llr, u_matrix, n, i);
        double cur_llr = llr[n][i];



        if (frozen_check(i, Bvec)) {
            std::pair<uint8_t, std::vector<uint8_t>> curpr;
            curpr = conv1bitEnc(0, cState, g);

            uint8_t ucap0 = curpr.first;
            u_matrix[n][i] = ucap0;
            cState = curpr.second;
            if (i % 2 == 1) {
                updatepartialsum(u_matrix, n, i);
            }


            if (i == 0) {
                mus[i] = actual_metric_cal(cur_llr, ucap0) + expected_metric;
            }
            else {
                mus[i] = mus[i - 1] + actual_metric_cal(cur_llr, ucap0);
            }
            mus[i] = mus[i] - cal_pe(i + 1, pe);
            vis[i] = ucap0;
            i++;
        }
        else {
            double exp_metric = cal_pe(i + 1, pe);
            std::pair<uint8_t, std::vector<uint8_t>> curpr0;
            std::pair<uint8_t, std::vector<uint8_t>> curpr1;
            curpr0 = conv1bitEnc(0, cState,g);
            curpr1 = conv1bitEnc(1, cState,g);
            double mu_0 = cal_PM(i - 1,llr,u_matrix) + actual_metric_cal(cur_llr, curpr0.first) - exp_metric;
            double mu_1 = cal_PM(i - 1, llr, u_matrix) + actual_metric_cal(cur_llr, curpr1.first) - exp_metric;
            std::vector<uint8_t> ucaps(2, 0);
            std::vector<std::vector<uint8_t>> cms_01;
            cms_01.resize(2, std::vector<uint8_t>(g.size() - 1, 0));
            cms_01[0] = curpr0.second;
            cms_01[1] = curpr1.second;

            ucaps[0] = curpr0.first;
            ucaps[1] = curpr1.first;

            double mu_max = MAX(mu_0, mu_1);
            double mu_min = MIN(mu_0, mu_1);

            int v_max = 0;
            int v_min = 0;
            if (mu_max == mu_0) {
                v_max = 0;
            }
            else {
                v_max = 1;
            }
            if (mu_min == mu_0) {
                v_min = 0;
            }
            else {
                v_min = 1;
            }

            if (mu_max > T) {
                if (B == -1) {
                    u_matrix[n][i] = ucaps[v_max];
                    if (i % 2 == 1) {
                        updatepartialsum(u_matrix, n, i);
                    }

                    //updatepartialsum(u_matrix, n, i);

                    kappa[j + 1] = 0;
                    betas[j + 1] = mu_max;
                    mus[i] = mu_max;
                    vis[i] = (uint8_t)v_max;
                    currState[j + 1] = cState;
                    cState = cms_01[v_max];
                    if (j == -1) {
                        mu = 0;
                    }
                    else if (j != -1) {
                        mu = betas[j];
                    }
                    if (mu < T + delta) {
                        T = update_Th(T, delta, betas[j + 1]);
                    }

                    j++;
                    i++;

                }
                else if (B == 1) {
                    if (mu_min > T) {
                        u_matrix[n][i] = ucaps[v_min];
                        vis[i] = (uint8_t)v_min;
                        if (i % 2 == 1) {
                            updatepartialsum(u_matrix, n, i);
                        }
                       // updatepartialsum(u_matrix, n, i);
                        betas[j + 1] = mu_min;
                        kappa[j + 1] = 1;
                        mus[i] = mu_min;
                        currState[j + 1] = cState;
                        cState = cms_01[v_min];
                        j++;
                        i++;
                        B = -1;
                    }
                    else if (mu_min <= T) {
                        if (j == -1) {
                            T -= delta;
                            B = -1;
                        }
                        else if (j != -1) {
                            int jp = j;
                            std::vector<uint8_t> csp = cState;
                            //BACK is to be called here ... with return arugment T,j,B
                           // std::cout << "back0_before:" << j << "\n";

                            std::vector<double> r_values = BACK(betas, j, T, kappa, delta);
                            T = r_values[0];
                            j = r_values[1];
                            B = r_values[2];
                            //std::cout << "back0_after:" << j << "\n";
                            //break;

                            if (jp == j) {
                                cState = csp;
                            }
                            else {
                                cState = currState[j + 1];
                            }
                            i = Bvec[j + 1];
                        }
                    }
                }
            }
            else if (mu_max <= T) {
                if (j == -1) {
                    //here
                    T -= delta;

                }
                else if (j != -1) {
                    // here
                    int jp = j;
                    std::vector<uint8_t> csp = cState;
                    // BACK is to be called here with return argurment in order T,j,B

                   // for (int i = 0; i < cState.size(); ++i) {
                    //    std::cout << (int)cState[i] << ",";
                   // }
                    //std::cout << "\n";

                    //std::cout << "back1_before:" << T << ","<< j << ","<< B << "\n";
                    std::vector<double> r_values = BACK(betas, j, T, kappa, delta);
                    T = r_values[0];
                    j = r_values[1];
                    B = r_values[2];
                   // std::cout << "back1_after:" << T << "," << j << "," << B << "\n";

                    //break;

                    i = Bvec[j + 1];
                    if (j == jp) {
                        cState = csp;
                    }
                    else {
                        cState = currState[j + 1];
                    }

                }
            }



        }

    }



    for (int i = 0; i < K; ++i) {
        msg_cap[i] = vis[Bvec_us[i]];

    }

    std::pair<int,std::vector<uint8_t>> r_pair = std::make_pair(v_nodes,msg_cap);
    return r_pair;
}

double cal_pe(int i, std::vector<double> pe) {
    double ac = 0;
    for (int j = 0; j < i; ++j) {
        ac += std::log(1 - pe[j]);
    }
    return ac;
}
double cal_PM(int i, std::vector<std::vector<double>> llr, std::vector<std::vector<uint8_t>> u_matrix) {
    double ac_0 = 0;
    int n = std::log2(N);
    for (int j = 0; j <= i; ++j) {
        ac_0 += actual_metric_cal(llr[n][j], u_matrix[n][j]);
    }
    return ac_0;
}


double actual_metric_cal(double lru, uint8_t ucap) {
    double a0 = -1 * std::log(1 + std::exp(-(1 - 2 * (double)ucap) * lru));
    return a0;
}

std::vector<double> BACK(std::vector<double> betas, int j, double T, std::vector<uint8_t> kappa, double delta) {
    double mu;
    std::vector<double> r_values(3, 0);
    while (true) {
        if (j == -1) {
            mu = -1e99;
        }
        if (j == 0) {
            mu = 0;
        }
        else if (j >= 1) {
            mu = betas[j - 1];
        }
        if (mu >= T) {
            j = j - 1;
            if (kappa[j + 1] == 0) {
                int B = 1;
                r_values[0] = T;
                r_values[1] = j;
                r_values[2] = B;
                return r_values;
                // needs to return T,j,B here
            }
        }
        else if (mu < T) {
            T = T - delta;
            int B = -1;
            r_values[0] = T;
            r_values[1] = j;
            r_values[2] = B;
            return r_values;
            // needs to return T,j,B here also...
        }

    }

}



double update_Th(double T, double delta, double tau) {
    while (T + delta <= tau) {
        T += delta;
    }
    return T;
}

std::vector<uint8_t> sc_decoding(std::vector<double> rx_signal, std::vector<int> B){
    int n = std::log2(N);
    std::vector<std::vector<double>> llr;
    llr.resize(n + 1, std::vector<double>(N, 0));
    std::vector<std::vector<uint8_t>> u_matrix;
    u_matrix.resize(n + 1, std::vector<uint8_t>(N, 0));
    std::vector<uint8_t> msg_cap(K, 0);
    for (int i = 0; i < rx_signal.size(); ++i) {
        llr[0][i] = rx_signal[i];
    }

    for (int phi = 0; phi < N; ++phi) {
        updatellrs(llr, u_matrix, n, phi);
        if (frozen_check(phi, B)) {
            //first part
            u_matrix[n][In_formation(phi, n, 0)] = 0;
        }
        else {
            // second part
            if (llr[n][In_formation(phi, n, 0)] >= 0) {
                u_matrix[n][In_formation(phi, n, 0)] = 0;
            }
            else {
                u_matrix[n][In_formation(phi, n, 0)] = 1;
            }
        }
        if (phi % 2 == 1) {
            updatepartialsum(u_matrix, n, phi);
        }

    }

    for (int i = 0; i < K; ++i) {
        msg_cap[i] = u_matrix[n][B[i]];
    }

    return msg_cap;
}

std::vector<int> reverse_order() {
    int n = std::log2(N);
    std::vector<int> rev_order(N,0);
    for (int i = 0; i < N; ++i) {
        int j = i;
        int r = 0;
        int k = 0;
        while (k < n) {
            r = (r << 1) | (j & 1);
            j = j >> 1;
            k++;
        }
        rev_order[i] = r;
    }
    return rev_order;
}

bool frozen_check(int ind, std::vector<int> B) {
    bool flag = true;
    for (int i = 0; i < B.size(); ++i) {
        if (ind == B[i]) {
            flag = false;
            break;
        }

    }
    return flag;
}







int In_formation(int phis, int lambda_index, int beta) {

    return phis + (std::pow(2, lambda_index) * beta);

}
void updatellrs(std::vector<std::vector<double>> &llr, std::vector<std::vector<uint8_t>> &u_matrix, int l, int phi) {
    if (l == 0) {
        return;
    }
    int si = phi / 2;
    int l_ = l - 1;
    if (phi % 2 == 0) {
        updatellrs(llr, u_matrix, l_, si);
    }
    for (int beta = 0; beta < std::pow(2, std::log2(N) - l); ++beta) {
        double a = llr[l_][In_formation(si, l_,2*beta)];
        double b = llr[l_][In_formation(si, l_, 2 * beta + 1)];
        if (phi % 2 == 0) {
            // first case here...
            llr[l][In_formation(phi, l, beta)] = (double)SIGN(a) * (double)SIGN(b) * (double) MIN(ABS(a), ABS(b));
        }
        else {
            // second case here
            int u_p = (int)u_matrix[l][In_formation(phi - 1, l, beta)];
            llr[l][In_formation(phi, l, beta)] = ((1 - 2 * (double)u_p)) * a + b;
        }
    }


}
void updatepartialsum(std::vector<std::vector<uint8_t>> &u_matrix, int l, int phi) {
    int si = phi / 2;
    int l_ = l - 1;

    for (int beta = 0; beta < std::pow(2, std::log2(N) - l); ++beta) {


        u_matrix[l_][In_formation(si, l_, 2 * beta)] = u_matrix[l][In_formation(phi - 1, l, beta)] ^ u_matrix[l][In_formation(phi, l, beta)];
        u_matrix[l_][In_formation(si, l_, 2 * beta + 1)] = u_matrix[l][In_formation(phi, l, beta)];
    }

    if (si % 2 == 1) {
        updatepartialsum(u_matrix, l_, si);
    }

}

