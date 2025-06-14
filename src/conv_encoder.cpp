
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include "conv_encoder.h"

std::pair<uint8_t, std::vector<uint8_t>> conv1bitEnc(uint8_t v, std::vector<uint8_t> currState, std::vector<uint8_t> g)
{
    uint8_t u = v & g[0];
    for (int j = 1; j < g.size(); ++j) {
        if (g[j] == 1) {
            u = u ^ currState[j - 1];
        }
    }

    std::vector<uint8_t> nextState(g.size() - 1, 0);
    nextState[0] = v;
    int n_i = 1;
    for (int i = 0; i < g.size() - 2; ++i) {
        nextState[n_i] = currState[i];
        n_i++;
    }
    return std::make_pair(u, nextState);
}

std::vector<uint8_t> convEncoder(std::vector<uint8_t> v, std::vector<uint8_t> g)
{
    std::vector<uint8_t> cState(g.size() - 1, 0);
    std::vector<uint8_t> u(v.size(), 0);
    std::pair<uint8_t, std::vector<uint8_t>> curpr;


    for (int i = 0; i < v.size(); ++i) {

        curpr = conv1bitEnc(v[i], cState, g);
        u[i] = curpr.first;
        cState = curpr.second;
    }

    return u;
}
