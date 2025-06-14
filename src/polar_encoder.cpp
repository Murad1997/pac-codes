
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

#include "polar_encoder.h"

int In_formation_enc(int phi, int lambda_index,int beta){
	return phi + (std::pow(2,lambda_index) * beta);
}


std::vector<uint8_t> rate_profiler(std::vector<uint8_t> u, std::vector<int> B)
{
    std::vector<uint8_t> v(N, 0);
    for (int i = 0; i < u.size(); ++i) {
        v[B[i]] = u[i];
    }

    return v;
}

std::vector<uint8_t> polarEncoder(std::vector<uint8_t> msg) {
    std::vector<uint8_t> code = msg;

    int n = std::log2(N);

    int inc = 2;
    for (int i = n - 1; i > -1; --i) {
        for (int j = 0; j < (int)std::pow(2, i); ++j) {

            std::vector<uint8_t>::const_iterator fst = code.begin() + inc * j;
            std::vector<uint8_t>::const_iterator lst = fst + inc;


            std::vector<uint8_t> a(fst, lst);
            int inds = (int)a.size() / 2;

            std::vector<uint8_t>::const_iterator first = a.begin() + inds;
            std::vector<uint8_t>::const_iterator second = a.end();
            std::vector<uint8_t> b(first, second);

            first = a.begin();
            second = a.end() - inds;
            std::vector<uint8_t> c(first, second);
            int up = 0;

            int st = inc * j;
            int end = b.size() + st;

            for (int t = st; t < end; ++t) {
                code[t] = b[up] ^ c[up];
                up++;
            }
            up = 0;
            for (int t = end; t < inc * (j + 1); ++t) {
                code[t] = b[up];
                up++;
            }


        }
        inc *= 2;
    }

    return code;
}


std::vector<uint8_t> polarEncoder_new(std::vector<uint8_t> msg){
	int n = std::log2(N);
	std::vector<std::vector<uint8_t>> P;
	P.resize(n+1,std::vector<uint8_t>(N,0));
	std::vector<uint8_t> code(N,0);
	for(int i=0; i<N; ++i){
		P[0][i] = msg[i];
	}
	for(int i=0; i<N; ++i){
		rpolarenc(P,n,i);
		code[i] = P[n][i];
	}
	return code;
}

void rpolarenc(std::vector<std::vector<uint8_t>> &P, int l, int phi){
	if(l==0){
		return ;
	}
	int si = phi/2;
	int l_ = l -1;
	if(phi % 2 ==0){
		rpolarenc(P,l_,si);
	}
	for(int beta=0; beta<std::pow(2,std::log2(N)-l); ++beta ){
		if(phi%2==0){

			P[l][In_formation_enc(phi,l,beta)] =P[l_][In_formation_enc(si,l_,2*beta)] ^ P[l_][In_formation_enc(si,l_,2*beta+1)];
		}
		else{
			P[l][In_formation_enc(phi,l,beta)] =P[l_][In_formation_enc(si,l_,2*beta + 1)];
		}

	}



}



