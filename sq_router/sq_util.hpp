#pragma once

#include "./sq_common.h"

namespace SQC = SQ_COMMON_SPACE;

namespace SQ_UTIL_SPACE {
    inline void get_zdec(int val, gvia_dt &z1, gvia_dt &z2) {
        #if SQR_CPU_ASSERT == 1
        SPDLOG_ASSERT(val >= 0 && val < Z1Z2, "{} {}", val, Z1Z2);
        #endif
        z1 = val/Z2; z2 = val % Z2;
    }
    void build_cost_mat() {
        for(int z = 0; z < Z; z++) {
            for(int x = 0; x < X; x++) {
                for(int y = 0; y < Y; y++) {
                    genc_dt cenc = UTIL::get_enc<CUBE_ENC>(z, x, y);
                    cost_mat[cenc] = UTIL::get_cpu_of_cost(cube_data[cenc]);
                }
            }
        }
    }
    void get_cost_mat_pfx() {
        int enc_id = 0;
        for(gvia_dt z = 0; z < Z; z++) {
            if(z & 1) {
                for(gcd_dt x = 0; x < X; x++) {
                    gcb_dt pfx = 0, prv = 0;
                    for(gcd_dt y = 0; y < Y; y++) {
                        enc_id = UTIL::get_enc<CUBE_ENC>(z, x, y);
                        prv=cost_mat[enc_id];
                        SQC::pth_mat[enc_id] = pfx;
                        pfx += prv;
                    }
                }
            } else {
                for(gcd_dt y = 0; y < Y; y++) {
                    gcb_dt pfx = 0, prv = 0;
                    for(gcd_dt x = 0; x < X; x++) {
                        enc_id = UTIL::get_enc<CUBE_ENC>(z, x, y);
                        prv=cost_mat[enc_id];
                        SQC::pth_mat[enc_id] = pfx;
                        pfx += prv;
                    }
                }
            }
        }
    }
}
