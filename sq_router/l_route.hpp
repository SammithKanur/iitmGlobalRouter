#pragma once
#include "./sq_common.h"
#include "./sq_util.hpp"

namespace SQC = SQ_COMMON_SPACE;
namespace SQU = SQ_UTIL_SPACE;

namespace L_ROUTE_SPACE {
    class L_ROUTE : public SQC::PRR {
        public:
        L_ROUTE() {
            return;
        }
        static gcb_dt get_viad(gvia_dt z1, gvia_dt z2, int x1, int y1) {
            gcb_dt res = 0;
            for(int i = min(z1, z2); i <= max(z1, z2); i++) {
                res += cost_mat[UTIL::get_enc<CUBE_ENC>(i, x1, y1)];
            }
            res -= cost_mat[UTIL::get_enc<CUBE_ENC>(z2, x1, y1)];
        }
        static gcb_dt get_wired(gvia_dt z1, gcd_dt &x1, gcd_dt &y1, gvia_dt z2, gcd_dt &x2, gcd_dt &y2) {
            gcb_dt res = 0;
            #if SQR_CPU_ASSERT == 1
            SPDLOG_ASSERT(z1==z2, "");
            #endif
            res = (z1 & 1) ? abs(SQC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x1, y2)] - SQC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x1, y1)]) : abs(SQC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x2, y1)] - SQC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x1, y1)]);
            if(z1 & 1) {
                if(y1 > y2) {
                    res += cost_mat[UTIL::get_enc<CUBE_ENC>(z1, x1, y1)] - cost_mat[UTIL::get_enc<CUBE_ENC>(z1, x1, y2)];
                }
                y1=y2;
            }
            else {
                if(x1 > x2) {
                    res += cost_mat[UTIL::get_enc<CUBE_ENC>(z1, x1, y1)] - cost_mat[UTIL::get_enc<CUBE_ENC>(z1, x2, y1)];
                }
                x1=x2;
            }
        }
        static void add_via(vector<SQC::interval_t> &res, gvia_dt z1, gvia_dt z2, gcd_dt x1, gcd_dt y1) {
            genc_dt enc1 = UTIL::get_enc<CUBE_ENC>(z1, x1, y1), enc2 = UTIL::get_enc<CUBE_ENC>(z1, x1, y1);
            res.push_back({enc1 : enc1, enc2 : enc2});
        }
        static void add_wire(vector<SQC::interval_t> &res, gvia_dt z1, gcd_dt &x1, gcd_dt &y1, gvia_dt z2, gcd_dt &x2, gcd_dt &y2) {
            #if SQR_CPU_ASSERT == 1
            SPDLOG_ASSERT(z1==z2, "");
            #endif
            genc_dt enc1, enc2;
            if(z1 & 1) {
                enc1 = UTIL::get_enc<CUBE_ENC>(z1, x1, y1), enc2 = UTIL::get_enc<CUBE_ENC>(z1, x1, y2);
                y1 = y2;
            } else {
                enc1 = UTIL::get_enc<CUBE_ENC>(z1, x1, y1), enc2 = UTIL::get_enc<CUBE_ENC>(z1, x2, y1);
                x1 = x2;
            }
            res.push_back(SQC::interval_t{enc1: enc1, enc2: enc2});
        }
        static vector<SQC::interval_t> get_path(genc_dt src, genc_dt dst, int cpat) {
            gvia_dt z1, z2;
            gcd_dt x1, y1, x2, y2;
            UTIL::get_dec_all<CUBE_ENC>(src, z1, x1, y1);
            UTIL::get_dec_all<CUBE_ENC>(dst, z2, x2, y2);
            gvia_dt via1, via2;
            vector<SQC::interval_t> res;
            res.clear();
            SQU::get_zdec(cpat, via1, via2);
            if(via1 != 0) {
                if(z1 != via1) {
                    add_via(res, via1, z1, x1, y1);
                    z1 = via1;
                }
            }
            add_wire(res, z1, x1, y1, z1, x2, y2);
            if(via2 != 0) {
                via2 = (via1 & 1) ? via2 * 2 : via2 * 2 + 1;
                if(z1 != via2) {
                    add_via(res, z1, via2, x1, y1);
                    z1 = via2;
                }
            }
            add_wire(res, z1, x1, y1, z1, x2, y2);
            if(z1 != z2) {
                add_via(res, z1, z2, x1, y1);
                z1 = z2;
            }
            #if SQR_CPU_ASSERT == 1
            SPDLOG_ASSERT(z1 == z2 && x1 == x2 && y1 == y2, "[zxy] [{}, {}, {}] [{}, {}, {}]", z1, x1, y1, z2, x2, y2);
            #endif
            return res;
        }
        SQC::pr_route_t route(genc_dt src, genc_dt dst) const override {
            // spdlog::info("in lshape route log");
            if(src == dst)  {
                spdlog::info("src and dst are equal");
                return SQC::pr_route_t{edge_dist: 0};
            }                
            gvia_dt z1, z2;
            gcd_dt x1, y1, x2, y2;
            UTIL::get_dec_all<CUBE_ENC>(src, z1, x1, y1);
            UTIL::get_dec_all<CUBE_ENC>(dst, z2, x2, y2);
            gcb_dt edist = MAX_OFGPU_COST;
            int npat = -1;
            for(int c_pat = 0; c_pat < 50; c_pat++) {
                gvia_dt via1, via2;
                gcb_dt c_res = 0;
                SQU::get_zdec(c_pat, via1, via2);
                // // exit(0);
                if(via1 != 0) {
                    if(z1 != via1) {
                        c_res += get_viad(z1, via1, x1, y1);
                    }
                    z1 = via1;
                }
                c_res += get_wired(z1, x1, y1, z1, x2, y2);
                // if(via2 != 0) {
                //     via2 = (via1 & 1) ? via2 * 2 : via2 * 2 + 1;
                //     if(z1 != via2) {
                //         c_res += get_viad(z1, via2, x1, y1);
                //     }
                //     z1 = via2;
                // }
                // c_res += get_wired(z1, x1, y1, z1, x2, y2);
                // if(z1 != z2) {
                //     c_res += get_viad(z1, z2, x1, y1);
                //     z1 = z2;
                // }
                if(z1 != z2 || x1 != x2 || y1 != y2) {
                    continue;
                }
                #if SQR_CPU_ASSERT == 1
                SPDLOG_ASSERT(z1 == z2 && x1 == x2 && y1 == y2, "[zxy] [{}, {}, {}] [{}, {}, {}]", z1, x1, y1, z2, x2, y2);
                #endif
                c_res += cost_mat[UTIL::get_enc<CUBE_ENC>(z1, x1, y1)];
                if(edist > c_res) {
                    npat = c_pat;
                    edist = c_res;
                }
            }
            // #if SQR_CPU_ASSERT == 1
            // SPDLOG_ASSERT(npat != -1, "{}", npat);
            // #endif
            // vector<SQC::interval_t> crt = get_path(src, dst, npat);
            vector<SQC::interval_t> crt;
            vector<genc_dt> pth_pts = {src, dst};
            // spdlog::info("the total distance is {}", edist);
            return SQC::pr_route_t{edge_dist: edist, route: crt, path_pts: pth_pts};
        }
    };
}
