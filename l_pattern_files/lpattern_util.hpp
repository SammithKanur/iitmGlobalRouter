#pragma once
#include "../common.h"
#include "./lpattern_common.h"
#include "../util.hpp"
#include "../mst_files/mst_common.h"

using namespace std;
using namespace COMMON_DATA;
namespace LPC = LPAT_COMMON;
namespace MSC = MST_COMMON;

namespace LPAT_UTIL {
    void get_djkst_path(in_t src, in_t dst) {
        in_t cval = dst;
        string res = "";
        while(cval != -1) {
            gcd_dt x, y;
            gvia_dt z;
            UTIL::get_dec_all<CUBE_ENC>(cval, z, x, y);
            res += "(" + to_string(z) + "," + to_string(x) + "," + to_string(y) +  "," + to_string(wcost_cube[cval]) + "),";
            SPDLOG_ASSERT(cval != LPC::lpat_par[cval], "");
            cval = LPC::lpat_par[cval];
        }
        spdlog::info("djkst lpattern path is {}", res);
    }
    int gdir(in_t z1, in_t x1, in_t y1, in_t z2, in_t x2, in_t y2) {
        int cng = 0, res = -1;
        if(z1 != z2) {
            cng++;
            res = ZDIR;
        }
        if(x1 != x2) {
            cng++;
            res = XDIR;
        }
        if(y1 != y2) {
            cng++;
            res = YDIR;
        }
        SPDLOG_ASSERT(cng == 1, "");
        return res;
    }
    inline in_t get_zenc(in_t &z1, in_t &z2) {
        SPDLOG_ASSERT(z1 < Z && z1 >= 0 && z2 < CEIL_DIV(Z, 2) && z2 > 0, "");
        return z1*Z2 + z2;
    }    
    inline void get_zdec(gvia_dt &val, gvia_dt &z1, gvia_dt &z2) {
        SPDLOG_ASSERT(val < Z1*Z2 && val >= 0, "");
        z1 = val/Z2; z2 = val % Z2;
    }
    __inline__ __device__ void get_gpu_zdec(gvia_dt &val, gvia_dt &z1, gvia_dt &z2) {
        z1 = val / gvia_dt(DZ2); z2 = val % gvia_dt(DZ2);
    }
    __inline__ __device__ void add_gpu_vias(gvia_dt zs, gvia_dt z2, gcd_dt x, gcd_dt y, gcb_dt &via_cost, gcb_dt *d_wcost_cube) {
        for(gvia_dt i = zs; i < z2; i++) {
            via_cost += d_vcost_cube[UTIL::get_gpu_enc<VIA_ENC>(i, x, y)];
        }
    }
    inline void add_vias(gvia_dt zs, gvia_dt z2, gcd_dt x, gcd_dt y, gcb_dt &via_cost) {
        for(gvia_dt i = zs; i < z2; i++) {
            via_cost += vcost_cube[UTIL::get_enc<VIA_ENC>(i, x, y)];
        }
    }
    __inline__ __device__ void add_gpu_wires(gcb_dt &wcost, gvia_dt &zs, gcd_dt &xs, gcd_dt &ys, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, gcb_dt *d_wcost_cube) {
        wcost += (zs & 1) ? abs(d_wcost_cube[UTIL::get_gpu_enc<WIRE_ENC>(z2, xs, y2)] - d_wcost_cube[UTIL::get_gpu_enc<WIRE_ENC>(zs, xs, ys)]) : abs(d_wcost_cube[UTIL::get_gpu_enc<WIRE_ENC>(z2, x2, ys)] - d_wcost_cube[UTIL::get_gpu_enc<WIRE_ENC>(zs, xs, ys)]);
        if(zs & 1) {ys=y2;}
        else {xs=x2;}        
    }
    inline void add_wires(gcb_dt &wcost, gvia_dt &zs, gcd_dt &xs, gcd_dt &ys, gvia_dt z2, gcd_dt x2, gcd_dt y2) {
        SPDLOG_ASSERT(zs==z2, "");
        wcost += (zs & 1) ? abs(wcost_cube[UTIL::get_enc<WIRE_ENC>(z2, xs, y2)] - wcost_cube[UTIL::get_enc<WIRE_ENC>(zs, xs, ys)]) : abs(wcost_cube[UTIL::get_enc<WIRE_ENC>(z2, x2, ys)] - wcost_cube[UTIL::get_enc<WIRE_ENC>(zs, xs, ys)]);
        if(zs & 1) {ys=y2;}
        else {xs=x2;}        
    }
    void print_lpat_path(genc_dt enc1, genc_dt enc2, gvia_dt pth) {
        gcd_dt  x1, y1, x2, y2, xs, ys;
        gvia_dt via1, via2, z1, z2, zs;
        get_zdec(pth, via1, via2);
        UTIL::get_dec_all<CUBE_ENC>(enc1, z1, x1, y1);
        UTIL::get_dec_all<CUBE_ENC>(enc2, z2, x2, y2);
        spdlog::info("checking lpat path for [z.x.y] ({},{},{}) ({},{},{})", z1, x1, y1, z2, x2, y2);
        xs = x1; ys = y1; zs = z1;
        gcd_dt res = 0;
        // if(z1 != via1) {
            
        // }
        // string ptres = "";
        // if(z1 != via1) {
        //     res += cost_cube[UTIL::get_enc(zs, xs, ys)];
        // }
        // for(int i = min(via1, z1) + 1; i < max(via1, z1); i++) {
        //     res += cost_cube[UTIL::get_enc(i, xs, ys)];
        //     ptres += "(" + to_string(i) + "," + to_string(xs) + "," + to_string(ys) + "),";
        // }
        // zs = via1;
        // //at via1 x1 y1
        // if(zs & 1) {
        //     for(int y = min(y1, y2); y <= max(y1, y2); y++) {
        //         res += cost_cube[UTIL::get_enc(zs, xs, y)];
        //         ptres += "(" + to_string(zs) + "," + to_string(xs) + "," + to_string(y) + "),";
        //     }
        //     ys = y2;
        // } else {
        //     for(int x = min(x1, x2); x <= max(x1, x2); x++) {
        //         res += cost_cube[UTIL::get_enc(zs, x, ys)];
        //         ptres += "(" + to_string(zs) + "," + to_string(x) + "," + to_string(ys) + "),";
        //     }
        //     xs = x2;
        // }
        // //add via2 now
        // via2 = (zs & 1) ? via2 * 2 : via2 * 2 + 1;
        // if(via2 > 0) {
        //     for(int i = min(via1, via2) + 1; i < max(via1, via2); i++) {
        //         res += cost_cube[UTIL::get_enc(i, xs, ys)];
        //         ptres += "(" + to_string(i) + "," + to_string(xs) + "," + to_string(ys) + "),";
        //     }
        //     zs = via2;
        //     if(zs & 1) {
        //         for(int y = min(y1, y2); y <= max(y1, y2); y++) {
        //             res += cost_cube[UTIL::get_enc(zs, xs, y)];
        //             ptres += "(" + to_string(zs) + "," + to_string(xs) + "," + to_string(y) + "),";
        //         }
        //         ys = y2;
        //     } else {
        //         for(int x = min(x1, x2); x <= max(x1, x2); x++) {
        //             res += cost_cube[UTIL::get_enc(zs, x, ys)];
        //             ptres += "(" + to_string(zs) + "," + to_string(x) + "," + to_string(ys) + "),";
        //         }
        //         xs = x2;
        //     }   
        // }
        // if(zs != z2) {
        //     for(int i = min(zs, z2) + 1; i < max(zs, z2); i++) {
        //         res += cost_cube[UTIL::get_enc(i, xs, ys)];
        //         ptres += "(" + to_string(i) + "," + to_string(xs) + "," + to_string(ys) + "),";
        //     }
        //     zs = z2;
        //     res += cost_cube[UTIL::get_enc(zs, xs, ys)];
        // }
        // SPDLOG_ASSERT(zs == z2 && xs == x2 && ys == y2, "");
        // spdlog::info("lpat path is {}", ptres);
        // spdlog::info("res of lpat dist {} via1 {} via2 {}", res, via1, via2);
    }
    __device__ gcb_dt tmin(gcb_dt &val, unsigned long long int &val1) {
        return min(val, __longlong_as_double(val1));
    }
    // __device__ void atomicMin_float(gcb_dt* address, gcb_dt &val) {
    //     unsigned  int *address_as_ull = (unsigned int *)address;
    //     unsigned int old = *address_as_ull, assumed;
    //     do {
    //         assumed = old;
    //         old = atomicCAS(address_as_ull, assumed, __float_as_int(fmin(val, __int_as_float(assumed))));
    //     } while (assumed != old);
    // }
    __device__ void atomicMin_double(gcb_dt* address, gcb_dt &val) {
        #if WT_FLOAT == 0
        unsigned long long int* address_as_ull = (unsigned long long int*)address;
        unsigned long long int old = *address_as_ull, assumed;
        do {
            assumed = old;
            // double assumed_val = __longlong_as_double(assumed);
            // double min_val = fmin(val, assumed_val);
            old = atomicCAS(address_as_ull, assumed, __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
            // old = atomicCAS(address_as_ull, assumed, __double_as_longlong(tmin(val, (assumed))));
        } while (assumed != old);
        #else
        unsigned  int *address_as_ull = (unsigned int *)address;
        unsigned int old = *address_as_ull, assumed;
        do {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed, __float_as_int(fmin(val, __int_as_float(assumed))));
        } while (assumed != old);
        #endif
    }
    void get_lpat_res1(in_t &enc1, in_t enc2, lpat_dt &pat, gcb_dt &res) {
        gvia_dt via1, via2, z1, z2;
        gcd_dt x1, x2, y1, y2;
        gvia_dt val = pat;
        res = 0;
        get_zdec(val, via1, via2);
        UTIL::get_dec_all<CUBE_ENC>(enc1, z1, x1, y1);
        UTIL::get_dec_all<CUBE_ENC>(enc2, z2, x2, y2);
        gcd_dt xs=x1, ys=y1;
        gvia_dt zs=z1;
        SPDLOG_ASSERT(via1 > 0 && via1 < Z1 && via2 >= 0 && via2 < Z2, "");
        map<in_t, bool> mp;
        mp.clear();
        auto add_test_vias = [&mp](gvia_dt &zs, gvia_dt &z2, gcd_dt &xs, gcd_dt &ys) {
            for(gvia_dt i = min(zs, z2); i <= max(zs, z2); i++) {
                mp[UTIL::get_enc<VIA_ENC>(i, xs, ys)] = 1;
            }
            zs = z2;
        };
        auto add_test_wires = [&mp](gvia_dt &zs, gcd_dt &xs, gcd_dt &ys, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2) {
            SPDLOG_ASSERT(zs == z2, "zs {} z {}", zs, z2);
            genc_dt enca = UTIL::get_enc<WIRE_ENC>(zs, xs, ys), encb;
            if(zs & 1) {
                encb = UTIL::get_enc<WIRE_ENC>(zs, xs, y2); ys = y2;
            } else {
                encb = UTIL::get_enc<WIRE_ENC>(zs, x2, ys); xs = x2;
            }
            for(genc_dt i = min(enca, encb); i <= max(enca, encb); i += Z) {
                mp[i] = 1;
            }
        };
        add_test_vias(zs, via1, xs, ys);
        add_test_wires(zs, xs, ys, zs, x2, y2);
        via2 = (zs & 1) ? 2*via2: via2*2+1;
        SPDLOG_ASSERT(via2 < Z, "via2 {} Z {}", via2, Z);
        if(zs != z2 || xs != x2 || ys != y2) {
            if(via2 > 0) {
                add_test_vias(zs, via2, xs, ys);
                add_test_wires(zs, xs, ys, zs, x2, y2);
            }   
            
        }
        if(zs != z2) {
            add_test_vias(zs, z2, xs, ys);
        }
        res = 0;
        SPDLOG_ASSERT(z2 == zs && ys == y2 && xs == x2, "[zs,xs,ys] [{},{},{}] [z2,x2,y2] [{},{},{}]", zs, xs, ys, z2, x2, y2);
        for(auto mp_nd : mp) {
            res += wcost_cube[mp_nd.first];
        }
        mp.clear();
    }
    // void get_lpat_res2(in_t &enc1, in_t enc2, lpat_dt &pat, gcb_dt &res) {
    //     gcd_dt via1, via2, x1, x2, y1, y2, z1, z2;
    //     gvia_dt val = pat;
    //     res = 0;
    //     fr_t wcost = 0, vcost = 0;
    //     get_zdec(val, via1, via2);
    //     UTIL::get_dec_all(enc1, z1, x1, y1);
    //     UTIL::get_dec_all(enc2, z2, x2, y2);
    //     in_t xs=x1, ys=y1, zs=z1;
    //     SPDLOG_ASSERT(via1 > 0 && via1 < Z1 && via2 >= 0 && via2 < Z2, "");
    //     if(via1 != z1) {
    //         wcost += cost_cube[UTIL::get_enc(zs, xs, ys)];
    //     }
    //     add_vias(zs, via1, xs, ys, vcost);
    //     add_wires(wcost, zs, xs, ys, zs, x2, y2);
    //     // if(x2 != xs || ys != y2) LUL::add_wires(wcost, zs, xs, ys, zs, x2, y2);
    //     via2 = (zs & 1) ? 2*via2: via2*2+1;
    //     SPDLOG_ASSERT(via2 < Z, "via2 {} Z {}", via2, Z);
    //     if(zs != z2 || xs != x2 || ys != y2) {
    //         if(via2 > 0) {
    //             add_vias(zs, via2, xs, ys, vcost);
    //             add_wires(wcost, zs, xs, ys, zs, x2, y2);
    //         }            
    //         // if(x2 != xs || ys != y2) LUL::add_wires(wcost, zs, xs, ys, zs, x2, y2);
    //         if(zs != z2) {
    //             add_vias(zs, z2, xs, ys, vcost);                
    //             wcost += cost_cube[UTIL::get_enc(z2, x2, y2)];
    //         }
    //     }
    //     SPDLOG_ASSERT((via1 & 1) != (via2 & 1), "");
    //     SPDLOG_ASSERT(z2 == zs && ys == y2 && xs == x2, "[zs,xs,ys] [{},{},{}] [z2,x2,y2] [{},{},{}]", zs, xs, ys, z2, x2, y2);
    //     res = vcost + wcost;
    //  }
     void get_edge_ids(int bt_id, vector<ledid_dt> &edge_ids) {
        edge_ids.clear();
        #if MST_TREE == 0
        int blen = npart_len[bt_id], bofst = npart_ofst[bt_id];
        long pofst = 0;
        long num_prs;
        UTIL::get_num_pts_batch(num_prs, bt_id);
        // for(int i = 0; i < num_prs; i++) {
        //     edge_ids.push_back(i);
        // }
        for(int i = bofst; i < bofst + blen; i++) {
            int cnet_id = npart_ids[i];
            net_t cnet = nets[cnet_id];
            int num_pins = cnet.pins.size();
            for(long j = pofst; j < pofst + num_pins - 1; j++) {
                edge_ids.push_back(j);
            }
            pofst += ((num_pins * (num_pins - 1)) >> 1);
        }
        #else
        edge_ids.resize(MSC::mst_edges.size());
        for(int i = 0; i < edge_ids.size(); i++) {
            edge_ids[i] = (MSC::mst_edges[i] >> MST_NETID_BITS);
        }
        #endif
    }
    __device__ void gpu_add_wire(pth_dt *pth_mat, gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2) {
        if((z1 & 1)) {
            if(y1 < y2) {
                // pth_mat[UTIL::get_gpu_enc(z1, x1, y1)] += 1;
                atomicAdd(pth_mat + UTIL::get_gpu_enc<WIRE_ENC>(z1, x1, y1), 1);
            } else {
                // pth_mat[UTIL::get_gpu_enc(z1, x1, y2)] += 1;
                atomicAdd(pth_mat + UTIL::get_gpu_enc<WIRE_ENC>(z1, x1, y2), 1);
            }
            y1 = y2;
        } else {
            if(x1 < x2) {
                // pth_mat[UTIL::get_gpu_enc(z1, x1, y1)] += 1;
                atomicAdd(pth_mat + UTIL::get_gpu_enc<WIRE_ENC>(z1, x1, y1), 1);
            } else {
                // pth_mat[UTIL::get_gpu_enc(z1, x2, y1)] += 1;
                atomicAdd(pth_mat + UTIL::get_gpu_enc<WIRE_ENC>(z1, x2, y1), 1);
            }
            x1 = x2;
        }
    }
    void add_wire(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2) {
        if((z1 & 1)) {
            if(y1 < y2) {
                LPC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x1, y1)] += 1;
            } else {
                LPC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x1, y2)] += 1;
            }
            y1 = y2;
        } else {
            if(x1 < x2) {
                LPC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x1, y1)] += 1;
            } else {
                LPC::pth_mat[UTIL::get_enc<WIRE_ENC>(z1, x2, y1)] += 1;
            }
            x1 = x2;
        }
    }
    void add_via(gvia_dt &z1, gvia_dt &z2, gcd_dt &x1, gcd_dt &y1) {
        for(gvia_dt i = min(z2, z1)+1; i < max(z2, z1); i++) {
            genc_dt cenc = UTIL::get_enc<VIA_ENC>(i, x1, y1);
            if(LPC::pth_mat[cenc] == 0) {
                LPC::pth_mat[cenc] = -1;
            }
        }
    }
    __device__ void gpu_add_via(pth_dt *pth_mat, gvia_dt &z1, gvia_dt &z2, gcd_dt &x1, gcd_dt &y1) {
        for(gvia_dt i = min(z2, z1)+1; i < max(z2, z1); i++) {
            genc_dt cenc = UTIL::get_gpu_enc<VIA_ENC>(i, x1, y1);
            if(pth_mat[cenc] == 0) {
                pth_mat[cenc] = -1;
            }
        }
    }
    void add_lwire(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, gvia_dt &via1) {
        z1 = via1;
        add_wire(z1, x1, y1, z2, x2, y2);
    }
    __device__ void gpu_add_lwire(pth_dt *pth_mat, gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, gvia_dt &via1) {
        z1 = via1;
        gpu_add_wire(pth_mat, z1, x1, y1, z2, x2, y2);
    }
    void add_lvias(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, gvia_dt &via1) {
        add_via(z1, via1, x1, y1);
        z1 = via1;
        if((z1 & 1)) {
            y1 = y2;
        } else {
            x1 = x2;
        }
    }
    __device__ void gpu_add_lvias(pth_dt *pth_mat, gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, gvia_dt &via1) {
        gpu_add_via(pth_mat, z1, via1, x1, y1);
        z1 = via1;
        if((z1 & 1)) {
            y1 = y2;
        } else {
            x1 = x2;
        }
    }
    void add_wires(vector<long> &edge_ids) {
        for(long i = 0; i < edge_ids.size(); i++) {
            long ed_id = edge_ids[i];
            genc_dt pt_enc1 = LPC::hd_ptprs[(ed_id << 1)], pt_enc2 = LPC::hd_ptprs[(ed_id << 1) + 1];
            gvia_dt lpt = LPC::hd_lpat_pat[ed_id], via1, via2;
            get_zdec(lpt, via1, via2);
            gcd_dt x1, y1, x2, y2, xs, ys;
            gvia_dt z1, z2, zs;
            UTIL::get_dec_all<CUBE_ENC>(pt_enc1, z1, x1, y1); UTIL::get_dec_all<CUBE_ENC>(pt_enc2, z2, x2, y2);
            zs = z1; xs = x1; ys = y1;
            zs = z1; xs = x1; ys = y1;
            add_lwire(zs, xs, ys, z2, x2, y2, via1);
            #if LGPU_ASSERT == 1
            assert(zs == via1 && (xs == x2 || ys == y2));
            #endif
            via2 = (zs & 1) ? 2*via2: via2*2+1;
            if(zs != z2 || xs != x2 || ys != y2) {
                if(via2 > 0) {
                    add_lwire(zs, xs, ys, z2, x2, y2, via2);
                }
                zs = z2;
            }
            #if LGPU_ASSERT == 1
            assert(xs == x2 && zs == z2 && ys == y2);
            #endif
        }
    }
    __global__ void gpu_add_wires(pth_dt *pth_mat, ledid_dt *ed_ids, lpat_dt *lpat, genc_dt *edges_arr, long num_edges) {
        long tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= num_edges) {
            return;
        }
        long ed_id = tid / LWIRE_FACTOR;
        ed_id = ed_ids[ed_id];
        genc_dt pt_enc1 = edges_arr[(ed_id << 1)], pt_enc2 = edges_arr[(ed_id << 1) + 1];
        gvia_dt lpt = lpat[ed_id], pat1, pat2;
        gvia_dt via1, via2;
        get_gpu_zdec(lpt, via1, via2);
        gcd_dt x1, y1, x2, y2, x3, y3, xs, ys;
        gvia_dt z1, z2, z3, zs;
        UTIL::get_gpu_dec_all<CUBE_ENC>(pt_enc1, z1, x1, y1); UTIL::get_gpu_dec_all<CUBE_ENC>(pt_enc2, z2, x2, y2);
        //add gpu wire
        zs = z1; xs = x1; ys = y1;
        gpu_add_lwire(pth_mat, zs, xs, ys, z2, x2, y2, via1);
        #if LGPU_ASSERT == 1
        assert(zs == via1 && (xs == x2 || ys == y2));
        #endif
        via2 = (zs & 1) ? 2*via2: via2*2+1;
        if(zs != z2 || xs != x2 || ys != y2) {
            if(via2 > 0) {
                gpu_add_lwire(pth_mat, zs, xs, ys, z2, x2, y2, via2);
            }
            zs = z2;
        }
        #if LGPU_ASSERT == 1
        assert(xs == x2 && zs == z2 && ys == y2);
        #endif
    }
    void add_vias(vector<long> &edge_ids) {
        for(long i = 0; i < edge_ids.size(); i++) {
            long ed_id = edge_ids[i];
            genc_dt pt_enc1 = LPC::hd_ptprs[(ed_id << 1)], pt_enc2 = LPC::hd_ptprs[(ed_id << 1) + 1];
            gvia_dt lpt = LPC::hd_lpat_pat[ed_id], via1, via2;
            get_zdec(lpt, via1, via2);
            gcd_dt x1, y1, x2, y2, xs, ys;
            gvia_dt z1, z2, zs;
            UTIL::get_dec_all<CUBE_ENC>(pt_enc1, z1, x1, y1); UTIL::get_dec_all<CUBE_ENC>(pt_enc2, z2, x2, y2);
            zs = z1; xs = x1; ys = y1;
            add_lvias(zs, xs, ys, z2, x2, y2, via1);
            #if LCPU_ASSERT == 1
            assert(zs == via1 && (xs == x2 || ys == y2));
            #endif
            via2 = (zs & 1) ? 2*via2: via2*2+1;
            if(zs != z2 || xs != x2 || ys != y2) {
                if(via2 > 0) {
                    add_lvias(zs, xs, ys, z2, x2, y2, via2);
                }
                #if LCPU_ASSERT == 1
                assert(xs == x2 && ys == y2);
                #endif
                if(zs != z2) {
                    add_via(zs, z2, xs, ys);
                    zs = z2;
                }
            }
            #if LCPU_ASSERT == 1
            assert(xs == x2 && zs == z2 && ys == y2);
            #endif
        }
    }
    __global__ void gpu_add_vias(pth_dt *pth_mat, ledid_dt *ed_ids, lpat_dt *lpat, genc_dt *edges_arr, long num_edges) {
        long tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= num_edges) {
            return;
        }
        long ed_id = tid / LWIRE_FACTOR;
        ed_id = ed_ids[ed_id];
        genc_dt pt_enc1 = edges_arr[(ed_id << 1)], pt_enc2 = edges_arr[(ed_id << 1) + 1];
        gvia_dt lpt = lpat[ed_id], pat1, pat2;
        gvia_dt via1, via2;
        get_gpu_zdec(lpt, via1, via2);
        gcd_dt x1, y1, x2, y2, x3, y3, xs, ys;
        gvia_dt z1, z2, z3, zs;
        UTIL::get_gpu_dec_all<CUBE_ENC>(pt_enc1, z1, x1, y1); UTIL::get_gpu_dec_all<CUBE_ENC>(pt_enc2, z2, x2, y2);
        //add gpu wire
        zs = z1; xs = x1; ys = y1;
        gpu_add_lvias(pth_mat, zs, xs, ys, z2, x2, y2, via1);
        #if LGPU_ASSERT == 1
        assert(zs == via1 && (xs == x2 || ys == y2));
        #endif
        via2 = (zs & 1) ? 2*via2: via2*2+1;
        if(zs != z2 || xs != x2 || ys != y2) {
            if(via2 > 0) {
                gpu_add_lvias(pth_mat, zs, xs, ys, z2, x2, y2, via2);
            }
            #if LGPU_ASSERT == 1
            assert(xs == x2 && ys == y2);
            #endif
            if(zs != z2) {
                gpu_add_via(pth_mat, zs, z2, xs, ys);
                zs = z2;
            }
        }
        #if LGPU_ASSERT == 1
        assert(xs == x2 && zs == z2 && ys == y2);
        #endif
    }
    void get_pat_pfx() {
        int enc_id = 0;
        for(gvia_dt z = 0; z < Z; z++) {
            if(z & 1) {
                for(gcd_dt x = 0; x < X; x++) {
                    in_t pfx = 0, prv = 0;
                    for(gcd_dt y = 0; y < Y; y++) {
                        enc_id = UTIL::get_enc<CUBE_ENC>(z, x, y);
                        prv=LPC::pth_mat[enc_id]; 
                        pfx += prv;
                        LPC::pth_mat[enc_id] = pfx; enc_id++;
                    }
                }
            } else {
                for(gcd_dt y = 0; y < Y; y++) {
                    in_t pfx = 0, prv = 0;
                    for(gcd_dt x = 0; x < X; x++) {
                        enc_id = UTIL::get_enc<CUBE_ENC>(z, x, y);
                        prv=LPC::pth_mat[enc_id];
                        pfx += prv;
                        LPC::pth_mat[enc_id] = pfx; enc_id++;                        
                    }
                }
            }
        }
    }
    __global__ void get_gpu_pat_epfx(pth_dt *pth_mat, int nth, int mem_size) {
        extern __shared__ pth_dt smarr[];
        __shared__ gpfx_dt clyr;
        __shared__ int ctrk, enc_ofst, rem_val, thload;
        __shared__ pth_dt pfx_ofst, tpfx_ofst;
        if(threadIdx.x >= nth) {
            return;
        }
        if(threadIdx.x == 0) {
            clyr = 0; ctrk = 0; pfx_ofst = 0; tpfx_ofst = 0;
            for(gpfx_dt i = 0; i < DZ && ctrk <= blockIdx.x; i++) {
                ctrk += (i & 1) ? DX : DY;
                clyr++;
            }
            clyr--;
            // assert(clyr >= 0 && clyr < DZ);
            ctrk -= (clyr&1) ? DX:DY; assert(ctrk <= blockIdx.x);
            ctrk = blockIdx.x - ctrk; enc_ofst = DX * DY * clyr + (clyr & 1) ? DX * ctrk : DY * ctrk; rem_val = (clyr & 1) ? DY : DX;
            enc_ofst = DX * DY * clyr + ((clyr & 1) ? DY * ctrk : DX * ctrk);
            thload = CEIL_DIV(mem_size, nth);
            enc_ofst = 0;
        }
        __syncthreads();
        while(rem_val > 0) {
            __syncthreads();
            //load values to smarr
            gvia_dt zcd = clyr;
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                smarr[i] = pth_mat[UTIL::get_gpu_enc<WIRE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //get the prefix sum for each partition parallely
            pth_dt pfx=smarr[max(0, min(threadIdx.x * thload - 1, mem_size-1))], prev;
            if(threadIdx.x == 0) pfx = 0;
            __syncthreads();
            for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload, mem_size); i++) {
                prev=smarr[i]; smarr[i]=pfx; pfx+=prev;
            }
            __syncthreads();
            //propagate the pfx end value of each partition to other partitions
            if(threadIdx.x == 0) {
                pfx = smarr[thload - 1];
                for(int i = thload + thload - 1; i < mem_size - 1; i += thload) {
                    pfx += smarr[i]; smarr[i] = pfx;
                }
                if(thload - 1 != mem_size - 1) {
                    pfx += smarr[mem_size - 1]; smarr[mem_size - 1] = pfx;
                }
                gcd_dt xcd = ((zcd&1)?ctrk:enc_ofst + min(rem_val, mem_size) - 1), ycd = ((zcd&1)?enc_ofst + min(rem_val, mem_size) - 1:ctrk);
                tpfx_ofst = pth_mat[UTIL::get_gpu_enc<WIRE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //spread the propagated value to each partition
            if(threadIdx.x > 0) {
                pfx = smarr[min(threadIdx.x * thload - 1, mem_size-1)];
                for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload-1, mem_size-1); i++) {
                    smarr[i] += pfx;
                }
            }
            __syncthreads();
            //write back values to dram
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                // printf("helolo gpu world [zxy] [%d,%d,%d]\n", zcd, xcd, ycd);
                pth_mat[UTIL::get_gpu_enc<WIRE_ENC>(zcd, xcd, ycd)] = smarr[i] + pfx_ofst;
            }
            __syncthreads();
            //update enc_ofst and rem_val
            if(threadIdx.x == 0) {
                pfx_ofst += pfx + tpfx_ofst;
                // pfx_ofst += smarr[min(mem_size, rem_val) - 1];
                enc_ofst += mem_size;
                rem_val -= mem_size;
            }
            __syncthreads();
        }
    }
    __global__ void get_gpu_pat_ipfx(pth_dt *pth_mat, int nth, int mem_size) {
        extern __shared__ pth_dt smarr[];
        __shared__ gpfx_dt clyr;
        __shared__ int ctrk, enc_ofst, rem_val, thload;
        __shared__ pth_dt pfx_ofst, tpfx_ofst;
        if(threadIdx.x >= nth) {
            return;
        }
        if(threadIdx.x == 0) {
            clyr = 0; ctrk = 0; pfx_ofst = 0; tpfx_ofst = 0;
            for(gpfx_dt i = 0; i < DZ && ctrk <= blockIdx.x; i++) {
                ctrk += (i & 1) ? DX : DY;
                clyr++;
            }
            clyr--;
            // assert(clyr >= 0 && clyr < DZ);
            ctrk -= (clyr&1) ? DX:DY; assert(ctrk <= blockIdx.x);
            ctrk = blockIdx.x - ctrk; enc_ofst = DX * DY * clyr + (clyr & 1) ? DX * ctrk : DY * ctrk; rem_val = (clyr & 1) ? DY : DX;
            enc_ofst = DX * DY * clyr + ((clyr & 1) ? DY * ctrk : DX * ctrk);
            thload = CEIL_DIV(mem_size, nth);
            enc_ofst = 0;
        }
        __syncthreads();
        while(rem_val > 0) {
            __syncthreads();
            //load values to smarr
            gvia_dt zcd = clyr;
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                smarr[i] = pth_mat[UTIL::get_gpu_enc<WIRE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //get the prefix sum for each partition parallely
            pth_dt pfx=smarr[max(0, min(threadIdx.x * thload - 1, mem_size-1))], prev;
            if(threadIdx.x == 0) pfx = 0;
            __syncthreads();
            for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload, mem_size); i++) {
                prev=smarr[i]; smarr[i]=pfx; pfx+=prev;
            }
            __syncthreads();
            //propagate the pfx end value of each partition to other partitions
            if(threadIdx.x == 0) {
                pfx = smarr[thload - 1];
                for(int i = thload + thload - 1; i < mem_size - 1; i += thload) {
                    pfx += smarr[i]; smarr[i] = pfx;
                }
                if(thload - 1 != mem_size - 1) {
                    pfx += smarr[mem_size - 1]; smarr[mem_size - 1] = pfx;
                }
                gcd_dt xcd = ((zcd&1)?ctrk:enc_ofst + min(rem_val, mem_size) - 1), ycd = ((zcd&1)?enc_ofst + min(rem_val, mem_size) - 1:ctrk);
                tpfx_ofst = pth_mat[UTIL::get_gpu_enc<WIRE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //spread the propagated value to each partition
            if(threadIdx.x > 0) {
                pfx = smarr[min(threadIdx.x * thload - 1, mem_size-1)];
                for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload-1, mem_size-1); i++) {
                    smarr[i] += pfx;
                }
            }
            __syncthreads();
            //write back values to dram
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                // printf("helolo gpu world [zxy] [%d,%d,%d]\n", zcd, xcd, ycd);
                pth_mat[UTIL::get_gpu_enc<WIRE_ENC>(zcd, xcd, ycd)] += smarr[i] + pfx_ofst;
            }
            __syncthreads();
            //update enc_ofst and rem_val
            if(threadIdx.x == 0) {
                pfx_ofst += pfx + tpfx_ofst;
                // pfx_ofst += smarr[min(mem_size, rem_val) - 1];
                enc_ofst += mem_size;
                rem_val -= mem_size;
            }
            __syncthreads();
        }
    }
    __global__ void update_gpu_cap_cube(pth_dt *pth_mat, gcb_dt *cost_cube, gcb_dt *cap_cube, int num_items) {
        int tid = threadIdx.x + blockDim.x * blockIdx.x;
        if(tid >= num_items) {
            return;
        }
        pth_dt pth_val = pth_mat[tid];
        if(pth_val == -1) {
            gcb_dt cap = cap_cube[tid];
            cap = max(0.0, cap-VIA_WT_PDT);
            cap_cube[tid] = cap;
            cost_cube[tid] = UTIL::get_gpu_of_cost(cap);
        } else if(pth_val > 0) {
            gcb_dt cap = cap_cube[tid];
            cap = max(0.0, cap-WIRE_WT_PDT);
            cap_cube[tid] = cap;
            cost_cube[tid] = UTIL::get_gpu_of_cost(cap);
        }
    }
    void add_edge_update(int bt_id) {
        vector<ledid_dt> edge_ids;
        get_edge_ids(bt_id, edge_ids);
        long num_edges = edge_ids.size();
        add_wires(edge_ids);
        get_pat_pfx();
        add_vias(edge_ids);
    }
    void gpu_add_edge_update(int bt_id, long num_prs) {
        // vector<ledid_dt> edge_ids;
        // get_edge_ids(bt_id, edge_ids);
        // long num_edges = edge_ids.size();
        // CUDA_MALLOC_1D(LPC::d_led_id, num_edges, ledid_dt);
        // CUDAMEMCPYHD_1D(LPC::d_led_id, edge_ids.data(), num_edges, ledid_dt);
        // num_edges *= LWIRE_FACTOR;
        // long num_blks = CEIL_DIV(num_edges, NUM_GPU_THREADS);
        // gpu_add_wires<<<num_blks, NUM_GPU_THREADS>>>(LPC::d_pth_mat, LPC::d_led_id, LPC::d_lpat_pat, LPC::d_ptprs, num_edges);
        // cudaDeviceSynchronize();
        // int num_pfx_blocks;
        // UTIL::get_gpu_pfx_sum_blks(num_pfx_blocks);
        // int mem_size = max(X, Y), nth = 128;
        // // cudaFuncSetAttribute(get_gpu_pat_pfx, cudaFuncAttributeMaxDynamicSharedMemorySize, mem_size * sizeof(pth_dt));
        // get_gpu_pat_ipfx<<<num_pfx_blocks, NUM_GPU_THREADS, mem_size * sizeof(pth_dt)>>>(LPC::d_pth_mat, nth, mem_size);
        // cudaDeviceSynchronize();
        // num_edges = edge_ids.size() * LVIA_FACTOR;
        // num_blks = CEIL_DIV(num_edges, NUM_GPU_THREADS);
        // gpu_add_vias<<<num_blks, NUM_GPU_THREADS>>>(LPC::d_pth_mat, LPC::d_led_id, LPC::d_lpat_pat, LPC::d_ptprs, num_edges);
        // cudaDeviceSynchronize();
        // num_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
        // update_gpu_cap_cube<<<num_blks, NUM_GPU_THREADS>>>(LPC::d_pth_mat, d_wcost_cube, d_cube_data, CSIZE);
        // cudaDeviceSynchronize();
        // UTIL::get_gpu_pfx_sum_blks(num_pfx_blocks);
        // UTIL::d_get_pfx_cost_cube<<<num_pfx_blocks, NUM_GPU_THREADS, sizeof(gcb_dt) * mem_size>>>(d_wcost_cube, d_pfx_cost_cube, nth, mem_size);
        // cudaDeviceSynchronize();
        // cudaFree(LPC::d_led_id);
        // edge_ids.clear();
    }
    void write_wire(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gcd_dt &x2, gcd_dt &y2) {
        if(z1 & 1) {
            UTIL::add_vertical_cpu(z1, UTIL::get_gran_cd(x1), UTIL::get_gran_cd(y1), UTIL::get_gran_cd(y2));
            y1 = y2;
        } else {
            UTIL::add_horizontal_cpu(z1, UTIL::get_gran_cd(y1), UTIL::get_gran_cd(x1), UTIL::get_gran_cd(x2));
            x1 = x2;
        }
    }
    void write_pattern(long &ed_id, int &nid) {
        pt_prs cptprs = rd_enc[ed_id];
        gvia_dt lpat = LPC::hd_lpat_pat[ed_id], via1, via2, enc1, enc2;
        gcd_dt x1, y1, x2, y2, xs, ys;
        gvia_dt z1, z2, zs;
        enc1 = cptprs.enc1; enc2 = cptprs.enc2;
        get_zdec(lpat, via1, via2);
        UTIL::get_dec_all<CUBE_ENC>(enc1, z1, x1, y1);
        UTIL::get_dec_all<CUBE_ENC>(enc2, z2, x2, y2);
        zs = z1; xs = x1; ys = y1;
        UTIL::add_vias_cpu(UTIL::get_gran_cd(xs), UTIL::get_gran_cd(ys), zs, via1);
        zs = via1;
        write_wire(zs, xs, ys, x2, y2);
        via2 = (zs & 1) ? 2*via2: via2*2+1;
        if(zs != z2 || xs != x2 || ys != y2) {
            if(via2 > 0) {
                UTIL::add_vias_cpu(UTIL::get_gran_cd(xs), UTIL::get_gran_cd(ys), zs, via2);
                zs = via2;
                write_wire(zs, xs, ys, x2, y2);
            }            
        }
        if(zs != z2) {
            UTIL::add_vias_cpu(UTIL::get_gran_cd(xs), UTIL::get_gran_cd(ys), zs, z2);
            zs = z2;
        }
        #if FWRITE_ASSERT == 1
        SPDLOG_ASSERT(zs == z2 && xs == x2 && ys == y2, "[zs,xs,ys] [{},{},{}] [z2,x2,y2] [{},{},{}]", zs, xs, ys, z2, x2, y2);
        #endif
    }
}
