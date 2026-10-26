#pragma once
#include "../common.h"
#include "../util.hpp"
#include "./lpattern_util.hpp"
#include "./lpattern_common.h"

using namespace std;
using namespace COMMON_DATA;
namespace LPC = LPAT_COMMON;
namespace LUL = LPAT_UTIL;

namespace LPATTERN {
    template<gvia_dt LPATS_PERBLK>
    __global__ void get_gpu_lpat_dist(lpat_dt *lpat, genc_dt *pt_prs, gcb_dt *lpat_res, gcb_dt *wcost_cube, gcb_dt *vcost_cube, int MAX_PTS) {
        __shared__ gcb_dt res[LPATS_PERBLK];
        __shared__ genc_dt encs[2 * LPATS_PERBLK];
        if((blockIdx.x * LPATS_PERBLK) + (threadIdx.x / LTHREAD_PER_EDGE) >= MAX_PTS) {
            return;
        }        
        if(threadIdx.x < 2 * LPATS_PERBLK) {
            int idx = (blockIdx.x * LPATS_PERBLK) * 2 + threadIdx.x;
            encs[threadIdx.x] = ((idx < MAX_PTS * 2)? pt_prs[idx]:-1);
            if(threadIdx.x < LPATS_PERBLK) res[threadIdx.x] = gcb_dt(MAX_OFGPU_COST);
        }
        __syncthreads();
        gvia_dt idx = threadIdx.x / DZ1Z2;
        gvia_dt cpat = threadIdx.x % DZ1Z2;
        gcb_dt c_cost = MAX_OFGPU_COST;
        while(idx < LPATS_PERBLK) {
            gvia_dt via1, via2, z1, z2, zs;
            gcd_dt x1, y1, x2, y2, xs, ys;
            UTIL::get_gpu_dec_all<CUBE_ENC>(encs[idx * 2], z1, x1, y1);
            UTIL::get_gpu_dec_all<CUBE_ENC>(encs[idx * 2 + 1], z2, x2, y2);
            LUL::get_gpu_zdec(cpat, via1, via2);
            // assert(via1 < 10 && via2 < 5);
            xs = x1; ys = y1; zs = z1;
            gcb_dt wcost = 0, vcost = 0;
            if(via1 == 0 || via1 >= DZ) {
                break;
            }
            if(zs != via1) {
                if(zs < via1) {
                    LUL::add_gpu_vias(zs, via1, xs, ys, vcost, vcost_cube);
                } else {
                    LUL::add_gpu_vias(via1, zs, xs, ys, vcost, vcost_cube);
                }                
                zs = via1;
            }
            LUL::add_gpu_wires(wcost, zs, xs, ys, zs, x2, y2, wcost_cube);
            via2 = ((zs & 1) ? 2 * via2 : via2 * 2 + 1);
            if(via2 >= DZ) {
                break;
            }
            if(zs != z2 || xs != x2 || ys != y2) {
                if(via2 > 0) {
                    if(zs != via2) {
                        if(zs < via2) {
                            LUL::add_gpu_vias(zs + 1, via2, xs, ys, vcost, vcost_cube);
                        } else {
                            LUL::add_gpu_vias(via2 + 1, zs, xs, ys, vcost, vcost_cube);
                        }
                        zs = via2;
                    }
                    LUL::add_gpu_wires(wcost, zs, xs, ys, zs, x2, y2, wcost_cube);
                }
                if(zs != z2) {
                    if(zs < z2) {
                        LUL::add_gpu_vias(zs + 1, z2, xs, ys, vcost, vcost_cube);
                    } else {
                        LUL::add_gpu_vias(z2 + 1, zs, xs, ys, vcost, vcost_cube);
                    }
                    zs = z2;
                }
                if(xs != x2 || ys != y2 || zs != z2) {
                    // #if LGPU_ASSERT == 1
                    // assert(xs == x2 && ys == y2 && zs == z2);
                    // #endif
                    break;
                }
            }
            #if LIVE_COST == 1
            c_cost = wcost + (vcost * VIA_WT_PDT);
            #else
            c_cost = abs(z1 - z2) + abs(x1 - x2) + abs(y1 - y2);
            #endif
            // c_cost = wcost;
            LUL::atomicMin_double(res+idx, c_cost);
            break;
        }
        __syncthreads();
        if(idx < LPATS_PERBLK && res[idx] == c_cost) {            
            idx = (blockIdx.x * LPATS_PERBLK) + idx;
            if(idx < MAX_PTS) {
                lpat_res[idx] = c_cost; lpat[idx] = cpat;
            }
        }
    }
    void get_lpat_dist(genc_dt &val1, genc_dt &val2, gcb_dt &res, lpat_dt &pat) {
        gcd_dt x1, y1, x2, y2;
        gvia_dt z1, z2;
        UTIL::get_dec_all<CUBE_ENC>(val1, z1, x1, y1); UTIL::get_dec_all<CUBE_ENC>(val2, z2, x2, y2);
        SPDLOG_ASSERT(z1<Z && z2 < Z && z1 >= 0 && z2 >= 0 && x1 < X && x1 >= 0 && y1 < Y && y1 >=0 && x2 < X && x2 >= 0 && y2 < Y && y2 >= 0, "[Z,X,Y]({},{},{}) ({},{},{}) ({},{},{})", Z, X, Y, z1, x1, y1, z2, x2, y2);
        pat = -1; res = MAX_OFCOST;
        for(gvia_dt i = 0; i < Z1*Z2; i++) {
            gvia_dt via1, via2;
            LUL::get_zdec(i, via1, via2);
            gcb_dt wcost = 0, vcost = 0;
            gcd_dt xs=x1, ys=y1;
            gvia_dt zs=z1;
            SPDLOG_ASSERT(via1 >= 0 && via1 < Z1 && via2 >= 0 && via2 < Z2, "");
            if(via1 == 0 || via1 >= Z) {
                continue;
            }
            if(zs != via1) {
                if(zs < via1) {
                    LUL::add_vias(zs, via1, xs, ys, vcost);
                } else {
                    LUL::add_vias(via1, zs, xs, ys, vcost);
                }                
                zs = via1;
            }
            LUL::add_wires(wcost, zs, xs, ys, zs, x2, y2);
            via2 = (zs & 1) ? 2*via2: via2*2+1;
            if(via2 >= Z) {
                continue;
            }
            if(zs != z2 || xs != x2 || ys != y2) {
                if(via2 > 0) {
                    if(via2 != zs) {
                        if(zs < via2) {
                            LUL::add_vias(zs + 1, via2, xs, ys, vcost);
                        } else {
                            LUL::add_vias(via2 + 1, zs, xs, ys, vcost);
                        }
                        zs = via2;
                    }
                    LUL::add_wires(wcost, zs, xs, ys, zs, x2, y2);
                }            
                // if(x2 != xs || ys != y2) LUL::add_wires(wcost, zs, xs, ys, zs, x2, y2);
                if(zs != z2) {
                    if(zs < z2) {
                        LUL::add_vias(zs + 1, z2, xs, ys, vcost);
                    } else {
                        LUL::add_vias(z2 + 1, zs, xs, ys, vcost);
                    }
                    zs = z2;
                }
            }
            if(xs != x2 || ys != y2 || zs != z2) {
                continue;
            }
            SPDLOG_ASSERT((via1 & 1) != (via2 & 1), "");
            SPDLOG_ASSERT(xs == x2 && ys == y2 && zs == z2, "[z,x,y] ({},{},{}) ({},{},{})", zs, xs, ys, z2, x2, y2);
            gcb_dt c_cost = wcost +  vcost;
            if(res > c_cost) {
                res = c_cost; pat = i;
            }
        }
        SPDLOG_ASSERT(pat != -1, "");
    }
    void get_lpat_djkst(in_t src, in_t dst, fr_t &res) {
        // memset(LPC::lpat_par, -1, sizeof(int) * CSIZE);
        // for(int i = 0; i < CSIZE; i++) {
        //     LPC::dcd_t cdt;
        //     cdt.dist = fr_t(MAX_OFCOST); cdt.cdir = 0; cdt.ndir_chng = 0;
        //     LPC::djkst_oval[i] = cdt;
        // }
        // in_t z1, x1, y1, z2, x2, y2;
        // UTIL::get_dec_all(src, z1, x1, y1); UTIL::get_dec_all(dst, z2, x2, y2);
        // auto cmp = [](const LPC::dcd_t &a, const LPC::dcd_t &b) {
        //     if(a.dist == b.dist) {
        //         if(a.z == b.z) {
        //             if(a.x == b.x) {
        //                 // if(a.y == b.y) {
        //                 //     if(a.ndir_chng == b.ndir_chng) {
        //                 //         return a.cdir < b.cdir;
        //                 //     }
        //                 //     return a.ndir_chng < b.ndir_chng;
        //                 // }
        //                 return a.y < b.y;
        //             }
        //             return a.x < b.x;
        //         }
        //         return a.z < b.z;
        //     }
        //     return a.dist < b.dist;
        // };
        // set<LPC::dcd_t, decltype(cmp)> pq(cmp);
        // pq.clear();
        // int num_pts = 0;
        // map<in_t, int> mp;
        // mp.clear();
        // auto add_xpts = [&pq](in_t xa, in_t xb, in_t ya) {            
        //     // spdlog::info("adding xpts xn {} xm {} ya {}", min(xa, xb), max(xa, xb), ya);
        //     for(in_t x = min(xa, xb); x <= max(xa, xb); x++) {
        //         for(in_t z = 0; z < Z; z++) {
        //             LPC::dcd_t dtt;
        //             dtt.x=x;dtt.y=ya;dtt.z=z;dtt.dist=fr_t(MAX_OFCOST);dtt.cdir=0;dtt.ndir_chng=0;
        //             pq.insert(dtt);
        //         }
        //     }            
        // };
        // auto add_ypts = [&pq](in_t ya, in_t yb, in_t xa) {            
        //     // spdlog::info("adding xpts ya {} yb {} xa {}", min(ya, yb), max(ya, yb), xa);
        //     for(in_t y = min(ya, yb); y <= max(ya, yb); y++) {
        //         for(in_t z = 0; z < Z; z++) {
        //             LPC::dcd_t dtt;
        //             dtt.x=xa;dtt.y=y;dtt.z=z;dtt.dist=fr_t(MAX_OFCOST);dtt.cdir=0;dtt.ndir_chng=0;
        //             pq.insert(dtt);
        //         }
        //     }            
        // };
        // auto is_safe = [&x1, &y1, &x2, &y2](in_t xa, in_t ya, in_t za) {
        //     return (xa == x1 || xa == x2 || ya == y1 || ya == y2) && za >= 0 && za < Z && xa >= min(x1, x2) && xa <= max(x1, x2) && ya >= min(y1, y2) && ya <= max(y1,y2);
        // };
        // auto add_pqpts = [&pq](LPC::dcd_t &pnd, in_t xa, in_t ya, in_t za) {
        //     in_t enc = UTIL::get_enc(za, xa, ya);
        //     // fr_t odist = LPC::lpat_dist[enc];
        //     LPC::dcd_t oval = LPC::djkst_oval[enc];
        //     fr_t cst = pnd.dist + cost_cube[enc];
        //     // spdlog::info("adding point ({},{},{})", xa, ya, za);
        //     LPC::dcd_t cdt = LPC::dcd_t{x:xa, y:ya, z:za, dist:oval.dist};
        //     int cdir = pnd.cdir;
        //     if(oval.dist > cst && pq.count(cdt)) {
        //         int ndir = LUL::gdir(pnd.z, pnd.x, pnd.y, za, xa, ya), ccng = pnd.ndir_chng;
        //         // LPC::dcd_t cdt = LPC::dcd_t{x:xa, y:ya, z:za, dist:odist};
        //         // SPDLOG_ASSERT(pq.count(cdt) == 1, "{} {}", odist, cst);
        //         if(ndir > 0 && cdir > 0) {
        //             if(ndir != cdir) {
        //                 ccng++;
        //             }
        //         } else {
        //             ndir = cdir + ndir;
        //         }
        //         if(ccng > 1) {
        //             return;
        //         }
        //         pq.erase(cdt); cdt.dist = cst; cdt.cdir = ndir; cdt.ndir_chng = ccng;
        //         pq.insert(cdt); LPC::djkst_oval[enc] = cdt;
        //         LPC::lpat_par[enc] = UTIL::get_enc(pnd.z, pnd.x, pnd.y);
        //     }
        // };
        // add_xpts(x1, x2, y1); add_xpts(x1, x2, y2); add_ypts(y1, y2, x1); add_ypts(y1, y2, x2);
        // LPC::dcd_t src_dt = {x:x1, y:y1, z:z1, dist:MAX_OFCOST, cdir:0, ndir_chng:0};
        // pq.erase(src_dt);
        // src_dt.dist = cost_cube[src];
        // // SPDLOG_ASSERT(pq.count(LPC::dcd_t{x:x2, y:y2, z:z2, dist:fr_t(MAX_OFCOST), cdir:0, ndir_chng:0}) == 1, "");
        // pq.insert(LPC::dcd_t{x:x1, y:y1, z:z1, dist:cost_cube[src], cdir:0, ndir_chng:0}); pq.insert(LPC::dcd_t{x:x2, y:y2, z:z2, dist:fr_t(MAX_OFCOST), cdir:0, ndir_chng:0});
        // LPC::djkst_oval[src].dist = cost_cube[src];
        // // spdlog::info("pq size {}", pq.size());
        // while(pq.size()) {
        //     LPC::dcd_t cur_elm = *(pq.begin());
        //     pq.erase(cur_elm);
        //     in_t cx = cur_elm.x, cy = cur_elm.y, cz = cur_elm.z;
        //     fr_t cdist = cur_elm.dist;
        //     in_t par = UTIL::get_enc(cz, cx, cy);
        //     // spdlog::info("parent point ({},{},{})", cx, cy, cz);
        //     if(cz > 0) {
        //         if(cz & 1) {
        //             if(is_safe(cx, cy+1, cz)) add_pqpts(cur_elm, cx, cy + 1, cz);
        //             if(is_safe(cx, cy-1, cz)) add_pqpts(cur_elm, cx, cy - 1, cz);
        //         } else {
        //             if(is_safe(cx + 1, cy, cz)) add_pqpts(cur_elm, cx + 1, cy, cz);
        //             if(is_safe(cx - 1, cy, cz)) add_pqpts(cur_elm, cx - 1, cy, cz);
        //         }
        //     }
        //     for(int i = 4; i < 6; i++) {
        //         if((cx == x1 && cy == y1) || (cx == x2 && cy == y2) || (cx == x1 && cy == y2) || (cx == x2 && cy == y1)) {
        //             if(is_safe(cx, cy, cz+1)) add_pqpts(cur_elm, cx, cy, cz + 1);
        //             if(is_safe(cx, cy, cz-1)) add_pqpts(cur_elm, cx, cy, cz - 1);
        //         }                
        //     }
        // }
        // res = LPC::djkst_oval[dst].dist;
        // // spdlog::info("dist value for dst is {}", res);
    }
}
