#pragma once
#include "./common.h"
#include "./util.hpp"
#include "./l_pattern_files/lpattern_util.hpp"
#include "./l_pattern_files/lpattern_common.h"

using namespace std;
using namespace COMMON_DATA;
namespace LPC = LPAT_COMMON;
namespace LUL = LPAT_UTIL;

namespace BIT_UTIL_SPACE {
    __inline__ __device__ gcb_dt get_gpu_of_cost(gcb_dt cap) {
        return (expf(-cap * 0.5) * OFCOST_FACTOR);
    }
    __inline__ __device__ gcb_dt get_gpu_of_cost_inc(gcb_dt &cap, float inc_val) {
        return get_gpu_of_cost(cap - inc_val) - get_gpu_of_cost(cap);
    }
    __inline__ __device__ gcb_dt gpu_query_bit_cost(gvia_dt z, gcd_dt x, gcd_dt y) {
        gcb_dt res = 0;
        if(z & 1) {
            y++;
            while(y > 0) {
                res += d_cost_data_bit[UTIL::get_gpu_enc<WIRE_ENC>(z, x, y - 1)];
                y = y - (y & (-y));
            }
        } else {
            x++;
            while(x > 0) {
                res += d_cost_data_bit[UTIL::get_gpu_enc<WIRE_ENC>(z, x - 1, y)];
                x = x - (x & (-x));
            }
        }
        return res;
    }
    __inline__ __device__ gcb_dt gpu_bit_range_query(gvia_dt z1, gcd_dt x1, gcd_dt y1, gcd_dt x2, gcd_dt y2) {
        return gpu_query_bit_cost(z1, x2, y2) - gpu_query_bit_cost(z1, x1, y1);
    }
    __inline__ __device__ void gpu_update_bit_cost(gvia_dt z, gcd_dt x, gcd_dt y, gcb_dt &val) {
        if(z & 1) {
            y++;
            while(y <= DY) {
                genc_dt cenc = UTIL::get_gpu_enc<WIRE_ENC>(z, x, y - 1);
                atomicAdd(d_cost_data_bit + cenc, val);
                y = y + (y & (-y));
            }
        } else {
            x++;
            while(x <= DX) {
                genc_dt cenc = UTIL::get_gpu_enc<WIRE_ENC>(z, x - 1, y);
                atomicAdd(d_cost_data_bit + cenc, val);
                x = x + (x & (-x));
            }
        }
    }
    __global__ void gpu_init_bit_cost_cube(int max_gcells) {
        int tid = threadIdx.x + blockDim.x * blockIdx.x;
        if(tid >= max_gcells) {
            return;
        }
        d_cost_data_bit[tid] = 0;
    }
    __global__ void gpu_build_init_bit_cost_cube(int max_gcells) {
        int tid = threadIdx.x + blockDim.x * blockIdx.x;
        if(tid >= max_gcells) {
            return;
        }
        gvia_dt z1;
        gcd_dt x1, y1;
        genc_dt cenc = tid;
        UTIL::get_gpu_dec_all<WIRE_ENC>(cenc, z1, x1, y1);
        gcb_dt val = get_gpu_of_cost(d_cube_data_bit[tid]);
        gpu_update_bit_cost(z1, x1, y1, val);
    }
    __global__ void gpu_upd_bit_cost_cube(int max_gcells) {
        int tid = threadIdx.x + blockDim.x * blockIdx.x;
        if(tid >= max_gcells) {
            return;
        }
        if(d_bit_vis_data[tid] == false) {
            return;
        }
        gvia_dt z1;
        gcd_dt x1, y1;
        genc_dt cenc = tid;
        UTIL::get_gpu_dec_all<WIRE_ENC>(cenc, z1, x1, y1);
        gcb_dt val = get_gpu_of_cost_inc(d_cube_data_bit[tid], 1.0);
        d_cube_data_bit[tid]--;
        gpu_update_bit_cost(z1, x1, y1, val);
    }
    __inline__ __device__ gcb_dt gpu_get_bit_wire_cost(gcb_dt &wcost, gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gcd_dt &x2, gcd_dt &y2) {
        if(z1 & 1) {
            if(y1 < y2) {
                wcost += gpu_bit_range_query(z1, x1, y1 - 1, x1, y2);
            } else if(y1 > y2) {
                wcost += gpu_bit_range_query(z1, x1, y2 - 1, x1, y1);
            }
            y1 = y2;
        } else {
            if(x1 < x2) {
                wcost += gpu_bit_range_query(z1, x1 - 1, y1, x2, y1);
            } else if(x1 > x2) {
                wcost += gpu_bit_range_query(z1, x2 - 1, y1, x1, y1);
            }
            x1 = x2;
        }
    }
    __inline__ __device__ gcb_dt gpu_bit_get_via_cost(gcb_dt &wcost, gvia_dt &z1, gvia_dt &z2, gcd_dt &x1, gcd_dt &y1) {
        gvia_dt zmin = min(z1, z2), zmax = max(z1, z2);
        for(int i = zmin; i < zmax; i++) {
            if(i & 1) {
                wcost += gpu_bit_range_query(i, x1, y1 - 1, x1, y1);
            } else {
                wcost += gpu_bit_range_query(i, x1 - 1, y1, x1, y1);
            }
        }
    }
    template<gvia_dt LPATS_PERBLK>
    __global__ void get_gpu_lpat_dist(lpat_dt *lpat, genc_dt *pt_prs, gcb_dt *lpat_res, int MAX_PTS) {
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
            gcb_dt ccost = 0;
            if(via1 == 0) {
                break;
            }
            bool bth_diff = false;
            if(xs != x2 && ys != y2) {
                bth_diff = true;
            }
            if(zs != via1) {
                gpu_bit_get_via_cost(ccost, zs, via1, xs, ys);
                zs = via1;
            }
            if(zs) {
                gpu_get_bit_wire_cost(ccost, zs, xs, ys, x2, y2);
            }
            via2 = ((zs & 1) ? 2 * via2 : via2 * 2 + 1);
            if(bth_diff && via2 != 0) {
                if(zs != via2) {
                    gpu_bit_get_via_cost(ccost, zs, via2, xs, ys);
                    zs = via2;
                    gpu_get_bit_wire_cost(ccost, zs, xs, ys, x2, y2);
                }
            }
            if(zs != z2) {
                gpu_bit_get_via_cost(ccost, zs, z2, xs, ys);
                zs = z2;
            }
            // if(zs & 1) {
            //     ccost += gpu_bit_range_query(zs, x1, y1 - 1, x1, y1);
            // } else {
            //     ccost += gpu_bit_range_query(zs, x1 - 1, y1, x1, y1);
            // }
            if(xs != x2 || ys != y2 || zs != z2) {
                break;
            }
            c_cost = ccost;
            LUL::atomicMin_double(res+idx, c_cost);
            break;
        }
        __syncthreads();
        if(idx < LPATS_PERBLK && res[idx] == c_cost) {
            #if MST_GPU_ASSERT == 1
                assert(c_cost != MAX_OFGPU_COST);
            #endif
            idx = (blockIdx.x * LPATS_PERBLK) + idx;
            if(idx < MAX_PTS) {
                lpat_res[idx] = c_cost; lpat[idx] = cpat;
            }
        }
    }
    void init_bit_data() {
        TIME_();
        spdlog::info("initialising bit cube cost data");
        int num_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
        gpu_init_bit_cost_cube<<<num_blks, NUM_GPU_THREADS>>>(CSIZE);
        cudaDeviceSynchronize();
        gpu_build_init_bit_cost_cube<<<num_blks, NUM_GPU_THREADS>>>(CSIZE);
        cudaDeviceSynchronize();
    }
}
