#include "mst_common.h"
#include "gmst_util.hpp"
#include "../l_pattern_files/lpattern_common.h"
#include "../l_pattern_files/lpattern_util.hpp"
#include "../bit_util.hpp"

using namespace std;
using namespace COMMON_DATA;

namespace GMSU = GMST_UTIL;
namespace MSC = MST_COMMON;
namespace LPC = LPAT_COMMON;
namespace LPU = LPAT_UTIL;
namespace BULS = BIT_UTIL_SPACE;

namespace MST_BIT_COMMIT {
    __inline__ __device__ void gpu_add_edg_cnt(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, int &cnid) {
        if(z1 != z2) {
            int zmin = min(z1, z2), zmax = max(z1, z2);
            for(int i = zmin + 1; i < zmax; i++) {
                MSC::d_net_edges_gcell_cnt[cnid]++;
            }
            z1 = z2;
            return;
        }
        if((z1 & 1)) {
            gcd_dt ymin = min(y1, y2), ymax = max(y1, y2);
            MSC::d_net_edges_gcell_cnt[cnid] += ymax - ymin + 1;
            y1 = y2;
        } else {
            gcd_dt xmin = min(x1, x2), xmax = max(x1, x2);
            MSC::d_net_edges_gcell_cnt[cnid] += xmax - xmin + 1;
            x1 = x2;
        }
    }
    __inline__ __device__ void gpu_hash_gcell(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, int &cnid) {
        if(z1 != z2) {
            int zmin = min(z1, z2) + 1, zmax = max(z1, z2) - 1;
            gcd_dt dist = zmax - zmin + 1;
            if(cnid < dist) {
                int cenc = UTIL::get_gpu_enc<WIRE_ENC>(zmin + cnid, x1, y1);
                d_bit_vis_data[cenc] = true;
            }
            cnid -= dist;
            z1 = z2;
            return;
        }
        if((z1 & 1)) {
            gcd_dt ymin = min(y1, y2), ymax = max(y1, y2);
            gcd_dt dist = ymax - ymin + 1;
            if(cnid < dist) {
                int cenc = UTIL::get_gpu_enc<WIRE_ENC>(z1, x1, ymin + cnid);
                d_bit_vis_data[cenc] = true;
            }
            y1 = y2;
            cnid -= dist;
        } else {
            gcd_dt xmin = min(x1, x2), xmax = max(x1, x2);
            gcd_dt dist = xmax - xmin + 1;
            if(cnid < dist) {
                int cenc = UTIL::get_gpu_enc<WIRE_ENC>(z1, xmin + cnid, y1);
                d_bit_vis_data[cenc] = true;
            }
            x1 = x2;
            cnid -= dist;
        }
    }
    __inline__ __device__ void gpu_bit_commit_wire(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, int &cnid) {
        int enc1, enc2;
        if((z1 & 1)) {
            gcd_dt ymin = min(y1, y2), ymax = max(y1, y2);
            enc1 = UTIL::get_gpu_enc<WIRE_ENC>(z1, x1, ymin), enc2 = UTIL::get_gpu_enc<WIRE_ENC>(z2, x1, ymax);
        } else {
            gcd_dt xmin = min(x1, x2), xmax = max(x1, x2);
            int enc1 = UTIL::get_gpu_enc<WIRE_ENC>(z1, xmin, y1), enc2 = UTIL::get_gpu_enc<WIRE_ENC>(z2, xmax, y1);
        }
        for(int i = enc1; i <= enc2; i++) {
            d_bit_vis_data[i] = true;
        }
    }
    __global__ void gpu_get_edg_gcell_cnt(int max_edges, int bofst, int blen) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_edges) {
            return;
        }
        nenc_dt cur_pat = MSC::d_net_mst_edges[tid];
        MSC::d_net_edges_gcell_cnt[tid] = 0;
        if(cur_pat == -1) {
            return;
        }
        int n1 = (cur_pat & NODE_MASK), n2 = int((cur_pat >> NODE_BITS) & NODE_MASK), cnet_id = int(cur_pat >> NNODE_BITS);
        long ed_id;
        if(n1 < n2) {
            int tmp_nd = n1;
            n1 = n2; n2 = tmp_nd;
        }
        ed_id = MSC::d_net_block_pfx1[cnet_id] + ((n1 * (n1 - 1)) >> 1) + n2;
        gvia_dt clpat = LPC::d_lpat_pat[ed_id], via1, via2;
        MSC::d_mst_lpat_pat[tid] = clpat;
        #if MST_MIGRATE_LPAT == 1
        MSC::d_mst_lpat_pat[tid] = clpat;
        #endif
        LPU::get_gpu_zdec(clpat, via1, via2);
        #if MST_GPU_ASSERT == 1
        assert(clpat >= 0 && via1 < DZ && via2 < DZ);
        #endif
        int ocnet_id = d_npart_ids[bofst + cnet_id];
        long onpfx = d_net_pins_pfx[ocnet_id];
        #if MST_GPU_ASSERT == 1
        int onlen = d_net_pins_len[ocnet_id];
        assert(n1 < onlen && n2 < onlen);
        #endif
        genc_dt enc_src = d_net_pins[onpfx + n1], enc_dst = d_net_pins[onpfx + n2];
        gvia_dt z1, z2, zs;
        gcd_dt x1, x2, y1, y2, xs, ys;
        UTIL::get_gpu_dec_all<CUBE_ENC>(enc_src, z1, x1, y1);
        UTIL::get_gpu_dec_all<CUBE_ENC>(enc_dst, z2, x2, y2);
        #if MST_GPU_ASSERT == 1
        assert(z1 < DZ && x1 < DX && y1 < DY);
        assert(z2 < DZ && x2 < DX && y2 < DY);
        #endif
        bool bth_diff = false;
        if(x1 != x2 && y1 != y2) {
            bth_diff = true;
        }
        if(z1 != via1) {
            gpu_add_edg_cnt(z1, x1, y1, via1, x1, y1, tid);
        }
        if(z1) gpu_add_edg_cnt(z1, x1, y1, z1, x2, y2, tid);
        via2 = ((z1 & 1) ? 2 * via2 : via2 * 2 + 1);
        if(bth_diff && via2 != 0) {
            #if MST_GPU_ASSERT == 1
            assert(via2 < DZ && via2 != 0);
            #endif
            if(z1 != via2) {
                gpu_add_edg_cnt(z1, x1, y1, via2, x1, y1, tid);
                gpu_add_edg_cnt(z1, x1, y1, z1, x2, y2, tid);
            }
        }
        if(z1 != z2) {
            gpu_add_edg_cnt(z1, x1, y1, z2, x1, y1, tid);
        }
        #if MST_GPU_ASSERT == 1
        assert(z1 == z2 && x1 == x2 && y1 == y2);
        #endif
    }
    __global__ void gpu_commit_edge(gccnt_dt max_edges_gcell_cnt, int max_edges, int bofst, int blen) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_edges_gcell_cnt) {
            return;
        }
        int st_id = 0, en_id = max_edges, md_id;
        while(st_id < en_id) {
            md_id = ((st_id + en_id) >> 1);
            if(MSC::d_net_edges_gcell_cnt_pfx[md_id] <= tid) {
                st_id = md_id + 1;
            } else {
                en_id = md_id;
            }
        }
        int edg_id = st_id - 1;
        #if MST_GPU_ASSERT == 1
            assert(edg_id < max_edges);
            assert(edg_id >= 0);
        #endif
        nenc_dt cur_pat = MSC::d_net_mst_edges[edg_id];
        if(cur_pat == -1) {
            edg_id -= 1;
        }
        cur_pat = MSC::d_net_mst_edges[edg_id];
        #if MST_GPU_ASSERT == 1
            assert(edg_id >= 0);
            assert(tid < MSC::d_net_edges_gcell_cnt_pfx[edg_id] + MSC::d_net_edges_gcell_cnt[edg_id]);
            assert(cur_pat != -1);
        #endif        
        int n1 = (cur_pat & NODE_MASK), n2 = int((cur_pat >> NODE_BITS) & NODE_MASK), cnet_id = int(cur_pat >> NNODE_BITS);
        long ed_id;
        if(n1 < n2) {
            int tmp_nd = n1;
            n1 = n2; n2 = tmp_nd;
        }
        ed_id = MSC::d_net_block_pfx1[cnet_id] + ((n1 * (n1 - 1)) >> 1) + n2;
        gvia_dt clpat = LPC::d_lpat_pat[ed_id], via1, via2;
        #if MST_MIGRATE_LPAT == 1
        MSC::d_mst_lpat_pat[tid] = clpat;
        #endif
        LPU::get_gpu_zdec(clpat, via1, via2);
        #if MST_GPU_ASSERT == 1
        assert(clpat >= 0 && via1 < DZ && via2 < DZ);
        #endif
        int ocnet_id = d_npart_ids[bofst + cnet_id];
        long onpfx = d_net_pins_pfx[ocnet_id];
        #if MST_GPU_ASSERT == 1
        int onlen = d_net_pins_len[ocnet_id];
        assert(n1 < onlen && n2 < onlen);
        #endif
        genc_dt enc_src = d_net_pins[onpfx + n1], enc_dst = d_net_pins[onpfx + n2];
        gvia_dt z1, z2, zs;
        gcd_dt x1, x2, y1, y2, xs, ys;
        UTIL::get_gpu_dec_all<CUBE_ENC>(enc_src, z1, x1, y1);
        UTIL::get_gpu_dec_all<CUBE_ENC>(enc_dst, z2, x2, y2);
        tid -= MSC::d_net_edges_gcell_cnt_pfx[edg_id];
        #if MST_GPU_ASSERT == 1
        assert(z1 < DZ && x1 < DX && y1 < DY);
        assert(z2 < DZ && x2 < DX && y2 < DY);
        assert(tid < MSC::d_net_edges_gcell_cnt[edg_id]);
        assert(tid >= 0);
        #endif
        bool bth_diff = false;
        if(x1 != x2 && y1 != y2) {
            bth_diff = true;
        }
        if(z1 != via1 && tid >= 0) {
            gpu_hash_gcell(z1, x1, y1, via1, x1, y1, tid);
        }
        if(z1) gpu_hash_gcell(z1, x1, y1, z1, x2, y2, tid);
        via2 = ((z1 & 1) ? 2 * via2 : via2 * 2 + 1);
        if(bth_diff && via2 != 0) {
            #if MST_GPU_ASSERT == 1
            assert(via2 < DZ && via2 != 0);
            #endif
            if(z1 != via2 && tid >= 0) {
                gpu_hash_gcell(z1, x1, y1, via2, x1, y1, tid);
                if(tid >= 0)
                    gpu_hash_gcell(z1, x1, y1, z1, x2, y2, tid);
            }
        }
        if(z1 != z2 && tid >= 0) {
            gpu_hash_gcell(z1, x1, y1, z2, x1, y1, tid);
        }
        z1 = z2;
        #if MST_GPU_ASSERT == 1
        assert(tid < 0);
        #endif
    }
    void bit_commit_edges(int &max_nodes_len) {
        CUDA_MEMSET(d_bit_vis_data, 0, CSIZE, bool);
        int num_blks = CEIL_DIV(max_nodes_len, NUM_GPU_THREADS);
        gpu_get_edg_gcell_cnt<<<num_blks, NUM_GPU_THREADS>>>(max_nodes_len, MSC::BOFST, MSC::BLEN);
        cudaDeviceSynchronize();
        thrust::device_ptr<gccnt_dt> d_edge_gcnt_ptr = thrust::device_pointer_cast(MSC::d_net_edges_gcell_cnt);
        thrust::device_ptr<gccnt_dt> d_edge_gcnt_optr = thrust::device_pointer_cast(MSC::d_net_edges_gcell_cnt_pfx);
        thrust::exclusive_scan(d_edge_gcnt_ptr, d_edge_gcnt_ptr + max_nodes_len, d_edge_gcnt_optr);
        cudaDeviceSynchronize();
        gccnt_dt gcnt_pfx, gcnt;
        CUDAMEMCPYDH_1D(&gcnt, MSC::d_net_edges_gcell_cnt + max_nodes_len - 2, 1, gccnt_dt);
        CUDAMEMCPYDH_1D(&gcnt_pfx, MSC::d_net_edges_gcell_cnt_pfx + max_nodes_len - 2, 1, gccnt_dt);
        gccnt_dt max_edg_cell_cnt = gcnt_pfx + gcnt;
        num_blks = CEIL_DIV(max_edg_cell_cnt, NUM_GPU_THREADS);
        gpu_commit_edge<<<num_blks, NUM_GPU_THREADS>>>(max_edg_cell_cnt, max_nodes_len, MSC::BOFST, MSC::BLEN);
        cudaDeviceSynchronize();
        num_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
        BULS::gpu_upd_bit_cost_cube<<<num_blks, NUM_GPU_THREADS>>>(CSIZE);
        cudaDeviceSynchronize();
    }
}
