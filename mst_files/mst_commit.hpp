#include "mst_common.h"
#include "gmst_util.hpp"
#include "../l_pattern_files/lpattern_common.h"
#include "../l_pattern_files/lpattern_util.hpp"
#include "./mst_bit_commit.hpp"

using namespace std;
using namespace COMMON_DATA;

namespace GMSU = GMST_UTIL;
namespace MSC = MST_COMMON;
namespace LPC = LPAT_COMMON;
namespace LPU = LPAT_UTIL;
namespace MBC = MST_BIT_COMMIT;

namespace MST_COMMIT {
    __inline__ __device__ void gpu_commit_via(gvia_dt &z1, gvia_dt &z2, gcd_dt &x1, gcd_dt &y1, int &cnid) {
        for(int i = min(z1, z2); i < max(z1, z2); i++) {
            int enc = UTIL::get_gpu_enc<WIRE_ENC>(i, x1, y1);
            // if((d_pre_demand_data[enc]) > 0)
            //     continue;
            if((i & 1)) {
                if(y1 == 0) {
                    atomicAdd(d_cube_data + enc, -MST_VIA_FULL_WT);
                    #if EVAL_OFCOST == 1 && EVAL_VIA_OFCOST == 1
                    atomicAdd(d_eval_via_cnt + enc, 2);
                    #endif
                } else if(y1 == DY - 1) {
                    atomicAdd(d_cube_data + enc - 1, -MST_VIA_FULL_WT);
                    #if EVAL_OFCOST == 1 && EVAL_VIA_OFCOST == 1
                    atomicAdd(d_eval_via_cnt + enc - 1, 2);
                    #endif
                } else {
                    atomicAdd(d_cube_data + enc - 1, -MST_VIA_FRAC_WT);
                    atomicAdd(d_cube_data + enc, -MST_VIA_FRAC_WT);
                    #if EVAL_OFCOST == 1 && EVAL_VIA_OFCOST == 1
                    atomicAdd(d_eval_via_cnt + enc - 1, 1);
                    atomicAdd(d_eval_via_cnt + enc, 1);
                    #endif
                }
            } else {
                if(x1 == 0) {
                    atomicAdd(d_cube_data + enc, -MST_VIA_FULL_WT);
                    #if EVAL_OFCOST == 1 && EVAL_VIA_OFCOST == 1
                    atomicAdd(d_eval_via_cnt + enc, 2);
                    #endif
                } else if(x1 == DX - 1) {
                    atomicAdd(d_cube_data + enc - 1, -MST_VIA_FULL_WT);
                    #if EVAL_OFCOST == 1 && EVAL_VIA_OFCOST == 1
                    atomicAdd(d_eval_via_cnt + enc - 1, 2);
                    #endif
                } else {
                    atomicAdd(d_cube_data + enc - 1, -MST_VIA_FRAC_WT);
                    atomicAdd(d_cube_data + enc, -MST_VIA_FRAC_WT);
                    #if EVAL_OFCOST == 1 && EVAL_VIA_OFCOST == 1
                    atomicAdd(d_eval_via_cnt + enc - 1, 1);
                    atomicAdd(d_eval_via_cnt + enc, 1);
                    #endif
                }
            }
        }
        #if EVAL_OFCOST == 1 && EVAL_CONN_CHECK == 1
        for(gvia_dt i = min(z1, z2); i <= max(z1, z2); i++) {
            d_eval_vis_data[UTIL::get_gpu_enc<WIRE_ENC>(i, x1, y1)] = cnid+1;
        }
        #endif
    }
    __inline__ __device__ void gpu_commit_wire(gvia_dt &z1, gcd_dt &x1, gcd_dt &y1, gvia_dt &z2, gcd_dt &x2, gcd_dt &y2, int &cnid) {
        #if MST_GPU_ASSERT == 1
        assert(z1 == z2);
        #endif
        if((z1 & 1)) {
            gcd_dt ymin = min(y1, y2), ymax = max(y1, y2);
            int enc1 = UTIL::get_gpu_enc<WIRE_ENC>(z1, x1, ymin), enc2 = UTIL::get_gpu_enc<WIRE_ENC>(z2, x1, ymax);
            #if MST_GPU_ASSERT == 1
            assert(enc1 <= enc2);
            #endif
            atomicAdd(d_pre_demand_data + enc1, 1);
            // if(ymin + 1 < DY) {
            //     atomicAdd(d_pre_demand_data + enc1 + 1, 1);
            // }
            atomicAdd(d_pre_demand_data + enc2, -1);
            // if(enc2 + 1 < DY) {
            //     atomicAdd(d_pre_demand_data + enc2 + 1, -1);
            // }            
            #if EVAL_OFCOST == 1 && EVAL_CONN_CHECK == 1
            for(int i = enc1; i <= enc2; i++) {
                d_eval_vis_data[i] = cnid + 1;
            }
            #endif
            y1 = y2;
        } else {
            gcd_dt xmin = min(x1, x2), xmax = max(x1, x2);
            int enc1 = UTIL::get_gpu_enc<WIRE_ENC>(z1, xmin, y1), enc2 = UTIL::get_gpu_enc<WIRE_ENC>(z2, xmax, y1);
            #if MST_GPU_ASSERT == 1
            assert(enc1 <= enc2);
            #endif
            atomicAdd(d_pre_demand_data + enc1, 1);
            // if(xmin + 1 < DX) {
            //     atomicAdd(d_pre_demand_data + enc1 + 1, 1);
            // }
            atomicAdd(d_pre_demand_data + enc2, -1);
            // if(xmax + 1 < DX) {
            //     atomicAdd(d_pre_demand_data + enc2 + 1, -1);
            // }
            #if EVAL_OFCOST == 1 && EVAL_CONN_CHECK == 1
            for(int i = enc1; i <= enc2; i++) {
                d_eval_vis_data[i] = cnid + 1;
            }
            #endif
            x1 = x2;
        }
    }
    __global__ void gpu_add_vias(int max_edges, int bofst, int blen) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_edges) {
            return;
        }
        nenc_dt cur_pat = MSC::d_net_mst_edges[tid];
        if(cur_pat == -1) {
            return;
        }
        int n1 = (cur_pat & NODE_MASK), n2 = int((cur_pat >> NODE_BITS) & NODE_MASK), cnet_id = int(cur_pat >> NNODE_BITS);
        int npfx = MSC::d_net_nodes_cnt_pfx[cnet_id], nlen = MSC::d_net_nodes_cnt[cnet_id];
        assert(n1 < nlen && n2 < nlen);
        long ed_id;
        if(n1 < n2) {
            int tmp_nd = n1;
            n1 = n2; n2 = tmp_nd;
        }
        ed_id = MSC::d_net_block_pfx1[cnet_id] + ((n1 * (n1 - 1)) >> 1) + n2;
        #if MST_GPU_ASSERT == 1
            assert(ed_id < MSC::d_net_block_pfx1[cnet_id] + MSC::d_net_block_len1[cnet_id]);
        #endif
        #if MST_MIGRATE_LPAT == 0
        gvia_dt clpat = LPC::d_lpat_pat[ed_id], via1, via2;
        #else
        gvia_dt clpat = MSC::d_mst_lpat_pat[tid], via1, via2;
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
        gpu_commit_via(z1, via1, x1, y1, cnet_id);
        z1 = via1;
        if(z1 & 1) {
            y1 = y2;
        } else {
            x1 = x2;
        }
        via2 = ((z1 & 1) ? 2 * via2 : via2 * 2 + 1);
        #if MST_GPU_ASSERT == 1
        assert(via2 < DZ);
        #endif
        if(via2 > 0) {
            gpu_commit_via(z1, via2, x1, y1, cnet_id);
            z1 = via2;
            if(z1 & 1) {
                y1 = y2;
            } else {
                x1 = x2;
            }
        }
        if(z1 != z2) {
            gpu_commit_via(z1, z2, x1, y1, cnet_id);
            z1 = z2;
        }
        #if MST_GPU_ASSERT == 1
        assert(z1 == z2 && x1 == x2 && y1 == y2);
        #endif
    }

    __global__ void gpu_add_wires(int max_edges, int bofst, int blen) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_edges) {
            return;
        }
        nenc_dt cur_pat = MSC::d_net_mst_edges[tid];
        if(cur_pat == -1) {
            return;
        }
        int n1 = (cur_pat & NODE_MASK), n2 = int((cur_pat >> NODE_BITS) & NODE_MASK), cnet_id = int(cur_pat >> NNODE_BITS);
        #if MST_GPU_ASSERT == 1
        int npfx = MSC::d_net_nodes_cnt_pfx[cnet_id], nlen = MSC::d_net_nodes_cnt[cnet_id];
        assert(n1 < nlen && n2 < nlen);
        #endif
        long ed_id;
        if(n1 < n2) {
            int tmp_nd = n1;
            n1 = n2; n2 = tmp_nd;
        }
        ed_id = MSC::d_net_block_pfx1[cnet_id] + ((n1 * (n1 - 1)) >> 1) + n2;
        #if MST_GPU_ASSERT == 1
            // if(ed_id >= MSC::d_net_block_pfx1[cnet_id] + MSC::d_net_block_len1[cnet_id]) {
            //     printf("[nid,n1,n2,ed_id] [%d,%d,%d,%d]\n", cnet_id, n1, n2, ed_id);
            // }
            assert(ed_id < MSC::d_net_block_pfx1[cnet_id] + MSC::d_net_block_len1[cnet_id]);
        #endif
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
        #if MST_GPU_ASSERT == 1
        assert(z1 < DZ && x1 < DX && y1 < DY);
        assert(z2 < DZ && x2 < DX && y2 < DY);
        #endif
        z1 = via1;
        gpu_commit_wire(z1, x1, y1, z1, x2, y2, cnet_id);
        via2 = ((z1 & 1) ? 2 * via2 : via2 * 2 + 1);
        if(via2 > 0) {
            #if MST_GPU_ASSERT == 1
            assert(via2 < DZ && via2 != 0);
            #endif
            z1 = via2;
            gpu_commit_wire(z1, x1, y1, z1, x2, y2, cnet_id);
        }
        z1 = z2;
        #if MST_GPU_ASSERT == 1
        assert(z1 == z2 && x1 == x2 && y1 == y2);
        #endif
    }
    __global__ void get_gpu_pat_ipfx(int *pat_mat, int nth, int mem_size) {
        extern __shared__ int smarr[];
        __shared__ int clyr;
        __shared__ int ctrk, enc_ofst, rem_val, thload;
        __shared__ int pfx_ofst, tpfx_ofst;
        if(threadIdx.x >= nth) {
            return;
        }
        if(threadIdx.x == 0) {
            clyr = 0; ctrk = 0; pfx_ofst = 0; tpfx_ofst = 0;
            for(int i = 0; i < DZ && ctrk <= blockIdx.x; i++) {
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
                smarr[i] = pat_mat[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //get the prefix sum for each partition parallely
            int pfx=smarr[max(0, min(threadIdx.x * thload - 1, mem_size-1))], prev;
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
                tpfx_ofst = pat_mat[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)];
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
                pat_mat[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)] += smarr[i] + pfx_ofst;
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
    __global__ void update_cap_cube(int max_size) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        if((d_pre_demand_data[tid] > 0)) {
            // d_cube_data[tid] = max(float(-200), d_cube_data[tid] - 1);
            d_cube_data[tid] -= 1.0;
            #if EVAL_OFCOST == 1
            d_eval_wire_cnt[tid]++;
            #endif
        }
    }
    void gpu_commit_edges(int &max_nodes_len) {
        //try using a different approach for commiting vias, collect them when adding wires itself!        
        if(BIT_PMODEL) {
            MBC::bit_commit_edges(max_nodes_len);
        } else {
            CUDA_MEMSET(d_pre_demand_data, 0, CSIZE, int);
            #if EVAL_OFCOST == 1 && EVAL_CONN_CHECK == 1
            CUDA_MEMSET(d_eval_vis_data, 0, CSIZE, int);
            #endif
            int num_blks = CEIL_DIV(max_nodes_len, NUM_GPU_THREADS);
            gpu_add_wires<<<num_blks, NUM_GPU_THREADS>>>(max_nodes_len, MSC::BOFST, MSC::BLEN);
            cudaDeviceSynchronize();
            
            int mem_size = max(X, Y), nth = 128;
            UTIL::get_gpu_pfx_sum_blks(num_blks);
            get_gpu_pat_ipfx<<<num_blks, NUM_GPU_THREADS, mem_size * sizeof(int)>>>(d_pre_demand_data, nth, mem_size);
            cudaDeviceSynchronize();
            num_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
            update_cap_cube<<<num_blks, NUM_GPU_THREADS>>>(CSIZE);
            cudaDeviceSynchronize();
            num_blks = CEIL_DIV(max_nodes_len, NUM_GPU_THREADS);
            gpu_add_vias<<<num_blks, NUM_GPU_THREADS>>>(max_nodes_len, MSC::BOFST, MSC::BLEN);
            cudaDeviceSynchronize();
            #if EVAL_CONN_CHECK == 1 && EVAL_OFCOST == 1
            vector<int> eval_vis_data; eval_vis_data.clear(); eval_vis_data.resize(CSIZE);
            CUDAMEMCPYDH_1D(eval_vis_data.data(), d_eval_vis_data, CSIZE, int);
            queue<long> qu;
            while(qu.size()){qu.pop();}
            for(int i = 0; i < MSC::BLEN; i++) {
                int nid = npart_ids[i + MSC::BOFST];
                net_t cnet = nets[nid];
                #if MST_ASSERT == 1
                SPDLOG_ASSERT(cnet.pins.size() > 0, "{}", cnet.pins.size());
                #endif
                for(pin_t &cpin : cnet.pins) {
                    gvia_dt cz = cpin.pcd.z;
                    gcd_dt cx = cpin.pcd.x, cy = cpin.pcd.y;
                    int cenc = UTIL::get_enc<WIRE_ENC>(cz, cx, cy);
                    #if MST_ASSERT == 1
                    SPDLOG_ASSERT(eval_vis_data[cenc] == i + 1, "{} {} [batchid,net_name,z,x,y] [{},{},{},{},{}]", eval_vis_data[cenc], i+1, MSC::BT_ID, cnet.net_name, cz, cx, cy);
                    #endif
                }
                long qnode = UTIL::get_enc<WIRE_ENC>(cnet.pins[0].pcd.z, cnet.pins[0].pcd.x, cnet.pins[0].pcd.y);
                qnode = ((qnode << MST_NETID_BITS) | (i + 1));
                qu.push(qnode);
            }
            while(qu.size()) {
                long cq_node = qu.front();
                qu.pop();
                int net_id = (cq_node & MST_NETID_MASK);
                if((net_id >> MST_NETID_BITS) == 1)
                    continue;
                genc_dt enc = (cq_node >> MST_NETID_BITS);            
                gvia_dt z;
                gcd_dt x, y;
                UTIL::get_dec_all<WIRE_ENC>(enc, z, x, y); 
                // spdlog::info("new cds [zxy] [{},{},{}]", z, x, y);
                eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, x, y)] |= (1 << EVAL_CONN_VIS_BITS);
                if(z & 1) {
                    int ny = y + 1;
                    if(ny < Y && eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, x, ny)] == net_id) {
                        eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, x, ny)] |= (1 << EVAL_CONN_VIS_BITS);
                        long nq_node = UTIL::get_enc<WIRE_ENC>(z, x, ny);
                        nq_node = (nq_node << MST_NETID_BITS) | net_id;
                        qu.push(nq_node);
                    }
                    ny = y - 1;
                    if(ny >= 0 && eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, x, ny)] == net_id) {
                        eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, x, ny)] |= (1 << EVAL_CONN_VIS_BITS);
                        long nq_node = UTIL::get_enc<WIRE_ENC>(z, x, ny);
                        nq_node = (nq_node << MST_NETID_BITS) | net_id;
                        qu.push(nq_node);
                    }
                } else {
                    int nx = x + 1;
                    if(nx < X && eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, nx, y)] == net_id) {
                        eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, nx, y)] |= (1 << EVAL_CONN_VIS_BITS);
                        long nq_node = UTIL::get_enc<WIRE_ENC>(z, nx, y);
                        nq_node = (nq_node << MST_NETID_BITS) | net_id;
                        qu.push(nq_node);
                    }
                    nx = x - 1;
                    if(nx >= 0 && eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, nx, y)] == net_id) {
                        eval_vis_data[UTIL::get_enc<WIRE_ENC>(z, nx, y)] |= (1 << EVAL_CONN_VIS_BITS);
                        long nq_node = UTIL::get_enc<WIRE_ENC>(z, nx, y);
                        nq_node = (nq_node << MST_NETID_BITS) | net_id;
                        qu.push(nq_node);
                    }
                }
                int nz = z + 1;
                if(nz < Z && eval_vis_data[UTIL::get_enc<WIRE_ENC>(nz, x, y)] == net_id) {
                    eval_vis_data[UTIL::get_enc<WIRE_ENC>(nz, x, y)] |= (1 << EVAL_CONN_VIS_BITS);
                    long nq_node = UTIL::get_enc<WIRE_ENC>(nz, x, y);
                    nq_node = (nq_node << MST_NETID_BITS) | net_id;
                    qu.push(nq_node);
                }
                nz = z - 1;
                if(nz >= 0 && eval_vis_data[UTIL::get_enc<WIRE_ENC>(nz, x, y)] == net_id) {
                    eval_vis_data[UTIL::get_enc<WIRE_ENC>(nz, x, y)] |= (1 << EVAL_CONN_VIS_BITS);
                    long nq_node = UTIL::get_enc<WIRE_ENC>(nz, x, y);
                    nq_node = (nq_node << MST_NETID_BITS) | net_id;
                    qu.push(nq_node);
                }
            }
            for(int i = 0; i < MSC::BLEN; i++) {
                int nid = npart_ids[i + MSC::BOFST];
                net_t cnet = nets[nid];
                #if MST_ASSERT == 1
                SPDLOG_ASSERT(cnet.pins.size() > 0, "{}", cnet.pins.size());
                #endif
                for(pin_t &cpin : cnet.pins) {
                    gvia_dt cz = cpin.pcd.z;
                    gcd_dt cx = cpin.pcd.x, cy = cpin.pcd.y;
                    int cenc = UTIL::get_enc<WIRE_ENC>(cz, cx, cy);
                    #if MST_ASSERT == 1
                    // if(eval_vis_data[cenc] != ((i + 1) | (1 << EVAL_CONN_VIS_BITS))) {
                    //     spdlog::critical("failing [zxy] [{},{},{}]", cz,cx,cy);
                    // }
                    SPDLOG_ASSERT(eval_vis_data[cenc] == ((i + 1) | (1 << EVAL_CONN_VIS_BITS)), "{} {} [batchid,net_name,z,x,y] [{},{},{},{},{}]", eval_vis_data[cenc], i+1, MSC::BT_ID, cnet.net_name, cz, cx, cy);
                    #endif
                }
            }
            #endif
        }
    }
}
