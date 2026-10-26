#pragma once

#include "mst_common.h"
#include "../l_pattern_files/lpattern_common.h" 
using namespace std;
namespace MSC = MST_COMMON;

namespace GMST_UTIL {
    __device__ void atomicMin_double(gcb_dt* address, gcb_dt &val) {
        #if WT_FLOAT == 0
        unsigned long long int* address_as_ull = (unsigned long long int*)address;
        unsigned long long int old = *address_as_ull, assumed;
        do {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed, __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
        } while (assumed != old);
        #else
        int *address_as_ull = (int *)address;
        int old = *address_as_ull, assumed;
        do {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed, __float_as_int(fmin(val, __int_as_float(assumed))));
        } while (assumed != old);
        #endif
    }
    __device__ void atomicAdd_long(long* address, long val) {
        unsigned long long* address_as_ull = (unsigned long long*)address;
        unsigned long long old_ull;
        unsigned long long assumed_ull;
        do {
            old_ull = *address_as_ull;
            assumed_ull = old_ull;
            old_ull = (unsigned long long)((long long)old_ull + val);
        } while (atomicCAS(address_as_ull, assumed_ull, old_ull) != assumed_ull);
    }
    __global__ void get_gpu_net_block_map(long max_size, long max_edges) {
        long tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        long bofst = MSC::d_net_block_pfx[tid], blen = MSC::d_net_block_len[tid];
        atomicAdd_long(MSC::d_net_block_map + bofst, (long)(tid));
        if(bofst + blen < max_edges)
            atomicAdd_long(MSC::d_net_block_map + bofst + blen, (long)(-tid));
        // atomicAdd((unsigned long long*)MSC::d_net_block_map + bofst, (unsigned long long)(tid)); 
        // atomicAdd((unsigned long long*)MSC::d_net_block_map + bofst + blen, (unsigned long long)(-tid));
    }
    __global__ void get_gpu_net_node_map(int max_size, int max_nodes) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        int bofst = MSC::d_net_nodes_cnt_pfx[tid], blen = MSC::d_net_nodes_cnt[tid];
        atomicAdd(MSC::d_net_node_id + bofst, tid); 
        if(bofst + blen < max_nodes)
            atomicAdd(MSC::d_net_node_id + bofst + blen, (-tid));
    }
    __global__ void get_gpu_net_block_map1(int max_size) {
        long tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        long bofst = MSC::d_net_block_pfx1[tid], blen = MSC::d_net_block_len1[tid];
        atomicAdd_long(MSC::d_net_block_map1 + bofst, (long)(tid)); 
        if(bofst + blen < max_size)
            atomicAdd_long(MSC::d_net_block_map1 + bofst + blen, (long)(-tid));
    }
    __device__ int gdsu(int &dsu_ofst, int &id) {
        int tmp_id = id;
        while(tmp_id != MSC::d_dsu_data[dsu_ofst + tmp_id]) {
            tmp_id = MSC::d_dsu_data[dsu_ofst + tmp_id];
        }
        MSC::d_dsu_data[dsu_ofst + id] = tmp_id;
        // MSC::d_dsu_data[dsu_ofst + id] = tmp_id;
        // int tmp_id1 = id, prv_id = id;
        // while(tmp_id1 != MSC::d_dsu_data[dsu_ofst + tmp_id1]) {
        //     tmp_id1 = MSC::d_dsu_data[dsu_ofst + tmp_id1];
        //     MSC::d_dsu_data[dsu_ofst + prv_id] = tmp_id;
        //     prv_id = tmp_id1;
        // }
        return tmp_id;     
    }
    __global__ void get_gpu_mst(int max_nets, gcb_dt *lpat_res) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_nets)  {
            return;
        }
        int npfx = MSC::d_net_nodes_cnt_pfx[tid], nlen = MSC::d_net_nodes_cnt[tid];
        for(int i = npfx; i < npfx + nlen; i++) {
            MSC::d_dsu_data[i] = i - npfx;
        }
        long epfx = MSC::d_net_block_pfx1[tid], elen = MSC::d_net_block_len1[tid];
        int num_edges_picked = 0;
        MSC::d_net_dist[tid] = 0.0;
        for(long i = epfx; i < epfx + elen; i++) {
            long edge_val = MSC::d_net_edge_id[i];
            int n2 = (edge_val & NODE_MASK), n1 = ((edge_val >> NODE_BITS) & NODE_MASK);
            int cnid = (edge_val >> NNODE_BITS);
            #if MST_GPU_ASSERT == 1
            assert(cnid == tid);
            assert(n1 < n2);
            assert(n1 < nlen && n2 < nlen);
            #endif
            n1 = gdsu(npfx, n1);
            n2 = gdsu(npfx, n2);
            if(n1 != n2) {
                num_edges_picked++;
                if(n1 < n2) {
                    MSC::d_dsu_data[npfx + n1] = n2;
                } else {
                    MSC::d_dsu_data[npfx + n2] = n1;
                }
                MSC::d_net_dist[tid] += lpat_res[i];
            }
            if(num_edges_picked == nlen - 1) {
                break;
            }
        }
    }
    __global__ void get_gpu_edge_map(int bofst, int blen, long max_size) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        int onet_id = MSC::d_net_block_map[tid], nnet_id = d_npart_ids[bofst + onet_id], deg = d_net_deg_data[nnet_id];
        #if MST_GPU_ASSERT == 1
        assert(onet_id < blen);
        #endif
        long pfx = MSC::d_net_block_pfx[onet_id], tofst = tid - pfx;
        #if MST_GPU_ASSERT == 1
        // if(pfx > tid) {
        //     printf("pfx is %ld and tid is %ld\n", pfx, tid);
        // }
        assert(pfx <= tid);
        #endif
        long n1 = tofst / deg, n2 = tofst % deg;
        long pfx1 = MSC::d_net_block_pfx1[onet_id];
        if(n1 >= n2) {
            return;
        }
        #if NET_ID_ORDER_REV == 0
        long oft = (1L * n1 * deg) - ((n1 * (n1 + 1)) >> 1) + (n2 - n1 - 1);
        #else
        long oft = ((n2 * (n2 - 1)) >> 1) + n1;
        #endif
        #if MST_GPU_ASSERT == 1
        assert(oft >= 0 && onet_id >= 0);
        #endif
        long cid = (pfx1) + oft;
        MSC::d_net_edge_id[cid] = ((1L * onet_id) << NNODE_BITS) | (n1 << NODE_BITS) | n2;
    }
    __global__ void gpu_dsu_init(int max_len) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_len) {
            return;
        }
        int cnid = MSC::d_net_node_id[tid], npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        #if MST_GPU_ASSERT == 1
        assert(tid >= npfx);
        #endif
        int ofst = tid - npfx;
        MSC::d_dsu_data[npfx + ofst] = ofst;
    }
    __global__ void gpu_init_comp_val(long max_size) {
        long tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        MSC::d_brvsk_node_val[tid] = gcb_dt(MAX_OFGPU_COST);
    }
    template<int tp>
    __global__ void get_gpu_brvsk_min_edge(long max_size, gcb_dt *lpat_res) {
        long tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        long edge_val = MSC::d_net_edge_id[tid];
        #if EDGE_NODE_MAP_STORE == 1
        int nid = (edge_val >> NNODE_BITS), n1 = ((edge_val >> NODE_BITS) & NODE_MASK), n2 = (edge_val & NODE_MASK);
        #else
        int nid = (edge_val >> NNODE_BITS), ofid = tid - MSC::d_net_block_pfx1[nid];
        // if(tp == 1)
        // int n1 = ((edge_val >> NODE_BITS) & NODE_MASK), n2 = (edge_val & NODE_MASK);
        if(tp == 1) {
            if(!MSC::d_failed_net[nid]) {
                return;
            }
        }
        #if MST_GPU_ASSERT == 1
        assert(ofid >= 0);
        #endif
        int n1 = (sqrt((1.0 + (8.0 * ofid))));
        n1 = ((1 + n1) / 2);
        int n2 = ofid - ((n1 * (n1 - 1)) >> 1);
        #if MST_GPU_ASSERT == 1
        if(n2 < 0) {
            printf("[n1,n2,ofid] [%d,%d,%d]\n", n1, n2, ofid);
        }
        assert(n2 >= 0);
        #endif
        #endif
        int nofst = MSC::d_net_nodes_cnt_pfx[nid], nlen = MSC::d_net_nodes_cnt[nid];
        #if MST_GPU_ASSERT == 1
        assert(n1 < nlen);
        #endif
        n1 = gdsu(nofst, n1); n2 = gdsu(nofst, n2);
        if(n1 == n2) {
            return;
        }
        gcb_dt edge_dist = lpat_res[tid];
        atomicMin_double(MSC::d_brvsk_node_val + nofst + n1, edge_dist);
        atomicMin_double(MSC::d_brvsk_node_val + nofst + n2, edge_dist);
    }
    template<int tp>
    __global__ void get_gpu_brvsk_min_node(long max_size, gcb_dt *lpat_res) {
        long tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        long edge_val = MSC::d_net_edge_id[tid];
        #if EDGE_NODE_MAP_STORE == 1
        int nid = (edge_val >> NNODE_BITS), n1 = ((edge_val >> NODE_BITS) & NODE_MASK), n2 = (edge_val & NODE_MASK);
        #else
        int nid = (edge_val >> NNODE_BITS), ofid = tid - MSC::d_net_block_pfx1[nid];
        // int n1 = ((edge_val >> NODE_BITS) & NODE_MASK), n2 = (edge_val & NODE_MASK);
        #if MST_GPU_ASSERT == 1
        assert(ofid >= 0);
        #endif
        if(tp == 1) {
            if(!MSC::d_failed_net[nid]) {
                return;
            }
        }
        int n1 = (sqrt((1.0 + 8.0 * ofid)));
        n1 = ((1 + n1) / 2);
        int n2 = ofid - ((n1 * (n1 - 1)) >> 1);
        #if MST_GPU_ASSERT == 1
        assert(n2 >= 0);
        #endif
        #endif
        int nofst = MSC::d_net_nodes_cnt_pfx[nid], nlen = MSC::d_net_nodes_cnt[nid];
        #if MST_GPU_ASSERT == 1
        assert(n1 < nlen);
        #endif
        n1 = gdsu(nofst, n1); n2 = gdsu(nofst, n2);
        if(n1 == n2) {
            return;
        }
        gcb_dt edge_dist = lpat_res[tid];
        if(MSC::d_brvsk_node_val[nofst + n1] == edge_dist) {
            atomicExch(MSC::d_brvsk_node + nofst + n1, ((1 << BRVSK_VIS_BITS) | n2));
            // if(n1 < n2) {
            //     MSC::d_brvsk_node[nofst + n1] = ((1 << BRVSK_VIS_BITS) | n2);
            // }
            // MSC::d_brvsk_node[n1] = n2
        }
        if(MSC::d_brvsk_node_val[nofst + n2] == edge_dist) {
            atomicExch(MSC::d_brvsk_node + nofst + n2, ((1 << BRVSK_VIS_BITS) | n1));
            // if(n2 < n1) {
            //     MSC::d_brvsk_node[nofst + n2] = ((1 << BRVSK_VIS_BITS) | n1);
            // }
        }
    }
    template<int tp>
    __global__ void get_gpu_brvsk_edges(int max_size) {
        int tid = threadIdx.x + blockDim.x * blockIdx.x;
        if(tid >= max_size) {
            return;
        }
        int cnid = MSC::d_net_node_id[tid], npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        #if MST_GPU_ASSERT == 1
        assert(tid >= npfx);
        #endif
        if(tp == 1) {
            if(!MSC::d_failed_net[cnid]) {
                return;
            }
        }
        int ofst = tid - npfx;
        #if MST_GPU_ASSERT == 1
        assert(ofst < nlen);
        #endif
        if((MSC::d_brvsk_node[tid] >> BRVSK_VIS_BITS) == 1) {
            int next_nd = (MSC::d_brvsk_node[tid] & BRVSK_VIS_MASK);            
            #if MST_GPU_ASSERT == 1
            assert(next_nd < nlen);
            assert(next_nd != -1);
            assert(next_nd != ofst);
            #endif
            int brvsk_edid = atomicAdd(MSC::d_brvsk_net_len + cnid, 1);
            // MSC::d_brvsk_edges[brvsk_edid + npfx] = ((1L * ofst) << NODE_BITS) | next_nd;
            MSC::d_brvsk_edges[(brvsk_edid & NODE_MASK) + npfx] = tid;
            MSC::d_brvsk_node[tid] &= BRVSK_VIS_MASK;
        }
    }
    template<int tp>
    __global__ void get_gpu_brvsk_edges_conn(int max_size) {
        int tid = threadIdx.x + blockDim.x * blockIdx.x;
        if(tid >= max_size) {
            return;
        }
        int cnid = tid, npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        if(tp == 1) {
            if(!MSC::d_failed_net[cnid]) {
                return;
            }
        }
        int brvsk_edge_len = (MSC::d_brvsk_net_len[cnid] & NODE_MASK);
        int brvsk_edge_ofst = (MSC::d_brvsk_net_len[cnid] >> NODE_BITS);
        for(int i = npfx; i < npfx + brvsk_edge_len; i++) {
            long brvsk_edge_id = MSC::d_brvsk_edges[i];
            int on1 = brvsk_edge_id - npfx, on2 = MSC::d_brvsk_node[brvsk_edge_id];
            int n1 = gdsu(npfx, on1), n2 = gdsu(npfx, on2);
            if(n1 < n2) {
                MSC::d_dsu_data[npfx + n1] = n2;
                MSC::d_net_dist[cnid] += MSC::d_brvsk_node_val[brvsk_edge_id];                    
                #if TEST_BRVSK_MST_FORM == 1
                MSC::d_brvsk_edge_cnt[cnid] += 1;
                // MSC::d_net_dist[cnid] += 
                #endif
                #if ADD_GPU_MST_EDGES == 1
                MSC::d_net_mst_edges[npfx + brvsk_edge_ofst++] = (((1L * cnid) << NNODE_BITS) | ((1L * on1) << NODE_BITS) | (1L * on2));
                #endif
            } else if(n2 < n1) {
                MSC::d_dsu_data[npfx + n2] = n1;
                MSC::d_net_dist[cnid] += MSC::d_brvsk_node_val[brvsk_edge_id];
                #if TEST_BRVSK_MST_FORM == 1
                MSC::d_brvsk_edge_cnt[cnid] += 1;
                #endif
                #if ADD_GPU_MST_EDGES == 1
                MSC::d_net_mst_edges[npfx + brvsk_edge_ofst++] = (((1L * cnid) << NNODE_BITS) | ((1L * on1) << NODE_BITS) | (1L * on2));
                #endif
            }
        }
        MSC::d_brvsk_net_len[cnid] = ((brvsk_edge_ofst) << NODE_BITS);
        #if MST_GPU_ASSERT == 1
        assert(brvsk_edge_ofst < nlen);
        #endif
        MSC::d_net_mst_edges[npfx + brvsk_edge_ofst] = -1;
        // MSC::d_brvsk_net_len[cnid] = 0;
    }
    __global__ void get_gpu_failed_nets(int max_len) {
        int tid = threadIdx.x + blockDim.x * blockIdx.x;
        if(tid >= max_len) {
            return;
        }
        int cnid = tid, npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        int brvsk_edge_len = (MSC::d_brvsk_net_len[cnid] >> NODE_BITS);
        #if MST_GPU_ASSERT == 1
        assert(brvsk_edge_len < nlen);
        #endif
        if(brvsk_edge_len < nlen - 1) {
            MSC::d_failed_net[tid] = 1;
            MSC::d_is_nets_failed[0] = 1;
        }
    }
    __global__ void get_gpu_brvsk_conn1(int max_size) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        int cnid = MSC::d_net_node_id[tid], npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        #if MST_GPU_ASSERT == 1
        assert(tid >= npfx);
        #endif
        int ofst = tid - npfx;
        #if MST_GPU_ASSERT == 1
        assert(ofst < nlen);
        #endif
        if((MSC::d_brvsk_node[tid] >> BRVSK_VIS_BITS) == 1) {
            int next_nd = (MSC::d_brvsk_node[tid] & BRVSK_VIS_MASK);            
            #if MST_GPU_ASSERT == 1
            assert(next_nd < nlen);
            assert(next_nd != -1);
            assert(next_nd != ofst);
            assert(ofst < next_nd);
            int nx = gdsu(npfx, ofst), ny = gdsu(npfx, next_nd);
            assert(nx < ny);
            assert(MSC::d_dsu_data[npfx + ofst] == ofst);
            #endif
            // if(nx < ny) {
            //     atomicExch(MSC::d_dsu_data + npfx + nx, ny);
            //     atomicAdd(MSC::d_brvsk_edge_cnt + cnid, 1);
            // }
            MSC::d_brvsk_node[tid] &= BRVSK_VIS_MASK;         
        }
    }
    __global__ void get_gpu_brvsk_conn2(int max_size) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        int cnid = MSC::d_net_node_id[tid], npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        #if MST_GPU_ASSERT == 1
        assert(tid >= npfx);
        #endif
        int ofst = tid - npfx;
        #if MST_GPU_ASSERT == 1
        assert(ofst < nlen);
        #endif
        if((MSC::d_brvsk_node[tid] >> BRVSK_VIS_BITS) == 1) {
            int next_nd = (MSC::d_brvsk_node[tid] & BRVSK_VIS_MASK);            
            #if MST_GPU_ASSERT == 1
            assert(next_nd < nlen);
            assert(next_nd != -1);
            assert(next_nd != ofst);
            assert(next_nd < ofst);
            #endif
            if(ofst > next_nd) {
                int nx = gdsu(npfx, ofst), ny = gdsu(npfx, next_nd);
                if(nx != ny) {
                    MSC::d_dsu_data[npfx + ofst] = next_nd;
                    atomicAdd(MSC::d_brvsk_edge_cnt + cnid, 1);
                }
                MSC::d_brvsk_node[tid] &= BRVSK_VIS_MASK;
            }
        }
    }
    __global__ void get_gpu_brvsk_maxx_conn_node(int max_size) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        int cnid = MSC::d_net_node_id[tid], npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        #if MST_GPU_ASSERT == 1
        assert(tid >= npfx);
        #endif
        int ofst = tid - npfx;
        #if MST_GPU_ASSERT == 1
        assert(ofst < nlen);
        #endif
        if((MSC::d_brvsk_node[tid] >> BRVSK_VIS_BITS) == 1) {
            int next_nd = (MSC::d_brvsk_node[tid] & BRVSK_VIS_MASK);
            #if MST_GPU_ASSERT == 1
            assert(next_nd < nlen);
            assert(next_nd != -1);
            assert(next_nd != ofst);
            #endif
            if(ofst < next_nd) {
                atomicMax(MSC::d_dsu_data + npfx + ofst, next_nd);
            } else {
                atomicMax(MSC::d_dsu_data + npfx + next_nd, ofst);
            }
        }
    }
    __global__ void get_gpu_net_edge_cnt(int bofst, int blen) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= blen) {
            return;
        }
        int deg = d_net_deg_data[d_npart_ids[bofst + tid]];
        MSC::d_net_nodes_cnt[bofst + tid] = deg;
        MSC::d_net_block_len[bofst + tid] = deg * deg;
        MSC::d_net_block_len1[bofst + tid] = ((deg * (deg - 1)) >> 1);
    }
    __global__ void get_gpu_ofst(int bofst, int num_nodes_ofst) {
        if(threadIdx.x != 0) {
            return;
        }
        MSC::d_net_nodes_cnt = MSC::td_net_nodes_cnt + bofst;
        MSC::d_net_nodes_cnt_pfx = MSC::td_net_nodes_cnt_pfx + bofst;
        MSC::d_net_block_len = MSC::td_net_block_len + bofst;
        MSC::d_net_block_pfx = MSC::td_net_block_pfx + bofst;
        MSC::d_net_block_len1 = MSC::td_net_block_len1 + bofst;
        MSC::d_net_block_pfx1 = MSC::td_net_block_pfx1 + bofst;
        MSC::d_net_mst_edges = MSC::td_net_mst_edges + num_nodes_ofst;
        MSC::d_mst_lpat_pat = MSC::td_mst_lpat_pat + num_nodes_ofst;
    }
    __global__ void get_gpu_brvsk_conn(int max_size) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        int cnid = MSC::d_net_node_id[tid], npfx = MSC::d_net_nodes_cnt_pfx[cnid], nlen = MSC::d_net_nodes_cnt[cnid];
        #if MST_GPU_ASSERT == 1
        assert(tid >= npfx);
        #endif
        int ofst = tid - npfx;
        #if MST_GPU_ASSERT == 1
        assert(ofst < nlen);
        #endif
        if((MSC::d_brvsk_node[tid] >> BRVSK_VIS_BITS) == 1) {
            int next_nd = MSC::d_brvsk_node[tid] & BRVSK_VIS_MASK;
            #if TEST_BRVSK_MST_FORM == 1
            if(ofst > next_nd) {
                if(ofst != (MSC::d_brvsk_node[npfx + next_nd] & BRVSK_VIS_MASK)) {
                    // printf("greater [ofst,next_nd,tid,arr[next_nd]] [%d,%d,%d,%d]\n", ofst, next_nd, tid, (MSC::d_brvsk_node[npfx + next_nd] & BRVSK_VIS_MASK));
                    atomicAdd(MSC::d_brvsk_edge_cnt + cnid, 1);
                }
            } else {
                // printf("[ofst,next_nd,tid] [%d,%d,%d]\n", ofst, next_nd, tid);
                atomicAdd(MSC::d_brvsk_edge_cnt + cnid, 1);
            }
            #endif
            #if MST_GPU_ASSERT == 1
            assert(next_nd < nlen);
            assert(next_nd != -1);
            // assert(next_nd != ofst);
            #endif
            ofst = gdsu(npfx, ofst); next_nd = gdsu(npfx, next_nd);
            // #if MST_GPU_ASSERT == 1
            // assert(next_nd != ofst);
            // #endif
            if(ofst < next_nd) {
                // atomicMax(MSC::d_dsu_data)
                MSC::d_dsu_data[npfx + ofst] = next_nd;
            } else {
                MSC::d_dsu_data[npfx + next_nd] = ofst;
            }
            MSC::d_brvsk_node[tid] &= BRVSK_VIS_MASK;
        }
    }
    void gpu_brvsk_init(long &max_num_edges, int &max_num_nodes) {
        int num_blks = CEIL_DIV(max_num_nodes, NUM_GPU_THREADS);
        gpu_init_comp_val<<<num_blks, NUM_GPU_THREADS>>>(max_num_nodes);
        cudaDeviceSynchronize();
    }
    template<int tp>
    void gpu_brvsk_solve(long &max_num_edges, int &max_num_nodes) {
        gpu_brvsk_init(max_num_edges, max_num_nodes);
        // spdlog::info("calling brvsk");
        int num_blks = CEIL_DIV(max_num_edges, NUM_GPU_THREADS);
        get_gpu_brvsk_min_edge<tp><<<num_blks, NUM_GPU_THREADS>>>(max_num_edges, LPC::d_lpat_res);
        cudaDeviceSynchronize();
        get_gpu_brvsk_min_node<tp><<<num_blks, NUM_GPU_THREADS>>>(max_num_edges, LPC::d_lpat_res);
        cudaDeviceSynchronize();
        num_blks = CEIL_DIV(max_num_nodes, NUM_GPU_THREADS);
        get_gpu_brvsk_edges<tp><<<num_blks, NUM_GPU_THREADS>>>(max_num_nodes);
        cudaDeviceSynchronize();
        num_blks = CEIL_DIV(MSC::BLEN, NUM_GPU_THREADS);
        get_gpu_brvsk_edges_conn<tp><<<num_blks, NUM_GPU_THREADS>>>(MSC::BLEN);
        cudaDeviceSynchronize();
    }
}
