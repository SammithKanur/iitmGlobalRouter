#pragma once
#include "mst_common.h"
#include "gmst_util.hpp"
#include "../l_pattern_files/lpattern_common.h"
#include "./mst_commit.hpp"

using namespace std;
using namespace COMMON_DATA;

namespace GMSU = GMST_UTIL;
namespace MSC = MST_COMMON;
namespace CMT = MST_COMMIT;

namespace MST_UTIL {
    void get_edge_data(int tid) {
        // spdlog::info("getting edge data in thread id {}", tid);
        for(int i = tid; i < MSC::BLEN; i += MSC::MST_THREADS) {
            int cnet_id = npart_ids[MSC::BOFST + i];
            long deg = nets[cnet_id].pins.size();
            MSC::MAX_DEGREE = max(MSC::MAX_DEGREE, int(deg));
            #if MST_TREE_GPU == 0 || MST_CPU_GPU_COMPARISON == 1
            MSC::thr_dsu_size[tid] = max(MSC::thr_dsu_size[tid], (int)deg);
            MSC::net_edges_cnt[i] = (deg * (deg - 1)) >> 1;
            MSC::net_res_size[i] = deg - 1;
            #endif
            #if MST_GPU_CONN_CHECK == 1 || MST_CPU_GPU_COMPARISON == 1 || EDGES_CNT_GPU == 0
            MSC::net_nodes_cnt[i] = deg;
            #endif
        }
    }
    int gdsu(int &t_dsu_ofst, int id) {
        return (MSC::thr_dsu_space[t_dsu_ofst + id] == id) ? id : MSC::thr_dsu_space[t_dsu_ofst + id] = gdsu(t_dsu_ofst, MSC::thr_dsu_space[t_dsu_ofst + id]);
    }
    void t_get_mst(int tid) {
        int t_dsu_ofst = MSC::thr_dsu_pfx[tid], t_dsu_len = MSC::thr_dsu_size[tid];
        for(int nid = tid; nid < MSC::BLEN; nid += MSC::MST_THREADS) {
            MSC::net_dist[nid] = 0;
            long eofst = MSC::net_edges_cnt_pfx[nid], elen = MSC::net_edges_cnt[nid], cnt = 0;
            int deg = nets[npart_ids[nid + MSC::BOFST]].pins.size(), cnid = npart_ids[nid + MSC::BOFST];
            for(int i = 0; i < deg; i++) {
                #if NET_ID_ORDER_REV == 1
                for(int j = 0; j < i; j++) {
                #else
                for(int j = i+1; j < deg; j++) {
                #endif
                    ednd_dt cednd = ((1L * i) << MST_NODE_BITS) | j;
                    MSC::edge_arr[eofst + cnt] = {ed_id: eofst + cnt, ed_nds: cednd, dist: MSC::edge_dist[eofst + cnt]};
                    cnt++;
                }
            }
            sort(MSC::edge_arr + eofst, MSC::edge_arr + eofst + elen, [](MSC::edge_t &a, MSC::edge_t &b) {
                return a.dist < b.dist;
            });
            #if MST_ASSERT == 1
            SPDLOG_ASSERT(cnt == elen, "cnt {} elen {}", cnt, elen);
            SPDLOG_ASSERT(deg <= t_dsu_len, "deg {} t_dsu_len {}", deg, t_dsu_len);
            #endif
            //initialise dsu array
            for(int i = t_dsu_ofst; i < t_dsu_ofst + deg; i++) {
                MSC::thr_dsu_space[i] = i - t_dsu_ofst;
            }
            //start collecting edges
            int c_edg_cnt = 0, res_ofst = MSC::net_res_pfx[nid], res_len = MSC::net_res_size[nid];
            for(long i = eofst; i < eofst + elen && c_edg_cnt < deg - 1; i++) {
                MSC::edge_t cedge = MSC::edge_arr[i];
                long nodes = cedge.ed_nds;
                int nd1 = nodes & MST_NODE_MASK, nd2 = (nodes >> MST_NODE_BITS) & MST_NODE_MASK;
                nd1 = gdsu(t_dsu_ofst, nd1); nd2 = gdsu(t_dsu_ofst, nd2);
                if(nd1 != nd2) {
                    MSC::thr_dsu_space[t_dsu_ofst + nd1] = nd2;
                    MSC::net_dist[nid] += cedge.dist;
                    MSC::mst_edges[res_ofst + c_edg_cnt] = (cedge.ed_id << MST_NETID_BITS) | cnid;
                    #if MST_CONN_CHECK == 1
                    MSC::mst_conn_check_edges[res_ofst + c_edg_cnt] = cedge;
                    #endif
                    c_edg_cnt++;
                }
            }
            #if MST_ASSERT == 1
            SPDLOG_ASSERT(c_edg_cnt == res_len, "c_edg_len {} res_len {}", c_edg_cnt, res_len);
            #endif
            #if MST_CONN_CHECK == 1
            for(int i = t_dsu_ofst; i < t_dsu_ofst + deg; i++) {
                MSC::thr_dsu_space[i] = i - t_dsu_ofst;
            }
            for(int i = res_ofst; i < res_ofst + res_len; i++) {
                MSC::edge_t cedge = MSC::mst_conn_check_edges[i];
                long nodes = cedge.ed_nds;
                int nd1 = nodes & MST_NODE_MASK, nd2 = (nodes >> MST_NODE_BITS) & MST_NODE_MASK;
                nd1 = gdsu(t_dsu_ofst, nd1); nd2 = gdsu(t_dsu_ofst, nd2);
                if(nd1 != nd2) {
                    MSC::thr_dsu_space[t_dsu_ofst + nd1] = nd2;
                }
            }
            int par_nd = gdsu(t_dsu_ofst, 0);
            for(int i = 0; i < deg; i++) {
                int nxt_nd = gdsu(t_dsu_ofst, i);
                SPDLOG_ASSERT(par_nd == nxt_nd, "mst tree is not connected par_nd {} nxt_nd {}", par_nd, nxt_nd);
            }
            #endif
        }
    }
    struct compare_nnode_bits {
        __host__ __device__
        bool operator()(const long &a, const long &b) const {
            return (a >> NNODE_BITS) < (b >> NNODE_BITS);
        }
    };
    int tgdsu(vector<int> &dsu_data, int id) {
        return (dsu_data[id] == id) ? id : dsu_data[id] = tgdsu(dsu_data, dsu_data[id]);
    }
    void check_gpu_mst_conn(int &max_num_nodes) {
        // spdlog::info("testing gpu conn");
        #if (MST_GPU_CONN_CHECK == 1) && ADD_GPU_MST_EDGES == 1
        vector<nenc_dt> mst_edges_cpu; mst_edges_cpu.clear(); mst_edges_cpu.resize(max_num_nodes);
        CUDAMEMCPYDH_1D(mst_edges_cpu.data(), MSC::d_net_mst_edges, max_num_nodes, nenc_dt);        
        vector<int> brvsk_net_len; brvsk_net_len.clear(); brvsk_net_len.resize(MSC::BLEN, 0);
        CUDAMEMCPYDH_1D(brvsk_net_len.data(), MSC::d_brvsk_net_len, MSC::BLEN, int);
        for(int i = 0; i < MSC::BLEN; i++) {
            net_t cnet = nets[npart_ids[MSC::BOFST + i]];
            int deg = cnet.pins.size();
            long nlen = MSC::net_nodes_cnt[i], npfx = MSC::net_nodes_cnt_pfx[i];
            vector<int> cpu_dsu_data;
            cpu_dsu_data.clear(); cpu_dsu_data.resize(deg);
            for(int j = 0; j < deg; j++) {
                cpu_dsu_data[j] = j;
            }
            SPDLOG_ASSERT(deg == MSC::net_nodes_cnt[i], "{} {}", deg, MSC::net_nodes_cnt[i]);
            int num_edges = (brvsk_net_len[i] >> NODE_BITS);
            for(int j = npfx; j < npfx + num_edges; j++) {
                nenc_dt edg = mst_edges_cpu[j];
                int n1 = edg & NODE_MASK, n2 = ((edg >> NODE_BITS) & NODE_MASK);
                // spdlog::info("checking conn for [zxy] [{},{},{}] [{},{},{}]", cnet.pins[n1].pcd.z, cnet.pins[n1].pcd.x, cnet.pins[n1].pcd.y, cnet.pins[n2].pcd.z, cnet.pins[n2].pcd.x, cnet.pins[n2].pcd.y);
                SPDLOG_ASSERT(n1 < nlen && n2 < nlen, "{} {} {}", n1, n2, nlen);
                n1 = tgdsu(cpu_dsu_data, n1); n2 = tgdsu(cpu_dsu_data, n2);
                if(n1 != n2) {
                    if(n1 < n2) {
                        cpu_dsu_data[n1] = n2;
                    } else {
                        cpu_dsu_data[n2] = n1;
                    }
                }
            }
            int par = tgdsu(cpu_dsu_data, 0);
            for(int j = 0; j < nlen; j++) {
                int cval = tgdsu(cpu_dsu_data, j);
                SPDLOG_ASSERT(par == cval, "{} {}", par, cval);
            }
        }
        #endif
    }
    void get_block_net_map() {
        CUDA_MALLOC_1D(MSC::d_net_dist, MSC::BLEN, gcb_dt); CUDA_MEMSET(MSC::d_net_dist, 0, MSC::BLEN, gcb_dt);
        CUDA_MALLOC_1D(MSC::d_failed_net, MSC::BLEN, bool); CUDA_MEMSET(MSC::d_failed_net, 0, MSC::BLEN, bool);
        CUDA_MALLOC_1D(MSC::d_is_nets_failed, 1, bool); CUDA_MEMSET(MSC::d_is_nets_failed, 0, 1, bool);                
        CUDA_MALLOC_1D(MSC::d_net_edge_id, MSC::NUM_EDGES, long);
        long num_edges, num_blks;
        int max_nodes_len;
        GMSU::get_gpu_ofst<<<1, NUM_GPU_THREADS>>>(MSC::BOFST, MSC::NUM_NODES_OFST);
        cudaDeviceSynchronize();
        #if EDGES_CNT_GPU == 0
        vector<long> block_net_len, block_net_pfx; block_net_len.clear(); block_net_len.resize(MSC::BLEN); block_net_pfx.clear(); block_net_pfx.resize(MSC::BLEN);
        vector<long> block_net_len1, block_net_pfx1; block_net_len1.clear(); block_net_len1.resize(MSC::BLEN); block_net_pfx1.clear(); block_net_pfx1.resize(MSC::BLEN);
        for(int i = MSC::BOFST; i < MSC::BOFST + MSC::BLEN; i++) {
            long deg = nets[npart_ids[i]].pins.size();
            block_net_len[i - MSC::BOFST] = deg * deg;
            block_net_len1[i - MSC::BOFST] = (deg * (deg - 1)) / 2;
        }
        UTIL::get_vec_epfx<long>(block_net_pfx, block_net_len, MSC::BLEN);
        UTIL::get_vec_epfx<long>(block_net_pfx1, block_net_len1, MSC::BLEN);
        CUDAMEMCPYHD_1D(MSC::d_net_nodes_cnt, MSC::net_nodes_cnt.data(), MSC::BLEN, int);
        CUDAMEMCPYHD_1D(MSC::d_net_nodes_cnt_pfx, MSC::net_nodes_cnt_pfx.data(), MSC::BLEN, int);
        CUDAMEMCPYHD_1D(MSC::d_net_block_len, block_net_len.data(), MSC::BLEN, long);
        CUDAMEMCPYHD_1D(MSC::d_net_block_pfx, block_net_pfx.data(), MSC::BLEN, long);
        CUDAMEMCPYHD_1D(MSC::d_net_block_len1, block_net_len1.data(), MSC::BLEN, long);
        CUDAMEMCPYHD_1D(MSC::d_net_block_pfx1, block_net_pfx1.data(), MSC::BLEN, long);
        num_edges = block_net_pfx[MSC::BLEN - 1] + block_net_len[MSC::BLEN - 1];
        num_blks = CEIL_DIV(MSC::BLEN, NUM_GPU_THREADS);
        CUDA_MALLOC_1D(MSC::d_dsu_data, MSC::net_nodes_cnt.back() + MSC::net_nodes_cnt_pfx.back(), int);
        max_nodes_len = MSC::net_nodes_cnt.back() + MSC::net_nodes_cnt_pfx.back();
        #else
        num_blks = CEIL_DIV(MSC::BLEN, NUM_GPU_THREADS);
        // GMSU::get_gpu_net_edge_cnt<<<num_blks, NUM_GPU_THREADS>>>(MSC::BOFST, MSC::BLEN);
        // cudaDeviceSynchronize();
        thrust::device_ptr<long> d_edge_cnt_ptr1 = thrust::device_pointer_cast(MSC::d_net_block_len1);                
        thrust::device_ptr<long> d_edge_cnt_optr1 = thrust::device_pointer_cast(MSC::d_net_block_pfx1);
        thrust::exclusive_scan(d_edge_cnt_ptr1, d_edge_cnt_ptr1 + MSC::BLEN, d_edge_cnt_optr1);        
        cudaDeviceSynchronize();
        long max_edge_sum_pf, max_edge_sum;
        int nodes_cnt, nodes_pfx;
        CUDAMEMCPYDH_1D(&max_edge_sum_pf, MSC::d_net_block_pfx + MSC::BLEN - 1, 1, long);
        CUDAMEMCPYDH_1D(&max_edge_sum, MSC::d_net_block_len + MSC::BLEN - 1, 1, long);
        CUDAMEMCPYDH_1D(&nodes_cnt, MSC::d_net_nodes_cnt + MSC::BLEN - 1, 1, int);
        CUDAMEMCPYDH_1D(&nodes_pfx, MSC::d_net_nodes_cnt_pfx + MSC::BLEN - 1, 1, int);
        num_edges = max_edge_sum_pf + max_edge_sum;
        max_nodes_len = nodes_cnt + nodes_pfx;
        CUDA_MALLOC_1D(MSC::d_dsu_data, max_nodes_len, int);
        #if MST_MIGRATE_LPAT == 1
        // CUDA_MALLOC_1D(MSC::d_mst_lpat_pat, max_nodes_len, gvia_dt);
        #endif
        #endif
        CUDA_MALLOC_1D(MSC::d_net_block_map, num_edges, long); CUDA_MEMSET(MSC::d_net_block_map, 0, num_edges, long);
        GMSU::get_gpu_net_block_map<<<num_blks, NUM_GPU_THREADS>>>(MSC::BLEN, num_edges);
        cudaDeviceSynchronize();
        thrust::device_ptr<long> d_ptr = thrust::device_pointer_cast(MSC::d_net_block_map);
        thrust::inclusive_scan(d_ptr, d_ptr + num_edges, d_ptr);
        cudaDeviceSynchronize();
        #if MST_BRVSK == 1
        CUDA_MALLOC_1D(MSC::d_net_node_id, max_nodes_len, int); CUDA_MEMSET(MSC::d_net_node_id, 0, max_nodes_len, int);
        GMSU::get_gpu_net_node_map<<<num_blks, NUM_GPU_THREADS>>>(MSC::BLEN, max_nodes_len);
        cudaDeviceSynchronize();
        thrust::device_ptr<int> d_node_ptr = thrust::device_pointer_cast(MSC::d_net_node_id);
        thrust::inclusive_scan(d_node_ptr, d_node_ptr + max_nodes_len, d_node_ptr);
        cudaDeviceSynchronize();
        num_blks = CEIL_DIV(max_nodes_len, NUM_GPU_THREADS);
        GMSU::gpu_dsu_init<<<num_blks, NUM_GPU_THREADS>>>(max_nodes_len);
        cudaDeviceSynchronize();
        #if TEST_MST_GPU_DSU_INIT == 1
        vector<int> dsu_cpu; dsu_cpu.clear(); dsu_cpu.resize(max_nodes_len);
        CUDAMEMCPYDH_1D(dsu_cpu.data(), MSC::d_dsu_data, max_nodes_len, int);
        int dsu_id = 0;
        for(int i = 0; i < MSC::BLEN; i++) {
            int deg = nets[npart_ids[MSC::BOFST + i]].pins.size();
            for(int j = 0; j < deg; j++) {
                int dsu_val = dsu_cpu[dsu_id];
                SPDLOG_ASSERT(dsu_val == j, "{} {}", dsu_val, j);
                dsu_id++;
            }
        }
        #endif
        #endif
        #if TEST_MST_GPU_MAP == 1
        vector<long> net_block_map_cpu; net_block_map_cpu.clear(); net_block_map_cpu.resize(num_edges);
        CUDAMEMCPYDH_1D(net_block_map_cpu.data(), MSC::d_net_block_map, num_edges, long);
        long eid = 0;
        for(int i = 0; i < MSC::BLEN; i++) {
            int cnid = npart_ids[MSC::BOFST + i];
            long deg = nets[cnid].pins.size();
            for(int j = 0; j < deg * deg; j++) {
                SPDLOG_ASSERT(net_block_map_cpu[eid] == i, "{} {}", i, net_block_map_cpu[eid]);
                eid++;
            }
        }
        net_block_map_cpu.clear();
        #endif
        num_blks = CEIL_DIV(num_edges, NUM_GPU_THREADS);
        GMSU::get_gpu_edge_map<<<num_blks, NUM_GPU_THREADS>>>(MSC::BOFST, MSC::BLEN, num_edges);
        cudaDeviceSynchronize();
        #if TEST_MST_GPU_MAP == 1
        vector<long> edge_map_cpu; edge_map_cpu.clear(); edge_map_cpu.resize(MSC::NUM_EDGES);
        CUDAMEMCPYDH_1D(edge_map_cpu.data(), MSC::d_net_edge_id, MSC::NUM_EDGES, long);
        long teid = 0;
        for(int i = 0; i < MSC::BLEN; i++) {
            int deg = nets[npart_ids[MSC::BOFST + i]].pins.size();
            for(int j = 0; j < deg; j++) {
                for(int k = j + 1; k < deg; k++) {
                    long edg = edge_map_cpu[teid++];
                    int n1 = edg & NODE_MASK, n2 = ((edg >> NODE_BITS) & NODE_MASK);
                    int nid = (edg >> NNODE_BITS);
                    swap(n1, n2);
                    SPDLOG_ASSERT(n1 == j && n2 == k && nid == i, "{} {} {} {} {} {}", n1, j, n2, k, nid, i);
                }
            }
        }
        #endif
        #if MST_TREE_GPU == 1
        #if MST_GPU_ZIP_SORT == 0
        // thrust::device_ptr node_ptr = thrust::device_pointer_cast(MSC::d_net_edge_id);
        // thrust::device_ptr edg_ptr = thrust::device_pointer_cast(LPC::d_lpat_res);
        // thrust::sort_by_key(edg_ptr, edg_ptr + MSC::NUM_EDGES, node_ptr);
        // thrust::sort_by_key(node_ptr, node_ptr + MSC::NUM_EDGES, edg_ptr, compare_nnode_bits());
        cudaDeviceSynchronize();
        #endif
        #if MST_BRVSK == 0
        num_blks = CEIL_DIV(MSC::BLEN, NUM_GPU_THREADS);
        GMSU::get_gpu_mst<<<num_blks, NUM_GPU_THREADS>>>(MSC::BLEN, LPC::d_lpat_res);
        cudaDeviceSynchronize();
        #else
        CUDA_MALLOC_1D(MSC::d_brvsk_node_val, max_nodes_len, gcb_dt);
        CUDA_MALLOC_1D(MSC::d_brvsk_edges, max_nodes_len, long);
        #if ADD_GPU_MST_EDGES == 1
        // CUDA_MALLOC_1D(MSC::d_net_mst_edges, max_nodes_len, nenc_dt);
        #endif
        CUDA_MALLOC_1D(MSC::d_brvsk_net_len, MSC::BLEN, int); CUDA_MEMSET(MSC::d_brvsk_net_len, 0, MSC::BLEN, int);
        CUDA_MALLOC_1D(MSC::d_brvsk_node, max_nodes_len, int); CUDA_MEMSET(MSC::d_brvsk_node, -1, max_nodes_len, int);
        #if TEST_BRVSK_MST_FORM == 1
        CUDA_MALLOC_1D(MSC::d_brvsk_edge_cnt, MSC::BLEN, int); CUDA_MEMSET(MSC::d_brvsk_edge_cnt, 0, MSC::BLEN, int);
        #endif
        GMSU::gpu_brvsk_init(MSC::NUM_EDGES, max_nodes_len);
        int num_itrs = ceil(log2(1.0 * MSC::MAX_DEGREE));
        #if MST_FEEDBACK == 0
        for(int i = 0; i < 2 * num_itrs; i++) {
        #else
        for(int i = 0; i < num_itrs; i++) {
        #endif
            GMSU::gpu_brvsk_solve<0>(MSC::NUM_EDGES, max_nodes_len);
        }
        #if MST_FEEDBACK == 1
        //try spawning only necessary amount of threads for this, instead of returning un-necessary threads
        num_blks = CEIL_DIV(MSC::BLEN, NUM_GPU_THREADS);
        GMSU::get_gpu_failed_nets<<<num_blks, NUM_GPU_THREADS>>>(MSC::BLEN);
        bool is_nets_failed = false;
        CUDAMEMCPYDH_1D(&is_nets_failed, MSC::d_is_nets_failed, 1, bool);
        while(is_nets_failed) {
            spdlog::warn("running brvsk once again for failed nets for batch {}", MSC::BT_ID);
            for(int i = 0; i < num_itrs; i++) {
                GMSU::gpu_brvsk_solve<1>(MSC::NUM_EDGES, max_nodes_len);
            }
            CUDA_MEMSET(MSC::d_failed_net, 0, MSC::BLEN, bool); CUDA_MEMSET(MSC::d_is_nets_failed, 0, 1, bool);
            GMSU::get_gpu_failed_nets<<<num_blks, NUM_GPU_THREADS>>>(MSC::BLEN);
            is_nets_failed = false;
            CUDAMEMCPYDH_1D(&is_nets_failed, MSC::d_is_nets_failed, 1, bool);
        }
        #endif
        #if MST_COMMIT_EDGES == 1
        CMT::gpu_commit_edges(max_nodes_len);
        CUDAMEMCPYDH_1D(MSC::net_mst_edges.data() + MSC::NUM_NODES_OFST, MSC::td_net_mst_edges + MSC::NUM_NODES_OFST, max_nodes_len, nenc_dt);
        CUDAMEMCPYDH_1D(MSC::mst_lpat_pat.data() + MSC::NUM_NODES_OFST, MSC::td_mst_lpat_pat + MSC::NUM_NODES_OFST, max_nodes_len, gvia_dt);
        write_route_bsem.release();
        #endif
        #if TEST_BRVSK_MST_FORM == 1
        vector<int> brvsk_edge_cnt; brvsk_edge_cnt.clear(); brvsk_edge_cnt.resize(MSC::BLEN, 0);
        // vector<int> brvsk_net_len; brvsk_net_len.clear(); brvsk_net_len.resize(MSC::BLEN, 0);
        CUDAMEMCPYDH_1D(brvsk_edge_cnt.data(), MSC::d_brvsk_edge_cnt, MSC::BLEN, int);
        // CUDAMEMCPYDH_1D(brvsk_net_len.data(), MSC::d_brvsk_net_len, MSC::BLEN, int);
        // for(int i = 0; i < MSC::BLEN; i++) {
        //     SPDLOG_ASSERT(brvsk_net_len[i] > 0, "{}", brvsk_net_len[i]);
        // }
        for(int i = 0; i < MSC::BLEN; i++) {
            net_t cnet = nets[npart_ids[MSC::BOFST + i]];
            int deg = cnet.pins.size();
            if(brvsk_edge_cnt[i] != deg - 1) {
                spdlog::critical("rerun mst for net [id,name,brvsk_edge_cnt,deg] [{},{},{},{}]", i, cnet.net_name, brvsk_edge_cnt[i], deg);
            }
            // SPDLOG_ASSERT(brvsk_edge_cnt[i] == deg - 1 || brvsk_edge_cnt[i] == deg, "{} {} {} {}", i, brvsk_edge_cnt[i], deg, cnet.net_name);
        }
        #endif
        // GMSU::gpu_brvsk_solve(MSC::NUM_EDGES, max_nodes_len);
        // num_blks = CEIL_DIV(MSC::BLEN, NUM_GPU_THREADS);
        // GMSU::get_gpu_mst<<<num_blks, NUM_GPU_THREADS>>>(MSC::BLEN, LPC::d_lpat_res);
        // cudaDeviceSynchronize();
        #if MST_GPU_CONN_CHECK == 1 && ADD_GPU_MST_EDGES == 1
        check_gpu_mst_conn(max_nodes_len);
        #endif
        #endif
        #endif
        #if TEST_MST_GPU_SORT == 1
        vector<long> keys; keys.clear(), keys.resize(MSC::NUM_EDGES);
        vector<gcb_dt> vals; vals.clear(); vals.resize(MSC::NUM_EDGES);
        CUDAMEMCPYDH_1D(keys.data(), MSC::d_net_edge_id, MSC::NUM_EDGES, long);
        CUDAMEMCPYDH_1D(vals.data(), LPC::d_lpat_res, MSC::NUM_EDGES, gcb_dt);
        int teid1 = 0;
        for(int i = 0; i < MSC::BLEN; i++) {
            int deg = nets[npart_ids[MSC::BOFST + i]].pins.size();
            int max_edges = (deg * (deg - 1)) >> 1;
            for(int j = 0; j < max_edges-1; j++) {
                long edg1 = keys[teid1], edg2 = keys[teid1+1];
                gcb_dt val1 = vals[teid1], val2 = vals[teid1+1];
                SPDLOG_ASSERT(val1 <= val2, "{} {}", val1, val2);
                int nid1 = (edg1 >> NNODE_BITS), nid2 = (edg2 >> NNODE_BITS);
                SPDLOG_ASSERT(nid1 == nid2, "{} {}", nid1, nid2);
                    // int n1 = edg & NODE_MASK, n2 = ((edg >> NODE_BITS) & NODE_MASK);
                    // int nid = (edg >> NNODE_BITS);
                    // swap(n1, n2);
                    // SPDLOG_ASSERT(n1 == j && n2 == k && nid == i, "{} {} {} {} {} {}", n1, j, n2, k, nid, i);
                teid1++;
            }
            teid1++;
        }
        #endif
        MSC::NUM_NODES_OFST += max_nodes_len;
    }
    void get_mst() {
        MSC::MAX_DEGREE = 0;
        // spdlog::info("doing for a batch");
        #if MST_TREE_GPU == 0 || MST_CPU_GPU_COMPARISON == 1
        MSC::net_edges_cnt.clear(); MSC::net_edges_cnt.resize(MSC::BLEN);
        MSC::net_edges_cnt_pfx.clear(); MSC::net_edges_cnt_pfx.resize(MSC::BLEN);
        MSC::thr_dsu_pfx.clear(); MSC::thr_dsu_pfx.resize(MSC::MST_THREADS, 0);
        MSC::thr_dsu_size.clear(); MSC::thr_dsu_size.resize(MSC::MST_THREADS, 0);
        MSC::net_res_size.clear(); MSC::net_res_size.resize(MSC::BLEN, 0);
        MSC::net_res_pfx.clear(); MSC::net_res_pfx.resize(MSC::BLEN);
        MSC::net_dist.clear(); MSC::net_dist.resize(MSC::BLEN);
        #endif
        #if MST_GPU_CONN_CHECK == 1 || MST_CPU_GPU_COMPARISON == 1 || EDGES_CNT_GPU == 0 
        MSC::net_nodes_cnt.clear(); MSC::net_nodes_cnt.resize(MSC::BLEN);
        MSC::net_nodes_cnt_pfx.clear(); MSC::net_nodes_cnt_pfx.resize(MSC::BLEN);
        #endif
        for(int i = 0; i < MSC::MST_THREADS; i++) {
            MSC::thrds[i] = thread(get_edge_data, i);
        }
        for(int i = 0; i < MSC::MST_THREADS; i++) {
            MSC::thrds[i].join();
        }
        #if MST_TREE_GPU == 0 || MST_CPU_GPU_COMPARISON == 1
        UTIL::get_vec_epfx<long>(MSC::net_edges_cnt_pfx, MSC::net_edges_cnt, MSC::BLEN);
        UTIL::get_vec_epfx<int>(MSC::thr_dsu_pfx, MSC::thr_dsu_size, MSC::MST_THREADS);
        UTIL::get_vec_epfx<int>(MSC::net_res_pfx, MSC::net_res_size, MSC::BLEN);
        #endif
        #if MST_GPU_CONN_CHECK == 1 || MST_CPU_GPU_COMPARISON == 1 || EDGES_CNT_GPU == 0
        UTIL::get_vec_epfx<int>(MSC::net_nodes_cnt_pfx, MSC::net_nodes_cnt, MSC::BLEN);
        #endif
        #if MST_TREE_GPU == 1
        // spdlog::info("number of net nodes {}", MSC::net_nodes_cnt_pfx.back() + MSC::net_nodes_cnt.back());
        get_block_net_map();
        #endif
        #if MST_TREE_GPU == 0 || MST_CPU_GPU_COMPARISON == 1
        MSC::thr_dsu_space.clear(); MSC::thr_dsu_space.resize(MSC::thr_dsu_size[MSC::MST_THREADS - 1] + MSC::thr_dsu_pfx[MSC::MST_THREADS - 1]);
        int num_mst_edges = MSC::net_res_pfx[MSC::BLEN - 1] + MSC::net_res_size[MSC::BLEN - 1];
        MSC::mst_edges.clear(); MSC::mst_edges.resize(num_mst_edges);
        #endif
        #if MST_CONN_CHECK == 1 && MST_TREE_GPU == 0
        MSC::mst_conn_check_edges.clear(); MSC::mst_conn_check_edges.resize(num_mst_edges);
        #endif        
        #if MST_ASSERT == 1 && MST_TREE_GPU == 0
        num_mst_edges = MSC::net_edges_cnt_pfx[MSC::BLEN - 1] + MSC::net_edges_cnt[MSC::BLEN-1];
        SPDLOG_ASSERT(num_mst_edges == MSC::NUM_EDGES, "num_mst_edges {} MSC::NUM_EDGES {}", num_mst_edges, MSC::NUM_EDGES);
        #endif
        #if MST_TREE_GPU == 0 || MST_CPU_GPU_COMPARISON == 1
        MSC::edge_arr = (MSC::edge_t *)malloc(sizeof(MSC::edge_t) * MSC::NUM_EDGES);
        //sort the edge arr for edges belonging to a net
        for(int i = 0; i < MSC::MST_THREADS; i++) {
            MSC::thrds[i] = thread(t_get_mst, i);
        }
        for(int i = 0; i < MSC::MST_THREADS; i++) {
            MSC::thrds[i].join();
        }
        free(MSC::edge_arr);
        #endif
        #if MST_CPU_GPU_COMPARISON == 1
        vector<gcb_dt> net_dist_cpu; net_dist_cpu.clear(); net_dist_cpu.resize(MSC::BLEN);
        CUDAMEMCPYDH_1D(net_dist_cpu.data(), MSC::d_net_dist, MSC::BLEN, gcb_dt);
        // spdlog::warn("testing brvsk and mst cpu");
        for(int i = 0; i < MSC::BLEN; i++) {
            net_t cnet = nets[npart_ids[MSC::BOFST + i]];
            int deg = cnet.pins.size();
            #if WT_FLOAT == 1
            if(abs(MSC::net_dist[i] - net_dist_cpu[i]) > FLOAT_ACC_VAL) {
            #else
            if(abs(MSC::net_dist[i] - net_dist_cpu[i]) > DOUBLE_ACC_VAL) {
            #endif
                spdlog::critical("rerun brvsk mst build for net [id,name,cpu_mst_res,brvsk_mst_res] [{},{},{},{}]", i, cnet.net_name, MSC::net_dist[i], net_dist_cpu[i]);
            }
            #if WT_FLOAT == 0
            SPDLOG_ASSERT(MSC::net_dist[i] == net_dist_cpu[i], "{} {}", MSC::net_dist[i], net_dist_cpu[i]);
            #endif
        }
        #endif
        #if MST_TREE_GPU == 1
        cudaFree(MSC::d_net_dist); cudaFree(MSC::d_net_edge_id); cudaFree(MSC::d_dsu_data); 
        cudaFree(MSC::d_net_block_map);
        #if MST_MIGRATE_LPAT == 1
        // cudaFree(MSC::d_mst_lpat_pat);
        #endif
        #if MST_BRVSK == 1
        cudaFree(MSC::d_net_node_id); cudaFree(MSC::d_brvsk_node_val); cudaFree(MSC::d_brvsk_node);
        cudaFree(MSC::d_brvsk_net_len); cudaFree(MSC::d_brvsk_edges);
        #if TEST_BRVSK_MST_FORM == 1
        cudaFree(MSC::d_brvsk_edge_cnt); 
        #endif
        #if ADD_GPU_MST_EDGES == 1
        // cudaFree(MSC::d_net_mst_edges);
        #endif
        #endif
        #endif
        //now do mst for every net, as the edges are sorted, use krsk algorithm
    }
}
