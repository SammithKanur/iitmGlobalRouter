#pragma once
#include "brvsk_common.h"

using namespace std;
using namespace COMMON_DATA;

namespace BRC = BRVSK_COMMON;

namespace BRVSK_UTIL {
    void get_edge_data(int tid) {
        // spdlog::info("getting edge data in thread id {}", tid);
        for(int i = tid; i < BRC::BLEN; i += BRC::NUM_THREADS) {
            int cnet_id = npart_ids[BRC::BOFST + i];
            long deg = nets[cnet_id].pins.size();
            BRC::thr_dsu_size[tid] = max(BRC::thr_dsu_size[tid], (int)deg);
            // BRC::net_nodes_cnt[i] = deg;
            BRC::net_edges_cnt[i] = (deg * (deg - 1)) >> 1;
            BRC::net_res_size[i] = deg - 1;
        }
    }
    int gdsu(int &t_dsu_ofst, int &id) {
        return (BRC::thr_dsu_space[t_dsu_ofst + id] == id) ? id : BRC::thr_dsu_space[t_dsu_ofst + id] = gdsu(t_dsu_ofst, BRC::thr_dsu_space[t_dsu_ofst + id]);
    }
    void t_get_mst(int tid) {
        int t_dsu_ofst = BRC::thr_dsu_pfx[tid], t_dsu_len = BRC::thr_dsu_size[tid];
        int t_nofst = t_dsu_ofst, t_nlen = t_dsu_len;
        for(int nid = tid; nid < BRC::BLEN; nid += BRC::NUM_THREADS) {
            BRC::net_dist[nid] = 0;
            long eofst = BRC::net_edges_cnt_pfx[nid], elen = BRC::net_edges_cnt[nid], cnt = 0;
            int deg = nets[npart_ids[nid + BRC::BOFST]].pins.size();
            for(int i = t_dsu_ofst; i < t_dsu_ofst + deg; i++) {
                BRC::thr_dsu_space[i] = i - t_dsu_ofst;
            }
            #if BRVSK_CPU_ASSERT == 1
            SPDLOG_ASSERT(deg <= t_nlen, "{} {}", deg, t_nlen);
            #endif            
            while(true) {
                for(int i = t_nofst; i < t_nofst + deg; i++) {
                    BRC::thr_net_node_dist[i] = MAX_OFGPU_COST;
                }
                for(int i = t_nofst; i < t_nofst + deg; i++) {
                    BRC::thr_node_min[i] = -1;
                }
                int eid = 0;
                for(int i = 0; i < deg; i++) {                    
                    for(int j = i + 1; j < deg; j++) {
                        gcb_dt cedist = BRC::edge_dist[eofst + eid++];
                        int pnd1 = gdsu(t_dsu_ofst, i), pnd2 = gdsu(t_dsu_ofst, j);
                        if(pnd1 == pnd2) continue;
                        if(BRC::thr_net_node_dist[pnd1 + t_nofst] > cedist) {
                            BRC::thr_net_node_dist[t_nofst + pnd1] = cedist;
                            BRC::thr_node_min[t_nofst + pnd1] = pnd2;
                        }
                        if(BRC::thr_net_node_dist[pnd2 + t_nofst] > cedist) {
                            BRC::thr_net_node_dist[pnd2 + t_nofst] = cedist;
                            BRC::thr_node_min[t_nofst + pnd2] = pnd1;
                        }
                    }
                }
                #if BRVSK_CPU_ASSERT == 1
                SPDLOG_ASSERT(eid == elen, "{} {}", eid, elen);
                #endif
                //now do dsu
                int nskips = 0;
                for(int i = 0; i < deg; i++) {
                    if(gdsu(t_dsu_ofst, i) != i) {
                        nskips++;
                        continue;
                    }
                    int nxt_nd = BRC::thr_node_min[i + t_dsu_ofst];
                    nxt_nd = gdsu(t_dsu_ofst, nxt_nd);
                    int cnd = gdsu(t_dsu_ofst, i);
                    if(nxt_nd == cnd) continue;
                    #if BRVSK_CPU_ASSERT == 1
                    int cnid = npart_ids[nid + BRC::BOFST];
                    SPDLOG_ASSERT(cnd == i, "{} {}", i, cnd);
                    SPDLOG_ASSERT(nxt_nd != cnd, "{} {} {} {}", nxt_nd, cnd, deg, nets[cnid].net_name);
                    #endif
                    BRC::thr_dsu_space[t_dsu_ofst + cnd] = nxt_nd;
                    BRC::net_dist[t_nofst + cnd] += BRC::thr_net_node_dist[i + t_dsu_ofst];
                }
                if(nskips == deg - 1) {
                    break;
                }
            }
        }
    }
    void get_mst() {
        BRC::net_nodes_cnt.clear(); BRC::net_nodes_cnt.resize(BRC::BLEN);
        BRC::net_nodes_cnt_pfx.clear(); BRC::net_nodes_cnt_pfx.resize(BRC::BLEN);
        BRC::net_edges_cnt.clear(); BRC::net_edges_cnt.resize(BRC::BLEN);
        BRC::net_edges_cnt_pfx.clear(); BRC::net_edges_cnt_pfx.resize(BRC::BLEN);
        BRC::thr_dsu_pfx.clear(); BRC::thr_dsu_pfx.resize(BRC::NUM_THREADS, 0);
        BRC::thr_dsu_size.clear(); BRC::thr_dsu_size.resize(BRC::NUM_THREADS, 0);
        BRC::net_res_size.clear(); BRC::net_res_size.resize(BRC::BLEN, 0);
        BRC::net_res_pfx.clear(); BRC::net_res_pfx.resize(BRC::BLEN);
        BRC::net_dist.clear(); BRC::net_dist.resize(BRC::BLEN);
        for(int i = 0; i < BRC::NUM_THREADS; i++) {
            BRC::thrds[i] = thread(get_edge_data, i);
        }
        for(int i = 0; i < BRC::NUM_THREADS; i++) {
            BRC::thrds[i].join();
        }
        UTIL::get_vec_epfx<long>(BRC::net_nodes_cnt_pfx, BRC::net_nodes_cnt, BRC::BLEN);
        UTIL::get_vec_epfx<long>(BRC::net_edges_cnt_pfx, BRC::net_edges_cnt, BRC::BLEN);
        UTIL::get_vec_epfx<int>(BRC::thr_dsu_pfx, BRC::thr_dsu_size, BRC::NUM_THREADS);
        UTIL::get_vec_epfx<int>(BRC::net_res_pfx, BRC::net_res_size, BRC::BLEN);
        BRC::thr_dsu_space.clear(); BRC::thr_dsu_space.resize(BRC::thr_dsu_size[BRC::NUM_THREADS - 1] + BRC::thr_dsu_pfx[BRC::NUM_THREADS - 1]);
        BRC::thr_net_node_dist.clear(); BRC::thr_net_node_dist.resize(BRC::thr_dsu_size[BRC::NUM_THREADS - 1] + BRC::thr_dsu_pfx[BRC::NUM_THREADS - 1]);
        BRC::thr_node_min.clear(); BRC::thr_node_min.resize(BRC::thr_dsu_size[BRC::NUM_THREADS - 1] + BRC::thr_dsu_pfx[BRC::NUM_THREADS - 1]);
        int num_mst_edges = BRC::net_res_pfx[BRC::BLEN - 1] + BRC::net_res_size[BRC::BLEN - 1];
        BRC::mst_edges.clear(); BRC::mst_edges.resize(num_mst_edges);
        for(int i = 0; i < BRC::NUM_THREADS; i++) {
            BRC::thrds[i] = thread(t_get_mst, i);
        }
        for(int i = 0; i < BRC::NUM_THREADS; i++) {
            BRC::thrds[i].join();
        }
    }
}
