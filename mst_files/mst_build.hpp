#pragma once
#include "mst_common.h"
#include "mst_util.hpp"
#include "../l_pattern_files/lpattern_common.h"
#include "../l_pattern_files/lpattern_util.hpp"

using namespace std;
using namespace COMMON_DATA;

namespace MSC = MST_COMMON;
namespace MSU = MST_UTIL;
namespace LPC = LPAT_COMMON;
namespace LPU = LPAT_UTIL;

namespace MST_BUILD {
    template<int pat_type>
    void get_mst_tree(int &bt_id, long &num_edges) {
        #if ADD_MST_PERF_DATA == 1
        auto start_time = chrono::high_resolution_clock::now();
        #endif
        MSC::BT_ID = bt_id; MSC::NUM_EDGES = num_edges;
        MSC::BOFST = npart_ofst[bt_id]; MSC::BLEN = npart_len[bt_id];
        MSC::edge_dist = (gcb_dt *)malloc(sizeof(gcb_dt) * MSC::NUM_EDGES);
        if(pat_type == 0) {
            CUDAMEMCPYDH_1D(MSC::edge_dist, LPC::d_lpat_res, MSC::NUM_EDGES, gcb_dt);
        } else {
            // CUDAMEMCPYDH_1D(MSC::edge_dist, ZPC::d_zpat_res, MSC::NUM_EDGES, gcb_dt);
        }
        MSC::thrds.clear(); MSC::thrds.resize(MSC::MST_THREADS);
        MSU::get_mst();
        #if MST_EDG_ORDER_CHECK == 1
        // SPDLOG_ASSERT(num_edges == MSC::mst_edges.size(), "num_edges {} mst_edges.size() {}", num_edges, MSC::mst_edges.size());
        int tncnt = 1;
        long ted1 = MSC::mst_edges[0], ted2;
        for(int i = 1; i < MSC::mst_edges.size(); i++) {
            ted2 = MSC::mst_edges[i];
            int n1 = ted1 & MST_NETID_MASK, n2 = ted2 & MST_NETID_MASK;
            if(n1 != n2) {
                SPDLOG_ASSERT(tncnt == nets[n1].pins.size() - 1, "ncnt {} actual_degree {}", tncnt, nets[n1].pins.size() - 1);
                tncnt = 0;
            }
            tncnt++;
            ted1 = ted2;
        }
        int tn1 = ted1 & MST_NETID_MASK;
        SPDLOG_ASSERT(tncnt == nets[tn1].pins.size() - 1, "ncnt {} actual_degree {}", tncnt, nets[tn1].pins.size() - 1);
        #endif
        #if COMMIT_LPAT == 1
        long ed1 = MSC::mst_edges[0], ed2,eeid;
        int n1 = ed1 & MST_NETID_MASK, n2;
        UTIL::write_net_name(nets[n1].net_name);
        for(int i = 1; i < MSC::mst_edges.size(); i++) {
            ed2 = MSC::mst_edges[i];
            n1 = ed1 & MST_NETID_MASK; n2 = ed2 & MST_NETID_MASK;
            if(pat_type == 0) {
                eeid = ed1 >> MST_NETID_BITS;
                LPU::write_pattern(eeid, n1);
            } else {
                eeid = ed1 >> MST_NETID_BITS;
                // ZPU::write_pattern(eeid, n1);
            }
            if(n1 != n2) {
                UTIL::write_closed_bracket();
                UTIL::write_net_name(nets[n2].net_name);
            }
            ed1 = ed2;
        }
        n1 = ed1 & MST_NETID_MASK;
        if(pat_type == 0) {
            eeid = ed1 >> MST_NETID_BITS;
            LPU::write_pattern(eeid, n1);
        } else {
            eeid = ed1 >> MST_NETID_BITS;
            // ZPU::write_pattern(eeid, n1);
        }
        UTIL::write_closed_bracket();
        #endif
        free(MSC::edge_dist);
        #if ADD_MST_PERF_DATA == 1 && false
        auto end_time = chrono::high_resolution_clock::now();
        long total_time_taken_ms = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
        #if MST_TREE_GPU == 1 && MST_BRVSK == 1
        DBI::insert_mst_perf_data(BRVSK_VERSION, bt_id, MSC::BLEN, total_time_taken_ms);
        #elif MST_TREE_GPU == 0
        DBI::insert_mst_perf_data(KRSK_CPU_VERSION, bt_id, MSC::BLEN, total_time_taken_ms);
        #endif
        #endif
    }
}
