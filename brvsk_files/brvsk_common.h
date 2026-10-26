#pragma once
#include "../common.h"
#include "../util.hpp"

#define BRVSK_CPU_ASSERT 0
#define BRVSK_TREE 0
namespace BRVSK_COMMON {
    int BOFST, BT_ID, BLEN;
    gcb_dt *edge_dist;
    long NUM_EDGES;
    long NUM_THREADS = 10;
    vector<long> net_nodes_cnt, net_nodes_cnt_pfx, net_edges_cnt, net_edges_cnt_pfx;
    vector<gcb_dt> net_dist, thr_net_node_dist;
    vector<thread> thrds;
    vector<int> thr_dsu_size, thr_dsu_pfx, thr_dsu_space, net_res_size, net_res_pfx, thr_node_min;
    vector<long> mst_edges;
}