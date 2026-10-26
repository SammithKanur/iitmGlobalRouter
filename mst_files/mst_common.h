#pragma once
#include "../common.h"
#include "../util.hpp"

using namespace std;
using namespace COMMON_DATA;

typedef long ednd_dt;
typedef long edid_dt;

#define MST_GPU_ZIP_SORT 1
#define MST_BRVSK 1
#define MST_NODE_BITS 30
#define MST_NODE_MASK ((1L << MST_NODE_BITS) - 1L)
#define MST_ASSERT 0
#define MST_GPU_ASSERT 0
#define MST_CPU_GPU_COMPARISON 0
#define TEST_MST_GPU_MAP 0
#define MST_TREE 1
#define MST_CONN_CHECK 0
#define MST_NETID_BITS 28
#define MST_NETID_MASK ((1 << MST_NETID_BITS) - 1)
#define MST_EDG_ORDER_CHECK 0
#define MST_TREE_GPU 1
#define NODE_BITS 15
#define NNODE_BITS (NODE_BITS + NODE_BITS)
#define NET_BITS 28
#define NODE_MASK ((1L << NODE_BITS) - 1L)
#define TEST_MST_GPU_SORT 0
#define TEST_MST_GPU_DSU_INIT 0
#define TEST_BRVSK_MST_FORM 0
#define BRVSK_VIS_BITS 22
#define BRVSK_VIS_MASK ((1 << BRVSK_VIS_BITS) - 1)
#define ADD_GPU_MST_EDGES 1
#define MST_GPU_CONN_CHECK 0
#define EDGES_CNT_GPU 1
#define MST_VIA_FULL_WT 1.0
#define MST_VIA_FRAC_WT 0.5
#define MST_WIRE_WT 1.0
#define MST_FEEDBACK 1

#define ADD_MST_PERF_DATA 0
#define ADD_MST_PERF_ALL_BATCH_DATA 0
#define BRVSK_VERSION 0
#define KRSK_CPU_VERSION 3
#define EDGE_NODE_MAP_STORE 0
#define TIME_STAMP_BITS 25
#define TIME_STAMP_MASK ((1 << TIME_STAMP_BITS) - 1)
#define EDGE_VAL ((1 << TIME_STAMP_BITS) | 1)
#define TIMESTAMP_VAL ((1 << TIME_STAMP_BITS))
#define MST_COMMIT_EDGES 1
#define MST_MIG 1

typedef long nenc_dt;
typedef int gccnt_dt;

namespace MST_COMMON {
    struct edge_t {
        edid_dt ed_id;
        ednd_dt ed_nds;
        gcb_dt dist;
    };
    edge_t *edge_arr;
    gcb_dt *edge_dist;
    int BT_ID, BOFST, BLEN, NUM_NODES_OFST = 0;
    long NUM_EDGES;
    int MAX_DEGREE;
    long MST_THREADS = 20;
    vector<long> net_edges_cnt, net_edges_cnt_pfx;
    vector<thread> thrds;
    vector<int> thr_dsu_size, thr_dsu_pfx, thr_dsu_space, net_res_size, net_res_pfx;
    vector<int> net_nodes_cnt, net_nodes_cnt_pfx;
    vector<long> mst_edges;
    vector<gcb_dt> net_dist;
    vector<edge_t> mst_conn_check_edges;
    __managed__ long *d_net_block_pfx, *d_net_block_len, *d_net_block_map, *d_net_block_map1, *d_net_block_pfx1, *d_net_block_len1;
    __managed__ long *d_bit_edge_cell_cnt;
    __managed__ long *d_net_edge_id, *d_brvsk_edges, *td_net_block_pfx, *td_net_block_len, *td_net_block_pfx1, *td_net_block_len1;
    __managed__ int *d_net_res_size, *d_net_res_pfx, *d_net_nodes_cnt, *d_net_nodes_cnt_pfx, *d_dsu_data, *td_net_nodes_cnt_pfx, *td_net_nodes_cnt;
    __managed__ int *d_net_node_id, *d_brvsk_node, *d_brvsk_edge_cnt, *d_brvsk_net_len;
    __managed__ gccnt_dt *d_net_edges_gcell_cnt, *d_net_edges_gcell_cnt_pfx;
    __managed__ bool *d_failed_net, *d_is_nets_failed;
    __managed__ gcb_dt *d_net_dist, *d_brvsk_node_val;
    __managed__ nenc_dt *d_net_mst_edges, *td_net_mst_edges;
    __managed__ gvia_dt *d_mst_lpat_pat, *td_mst_lpat_pat;
    vector<long> net_block_pfx, net_block_len, net_block_pfx1, net_block_len1;
    vector<nenc_dt> net_mst_edges;
    vector<gvia_dt> mst_lpat_pat;
}
