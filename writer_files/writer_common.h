#pragma once

#include "../common.h"
#include "../mst_files/mst_common.h"
#include "../l_pattern_files/lpattern_util.hpp"

#define ROUTE_WRITE_CPU_ASSERT 1

namespace LPU = LPAT_UTIL;

namespace ROUTE_WRITER_COMMON {
    long *net_block_pfx1, *net_block_pfx, *net_block_len, *net_block_len1;
    int *net_nodes_cnt, *net_nodes_cnt_pfx;
    nenc_dt *net_mst_edges;
    gvia_dt *mst_lpat_pat;
    int NUM_THREADS = 4, BT_ID = 0;
    vector<thread> thrds;
    int BOFST, BLEN, NUM_NODES_OFST = 0;
}
