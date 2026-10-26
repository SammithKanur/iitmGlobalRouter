#pragma once
#include "brvsk_common.h"
#include "../l_pattern_files/lpattern_common.h"
#include "../l_pattern_files/lpattern_util.hpp"
#include "brvsk_util.hpp"

using namespace COMMON_DATA;
using namespace std;

namespace BRC = BRVSK_COMMON;
namespace BRU = BRVSK_UTIL;
namespace LPC = LPAT_COMMON;
namespace LPU = LPAT_UTIL;

namespace BRVSK_BUILD {
    void get_mst_tree(int &bt_id, long &num_edges) {
        BRC::BT_ID = bt_id; BRC::NUM_EDGES = num_edges;
        BRC::BOFST = npart_ofst[bt_id]; BRC::BLEN = npart_len[bt_id];
        BRC::edge_dist = (gcb_dt *)malloc(sizeof(gcb_dt) * BRC::NUM_EDGES);
        CUDAMEMCPYDH_1D(BRC::edge_dist, LPC::d_lpat_res, BRC::NUM_EDGES, gcb_dt);
        BRC::thrds.clear(); BRC::thrds.resize(BRC::NUM_THREADS);
        BRU::get_mst();        
    }
}
