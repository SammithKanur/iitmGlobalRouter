#pragma once
#include "./writer_common.h"
#include "./writer_util.hpp"

using namespace std;
using namespace COMMON_DATA;
namespace RWC = ROUTE_WRITER_COMMON;
namespace MSC = MST_COMMON;
namespace RWU = ROUTE_WRITER_UTIL;

namespace ROUTE_WRITER {
    void write_route() {
        RWC::BT_ID = 0;
        RWC::thrds.clear(); RWC::thrds.resize(RWC::NUM_THREADS);
        write_route_bsem.acquire();
        while(RWC::BT_ID < PART_ID) {
            RWC::BOFST = npart_ofst[RWC::BT_ID];
            RWC::BLEN = npart_len[RWC::BT_ID];
            RWC::net_block_len1 = MSC::net_block_len1.data() + RWC::BOFST;
            RWC::net_block_len = MSC::net_block_len.data() + RWC::BOFST;
            RWC::net_block_pfx1 = MSC::net_block_pfx1.data() + RWC::BOFST;
            RWC::net_block_pfx = MSC::net_block_pfx.data() + RWC::BOFST;
            RWC::net_nodes_cnt = MSC::net_nodes_cnt.data() + RWC::BOFST;
            RWC::net_nodes_cnt_pfx = MSC::net_nodes_cnt_pfx.data() + RWC::BOFST;
            RWC::net_mst_edges = MSC::net_mst_edges.data() + RWC::NUM_NODES_OFST;
            RWC::mst_lpat_pat = MSC::mst_lpat_pat.data() + RWC::NUM_NODES_OFST;
            // RWU::write_route(0);
            for(int i = 0; i < RWC::NUM_THREADS; i++) {
                if(WRT_MRG) {
                    RWC::thrds[i] = thread(RWU::mrg_route_write, i);    
                } else {
                    RWC::thrds[i] = thread(RWU::write_route, i);
                }
            }
            for(int i = 0; i < RWC::NUM_THREADS; i++) {
                RWC::thrds[i].join();
            }
            int max_nodes_len = RWC::net_nodes_cnt_pfx[RWC::BLEN - 1] + RWC::net_nodes_cnt[RWC::BLEN - 1];
            if(RWC::BT_ID == 0) {
                for(int i = 0; i < RWC::NUM_THREADS; i++) {
                    RWC::thrds[i] = thread(RWU::write_odeg_route, i);
                }
                for(int i = 0; i < RWC::NUM_THREADS; i++) {
                    RWC::thrds[i].join();
                }   
            }
            RWC::BT_ID++;
            RWC::NUM_NODES_OFST += max_nodes_len;
            write_route_bsem.acquire();
        }
    }
}
