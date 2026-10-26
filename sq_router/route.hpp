#pragma once
#include "./sq_common.h"
#include "./sq_util.hpp"
#include "./mst_rsmt.hpp"

namespace SQC = SQ_COMMON_SPACE;
namespace SQU = SQ_UTIL_SPACE;
namespace MRS = MST_RSMT_SPACE;

namespace SQ_ROUTER {
    void sq_route() {
        net_ids.clear();
        net_ids.resize(NUM_NETS);
        SQC::pth_mat.clear();
        SQC::pth_mat.resize(X * Y * Z);
        for(int i = 0; i < NUM_NETS; i++) {
            net_ids[i] = i;
        }
        sort(net_ids.begin(), net_ids.begin() + NUM_NETS, [](int a, int b) {
            int apsize = nets[a].pins.size(), bpsize = nets[b].pins.size();
            if(apsize == bpsize) {
                return nets[a].net_name < nets[b].net_name;
            }
            return apsize < bpsize;
        });
        SQU::build_cost_mat();
        SQU::get_cost_mat_pfx();
        MRS::mst_rsmt rtr = MRS::mst_rsmt();
        spdlog::info("calling cnst rsmt");
        for(int i = 0; i < NUM_NETS; i++) {
            rtr.cnst_rsmt(i);
        }
    }
}
