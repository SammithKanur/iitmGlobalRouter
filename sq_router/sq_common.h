#pragma once
#include "../common.h"
#include "../util.hpp"

using namespace COMMON_DATA;
#define SQ_XDIR 1
#define SQ_YDIR 2
#define SQ_NDIR 0
#define SQR_CPU_ASSERT 1

namespace SQ_COMMON_SPACE {
    vector<gcb_dt> pth_mat;
    struct interval_t {
        genc_dt enc1, enc2;
    };
    struct pr_route_t {
        gcb_dt edge_dist;
        vector<interval_t> route;
        vector<genc_dt> path_pts;
    };
    struct rsmt_t {
        vector<interval_t> route;
    };
    class PRR {
        public:
        virtual pr_route_t route(genc_dt src, genc_dt dst) const {
            return pr_route_t{};
        }
    };
    class rsmt_cnst {
        public:
        virtual rsmt_t cnst_rsmt(int net_id) const {
            return rsmt_t{};
        }
    };
}
