#pragma once

#include "./sq_common.h"
#include "./sq_util.hpp"
#include "./l_route.hpp"
#include "./djkst_route_bends.hpp"

namespace SQC = SQ_COMMON_SPACE;
namespace SQU = SQ_UTIL_SPACE;
namespace LRS = L_ROUTE_SPACE;
namespace DRBS = DJKST_RB_SPACE;

namespace MST_RSMT_SPACE {
    class mst_rsmt : public SQC::rsmt_cnst {
        public:
        LRS::L_ROUTE lpr;
        mst_rsmt() {
            lpr = LRS::L_ROUTE();
        }
        SQC::rsmt_t cnst_rsmt(int net_id) const override {
            // spdlog::info("in mst rsmt cnst");
            int nnid = net_ids[net_id];
            net_t cnet = nets[nnid];
            for(int i = 0; i < cnet.pins.size(); i++) {
                coord_t ccd = cnet.pins[i].pcd;
                gvia_dt zs = ccd.z;
                gcd_dt xs = ccd.x, ys = ccd.y;
                genc_dt senc = UTIL::get_enc<CUBE_ENC>(zs, xs, ys);
                for(int j = 0; j < i; j++) {
                    coord_t ncd = cnet.pins[j].pcd;
                    gvia_dt zd = ncd.z;
                    gcd_dt xd = ncd.x, yd = ncd.y;
                    genc_dt denc = UTIL::get_enc<CUBE_ENC>(zd, xd, yd);
                    SQC::pr_route_t crt = lpr.route(senc, denc);
                    // gcb_dt cres = DRBS::get_dist(senc, denc, 1);
                    spdlog::info("result is {}", crt.edge_dist);
                }
            }
            return SQC::rsmt_t{};
        }
    };
}
