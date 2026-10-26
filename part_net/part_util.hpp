#pragma once
#include "../common.h"

using namespace std;
using namespace COMMON_DATA;

namespace PART_UTIL {
    void test_part(bbx_t *bx_arr, int num_parts) {
        TIME_();
        for(int i = 0; i < num_parts; i++) {
            int ofst = npart_ofst[i], len = npart_len[i];
            bxtree_t crtree = bxtree_t(); crtree.clear();
            vector<box_t> res;
            for(int j = ofst; j < ofst + len; j++) {
                int net_id = npart_ids[j];
                bbx_t cbx = bx_arr[net_id];
                // box_t bx = box_t(point_t(0, 1), point_t(2, 3));
                box_t bx = box_t(point_t(cbx.minx, cbx.miny), point_t(cbx.maxx, cbx.maxy));
                res.clear();
                crtree.query(bgi::intersects(bx), back_inserter(res));
                if(res.size() > 0) {
                    spdlog::critical("part_id {} net_id {} current rectangle is intersecting {} {} {} {}", i, net_id, cbx.minx, cbx.miny, cbx.maxx, cbx.maxy);
                    // for(box_t crt : res) {
                    //     int inet_id = crt.second;
                    //     bbx_t ibx = bx_arr[inet_id];
                    //     spdlog::info("inet_id {} intersecting rectangle {} {} {} {}", inet_id, ibx.minx, ibx.miny, ibx.maxx, ibx.maxy);
                    // }
                    exit(0);
                }
                // spdlog::info("inserting {} {} {} {}", cbx.minx, cbx.miny, cbx.maxx, cbx.maxy);
                crtree.insert(bx);
            }
        }
        spdlog::info("no rectangles are overlapping in any partition");
    }
}
