#pragma once
#include "writer_common.h"

using namespace std;
using namespace COMMON_DATA;
namespace RWC = ROUTE_WRITER_COMMON;

namespace ROUTE_WRITER_UTIL {
    __inline__ void write_net_name(string &res_str, string &nname) {
        nname.erase(std::remove(nname.begin(), nname.end(), '\n'), nname.end());
        nname.erase(std::remove(nname.begin(), nname.end(), ' '), nname.end());
        res_str += nname + "\n(\n";
    }
    __inline__ void write_closed_bracket(string &res_str) {
        res_str += ")";
    }
    __inline__ long get_gran_cd(long cd) {
        return 4200L * cd + 2100L;
    }
    void add_vias_cpu(string &res_str, long x, long y, long z1, long z2) {
      if(z1 > z2) {swap(z1, z2);}
      for(int i = z1;i < z2; i++) {
        if(!RSEG_FMT) {
            res_str += to_string(i) + " " + to_string(x) + " " + to_string(y) + " " + to_string(i) + " " + to_string(x) + " " + to_string(y) + "\n";
        } else {
            res_str += to_string(get_gran_cd(x)) + " " + to_string(get_gran_cd(y)) + " metal" + to_string(i + 1) + " " + to_string(get_gran_cd(x)) + " " + to_string(get_gran_cd(y)) + " metal" + to_string(i + 2) + "\n";
        }
      }
    }
    void add_horizontal_cpu(string &res_str, long z, long y, long x1, long x2) {
        if(x1 > x2) swap(x1, x2);
        // SPDLOG_ASSERT(x1 < x2, "{} {}", x1, x2);
        if(x1 == x2) {
            return;
        }
        if(!RSEG_FMT) {
            res_str += to_string(z) + " " + to_string(x1) + " " + to_string(y) + " " + to_string(z) + " " + to_string(x2) + " " + to_string(y) + "\n";
        } else {
            res_str += to_string(get_gran_cd(x1)) + " " + to_string(get_gran_cd(y)) + " metal" + to_string(z + 1) + " " + to_string(get_gran_cd(x2)) + " " + to_string(get_gran_cd(y)) + " metal" + to_string(z + 2) + "\n";
        }
    }
    void add_vertical_cpu(string &res_str, long z, long x, long y1, long y2) {
        if(y1 > y2) swap(y1, y2);
        // SPDLOG_ASSERT(y1 < y2, "{} {}", y1, y2);
        if(y1 == y2) {
            return;
        }
        if(!RSEG_FMT) {
            res_str += to_string(z) + " " + to_string(x) + " " + to_string(y1) + " " + to_string(z) + " " + to_string(x) + " " + to_string(y2) + "\n";
        } else {
            res_str += to_string(get_gran_cd(x)) + " " + to_string(get_gran_cd(y1)) + " metal" + to_string(z + 1) + " " + to_string(get_gran_cd(x)) + " " + to_string(get_gran_cd(y2)) + " metal" + to_string(z + 2) + "\n";
        }
    }
    void add_vias_cpu_vec(vector<minterval_t> &edgs, long x, long y, long z1, long z2) {
      if(z1 > z2) {swap(z1, z2);}
      for(int i = z1;i < z2; i++) {
        genc_dt enc1 = UTIL::get_enc<WIRE_ENC>(i, x, y);
        edgs.push_back(minterval_t(enc1, enc1));
      }
    }
    void add_horizontal_cpu_vec(vector<minterval_t> &edgs, long z, long y, long x1, long x2) {
        if(x1 > x2) swap(x1, x2);
        // SPDLOG_ASSERT(x1 < x2, "{} {}", x1, x2);
        if(x1 == x2) {
            return;
        }
        genc_dt enc1 = UTIL::get_enc<WIRE_ENC>(z, x1, y), enc2 = UTIL::get_enc<WIRE_ENC>(z, x2, y);
        edgs.push_back(minterval_t(enc1, enc2));
    }
    void add_vertical_cpu_vec(vector<minterval_t> &edgs, long z, long x, long y1, long y2) {
        if(y1 > y2) swap(y1, y2);
        // SPDLOG_ASSERT(y1 < y2, "{} {}", y1, y2);
        if(y1 == y2) {
            return;
        }
        genc_dt enc1 = UTIL::get_enc<WIRE_ENC>(z, x, y1), enc2 = UTIL::get_enc<WIRE_ENC>(z, x, y2);
        edgs.push_back(minterval_t(enc1, enc2));
    }
    void write_odeg_route(int tid) {
        for(int i = tid; i < one_deg_nets.size(); i += RWC::NUM_THREADS) {
            #if ROUTE_WRITE_CPU_ASSERT == 1
            SPDLOG_ASSERT(one_deg_nets[i].pins.size() == 1, "size {}", one_deg_nets.size());
            #endif
            pin_t cpin = one_deg_nets[i].pins[0];
            long z, x, y;
            z = cpin.pcd.z; x = cpin.pcd.x; y = cpin.pcd.y;
            string res_str = "";
            write_net_name(res_str, one_deg_nets[i].net_name);
            add_vias_cpu(res_str, x, y, z, z + 1);
            write_closed_bracket(res_str);
            gr_file_logger->info(res_str);
        }
    }
    void write_route(int tid) {
        for(int i = tid; i < RWC::BLEN; i += RWC::NUM_THREADS) {
            int deg = RWC::net_nodes_cnt[i];
            int deg_pfx = RWC::net_nodes_cnt_pfx[i];
            int cnid = npart_ids[i + RWC::BOFST];
            net_t cnet = nets[cnid];
            string net_gr = "";
            write_net_name(net_gr, cnet.net_name);
            // spdlog::info("[deg, deg_pfx] [{},{}]", deg, deg_pfx);
            for(int j = deg_pfx; j < deg_pfx + deg; j++) {
                nenc_dt mst_cur_edge = RWC::net_mst_edges[j];
                if(mst_cur_edge == -1)
                    continue;
                int n1 = (mst_cur_edge & NODE_MASK), n2 = int((mst_cur_edge >> NODE_BITS) & NODE_MASK), cnet_id = int(mst_cur_edge >> NNODE_BITS);
                if(n1 < n2) swap(n1, n2);
                #if ROUTE_WRITE_CPU_ASSERT == 1
                SPDLOG_ASSERT(cnet_id == i, "{} {} {}", cnet_id, i, RWC::BT_ID);
                #endif
                gvia_dt clpat = RWC::mst_lpat_pat[j], via1, via2;
                LPU::get_zdec(clpat, via1, via2);
                #if ROUTE_WRITE_CPU_ASSERT == 1
                SPDLOG_ASSERT(n1 < cnet.pins.size() && n2 < cnet.pins.size(), "[n1,n2,pins_size,net_name] [{},{},{},{}]", n1, n2, cnet.pins.size(), cnet.net_name);
                #endif
                gvia_dt z1, z2;
                gcd_dt x1, x2, y1, y2;
                z1 = cnet.pins[n1].pcd.z; z2 = cnet.pins[n2].pcd.z; 
                x1 = cnet.pins[n1].pcd.x; x2 = cnet.pins[n2].pcd.x;
                y1 = cnet.pins[n1].pcd.y; y2 = cnet.pins[n2].pcd.y;
                bool bth_diff = false;
                if(x1 != x2 && y1 != y2) {
                    bth_diff = true;
                }
                if(z1 != via1) {
                    add_vias_cpu(net_gr, (x1), (y1), z1, via1);
                    z1 = via1;
                }
                if(z1 > 0) {
                    if(z1 & 1) {
                        add_vertical_cpu(net_gr, z1, (x1), (y1), (y2));
                        y1 = y2;
                    } else {
                        add_horizontal_cpu(net_gr, z1, (y1), (x1), (x2));
                        x1 = x2;
                    }
                }
                via2 = ((z1 & 1) ? 2 * via2 : via2 * 2 + 1);
                if(bth_diff && via2 != 0) {
                    if(z1 != via2) {
                        add_vias_cpu(net_gr, (x1), (y1), z1, via2);
                        z1 = via2;
                        if(z1 & 1) {
                            add_vertical_cpu(net_gr, z1, (x1), (y1), (y2));
                            y1 = y2;
                        } else {
                            add_horizontal_cpu(net_gr, z1, (y1), (x1), (x2));
                            x1 = x2;
                        }
                    }
                }
                if(z1 != z2) {
                    add_vias_cpu(net_gr, (x1), (y1), z1, z2);
                    z1 = z2;
                }
                #if ROUTE_WRITE_CPU_ASSERT == 1
                SPDLOG_ASSERT(z1 == z2 && x1 == x2 && y1 == y2, "[zxy] [{},{},{}] [{},{},{}] {}", z1, x1, y1, z2, x2, y2, clpat);
                #endif
                //start adding routes now
            }
            write_closed_bracket(net_gr);
            gr_file_logger->info(net_gr);
        }
    }
    void merge_st_edges_nnew(vector<minterval_t> &st_edges) {
        if(st_edges.size() == 0) {
            return;
        }
        vector<minterval_t> nst_edges;
        nst_edges.clear();
        map<genc_dt, bool> cmp; cmp.clear();
        for(minterval_t &cintv : st_edges) {
            genc_dt enc1 = cintv.enc1, enc2 = cintv.enc2;
            if(enc1 == enc2) {
                cmp[enc1] = true;
            } else {
                nst_edges.push_back(cintv);
            }
        }
        if(nst_edges.size()) {
            sort(nst_edges.begin(), nst_edges.end(), [](minterval_t &a, minterval_t &b) {
                if(a.enc1 == b.enc1) {
                    return a.enc2 < b.enc2;
                }
                return a.enc1 < b.enc1;
            });
            minterval_t pintv = nst_edges[0];
            int fid = 0;
            for(int i = 1; i < nst_edges.size(); i++) {
                minterval_t cintv = nst_edges[i];
                if(pintv.enc2 >= cintv.enc1) {
                    pintv.enc2 = max(cintv.enc2, pintv.enc2);
                } else {
                    nst_edges[fid++] = pintv;
                    pintv = cintv;
                }
            }
            nst_edges[fid++] = pintv;
            while(nst_edges.size() > fid) {
                nst_edges.pop_back();
            }
        }
        for(auto mp_nd : cmp) {
            genc_dt cenc = mp_nd.first;
            nst_edges.push_back(minterval_t(cenc, cenc));
        }
        st_edges.clear();
        st_edges = move(nst_edges);
    }
    void add_segment(string &res_str, int z1, int x1, int y1, int z2, int x2, int y2) {
        if(RSEG_FMT) {
            // res_str += to_string(z1) + " " + to_string(x1) + " " + to_string(y1) + " " + to_string(z2) + " " + to_string(x2) + " " + to_string(y2) + "\n";
            if(x1 == x2 && y1 == y2 && z1 == z2) {
                SPDLOG_ASSERT(z1 < 9, "z1 = {}", z1);
                res_str += to_string(get_gran_cd(x1)) + " " + to_string(get_gran_cd(y1)) + " metal" + to_string(z1 + 1) + " " + to_string(get_gran_cd(x2)) + " " + to_string(get_gran_cd(y2)) + " metal" + to_string(z1 + 2) + "\n";
            } else {
                res_str += to_string(get_gran_cd(x1)) + " " + to_string(get_gran_cd(y1)) + " metal" + to_string(z1 + 1) + " " + to_string(get_gran_cd(x2)) + " " + to_string(get_gran_cd(y2)) + " metal" + to_string(z1 + 1) + "\n";
            }
        } else {
            res_str += to_string(z1) + " " + to_string(x1) + " " + to_string(y1) + " " + to_string(z2) + " " + to_string(x2) + " " + to_string(y2) + "\n";
        }
    }
    void write_segment(string &res_str, minterval_t &cintv) {
        genc_dt enc1 = cintv.enc1, enc2 = cintv.enc2;
        gvia_dt z1, z2;
        gcd_dt x1, x2, y1, y2;
        UTIL::get_dec_all<CUBE_ENC>(enc1, z1, x1, y1);
        UTIL::get_dec_all<CUBE_ENC>(enc2, z2, x2, y2);
        add_segment(res_str, z1, x1, y1, z2, x2, y2);
    }
    void add_net_route(string &res_str, vector<minterval_t> &edgs) {
        for(minterval_t &cintv : edgs) {
            write_segment(res_str, cintv);
        }
        write_closed_bracket(res_str);
        gr_file_logger->info(res_str);
    }
    void oprd_ppatch(net_t &cnet, vector<minterval_t> &edgs) {
        map<genc_dt, bool> cmp;
        cmp.clear();
        auto cmp1 = [](minterval_t a, minterval_t b) {
            if(a.enc1 == b.enc1) {
                return a.enc2 < b.enc2;
            }
            return a.enc1 < b.enc1;
        };
        set<minterval_t, decltype(cmp1)> cst(cmp1); cst.clear();
        for(minterval_t cintv : edgs) {
            if(cintv.enc1 == cintv.enc2) {
                gvia_dt cz;
                gcd_dt cx, cy;
                UTIL::get_dec_all<CUBE_ENC>(cintv.enc1, cz, cx, cy);
                genc_dt enc1 = UTIL::get_enc<WIRE_ENC>(cz, cx, cy), enc2 = UTIL::get_enc<WIRE_ENC>(cz + 1, cx, cy);
                cmp[enc1] = true; cmp[enc2] = true;
            } else {
                cmp[cintv.enc1] = true; cmp[cintv.enc2] = true;
                cst.insert(cintv);
            }
        }
        for(pin_t &cpin : cnet.pins) {
            int cx = cpin.pcd.x, cy = cpin.pcd.y, cz = cpin.pcd.z;
            genc_dt cenc = UTIL::get_enc<WIRE_ENC>(cz, cx, cy);
            if(cz == 0) {
                continue;
            }
            if(cmp.count(cenc) == 0) {
                //endpoint does not exist on the interval
                spdlog::info("repairing endpoint {}", cenc);
                auto it = cst.upper_bound({minterval_t(cenc, cenc)});
                SPDLOG_ASSERT(it != cst.begin(), "this is not acceptable");
                it--;
                minterval_t cintv = *it;
                SPDLOG_ASSERT(cintv.enc1 < cenc && cenc < cintv.enc2, "[{} {}] {}", cintv.enc1, cintv.enc2, cenc);
                cst.erase(it);
                cst.insert(minterval_t(cintv.enc1, cenc)); cst.insert(minterval_t(cenc, cintv.enc2));
            }
        }
        vector<minterval_t> nedgs; nedgs.clear();
        for(minterval_t cintv : edgs) {
            if(cintv.enc1 == cintv.enc2) {
                nedgs.push_back(cintv);
            }
        }
        for(minterval_t cintv : cst) {
            nedgs.push_back(cintv);
        }
        edgs = move(nedgs);
    }
    void mrg_route_write(int tid) {
        for(int i = tid; i < RWC::BLEN; i += RWC::NUM_THREADS) {
            int deg = RWC::net_nodes_cnt[i];
            int deg_pfx = RWC::net_nodes_cnt_pfx[i];
            int cnid = npart_ids[i + RWC::BOFST];
            net_t cnet = nets[cnid];
            vector<minterval_t> edgs;
            edgs.clear();
            string net_gr = "";
            // spdlog::info("[deg, deg_pfx] [{},{}]", deg, deg_pfx);
            write_net_name(net_gr, cnet.net_name);
            for(int j = deg_pfx; j < deg_pfx + deg; j++) {
                nenc_dt mst_cur_edge = RWC::net_mst_edges[j];
                if(mst_cur_edge == -1)
                    continue;
                int n1 = (mst_cur_edge & NODE_MASK), n2 = int((mst_cur_edge >> NODE_BITS) & NODE_MASK), cnet_id = int(mst_cur_edge >> NNODE_BITS);
                if(n1 < n2) swap(n1, n2);
                #if ROUTE_WRITE_CPU_ASSERT == 1
                SPDLOG_ASSERT(cnet_id == i, "{} {} {}", cnet_id, i, RWC::BT_ID);
                #endif
                gvia_dt clpat = RWC::mst_lpat_pat[j], via1, via2;
                LPU::get_zdec(clpat, via1, via2);
                #if ROUTE_WRITE_CPU_ASSERT == 1
                SPDLOG_ASSERT(n1 < cnet.pins.size() && n2 < cnet.pins.size(), "[n1,n2,pins_size,net_name] [{},{},{},{}]", n1, n2, cnet.pins.size(), cnet.net_name);
                #endif
                gvia_dt z1, z2;
                gcd_dt x1, x2, y1, y2;
                z1 = cnet.pins[n1].pcd.z; z2 = cnet.pins[n2].pcd.z; 
                x1 = cnet.pins[n1].pcd.x; x2 = cnet.pins[n2].pcd.x;
                y1 = cnet.pins[n1].pcd.y; y2 = cnet.pins[n2].pcd.y;
                bool bth_diff = false;
                if(x1 != x2 && y1 != y2) {
                    bth_diff = true;
                }
                if(z1 != via1) {
                    add_vias_cpu_vec(edgs, (x1), (y1), z1, via1);
                    z1 = via1;
                }
                if(z1 > 0) {
                    if(z1 & 1) {
                        add_vertical_cpu_vec(edgs, z1, (x1), (y1), (y2));
                        y1 = y2;
                    } else {
                        add_horizontal_cpu_vec(edgs, z1, (y1), (x1), (x2));
                        x1 = x2;
                    }
                }
                via2 = ((z1 & 1) ? 2 * via2 : via2 * 2 + 1);
                if(bth_diff && via2 != 0) {
                    if(z1 != via2) {
                        add_vias_cpu_vec(edgs, (x1), (y1), z1, via2);
                        z1 = via2;
                        if(z1 & 1) {
                            add_vertical_cpu_vec(edgs, z1, (x1), (y1), (y2));
                            y1 = y2;
                        } else {
                            add_horizontal_cpu_vec(edgs, z1, (y1), (x1), (x2));
                            x1 = x2;
                        }
                    }
                }
                if(z1 != z2) {
                    add_vias_cpu_vec(edgs, (x1), (y1), z1, z2);
                    z1 = z2;
                }
                #if ROUTE_WRITE_CPU_ASSERT == 1
                SPDLOG_ASSERT(z1 == z2 && x1 == x2 && y1 == y2, "[zxy] [{},{},{}] [{},{},{}] {}", z1, x1, y1, z2, x2, y2, clpat);
                #endif
            }
            merge_st_edges_nnew(edgs);
            if(OPRD_PPATCH)
                oprd_ppatch(cnet, edgs);
            add_net_route(net_gr, edgs);
        }
    }
}
