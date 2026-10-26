#pragma once
#include "../common.h"
#include "../util.hpp"
#include "./part_util.hpp"
#include "../mst_files/mst_common.h"

using namespace std;
using namespace COMMON_DATA;
namespace CDT = COMMON_DATA;
namespace MSC = MST_COMMON;

namespace NET_PART {
    bool *vis_data;
    int TEST_DNM = 1;
    int REM_NETS, TKN, GP_ID;
    int *net_ids;
    bbx_t *bx_arr;
    gcd_dt CMAXX;
    void init_data() {
        vis_data = (bool *)malloc(sizeof(bool) * NUM_NETS); memset(vis_data, 0, sizeof(bool) * NUM_NETS);
        npart_len = (int *)malloc(sizeof(int) * NUM_NETS); memset(npart_len, 0, sizeof(int) * NUM_NETS);
        npart_ofst = (int *)malloc(sizeof(int) * NUM_NETS); memset(npart_ofst, 0, sizeof(int) * NUM_NETS);
        npart_ids = (int *)malloc(sizeof(int) * NUM_NETS);
        net_ids = (int *)malloc(sizeof(int) * NUM_NETS);
        CUDA_MALLOC_1D(MSC::td_net_nodes_cnt, NUM_NETS, int);
        CUDA_MALLOC_1D(MSC::td_net_nodes_cnt_pfx, NUM_NETS, int);
        CUDA_MALLOC_1D(MSC::td_net_block_pfx, NUM_NETS, long);
        CUDA_MALLOC_1D(MSC::td_net_block_len, NUM_NETS, long);
        CUDA_MALLOC_1D(MSC::td_net_block_pfx1, NUM_NETS, long);
        CUDA_MALLOC_1D(MSC::td_net_block_len1, NUM_NETS, long);
        if(BIT_PMODEL) {
            CUDA_MALLOC_1D(MSC::d_net_edges_gcell_cnt, NUM_NODES_LEN, gccnt_dt);
            CUDA_MALLOC_1D(MSC::d_net_edges_gcell_cnt_pfx, NUM_NODES_LEN, gccnt_dt);
        }
        CUDA_MALLOC_1D(MSC::td_net_mst_edges, NUM_NODES_LEN, nenc_dt);
        CUDA_MALLOC_1D(MSC::td_mst_lpat_pat, NUM_NODES_LEN, gvia_dt);
        MSC::net_mst_edges.clear(); MSC::net_mst_edges.resize(NUM_NODES_LEN);
        MSC::mst_lpat_pat.clear(); MSC::mst_lpat_pat.resize(NUM_NODES_LEN);
        MSC::net_block_pfx.clear(); MSC::net_block_pfx1.clear(); MSC::net_block_len.clear(); MSC::net_block_len1.clear();
        MSC::net_nodes_cnt.clear(); MSC::net_nodes_cnt_pfx.clear();
        MSC::net_block_pfx.resize(NUM_NETS); MSC::net_block_pfx1.resize(NUM_NETS); MSC::net_block_len.resize(NUM_NETS); MSC::net_block_len1.resize(NUM_NETS);
        MSC::net_nodes_cnt.resize(NUM_NETS); MSC::net_nodes_cnt_pfx.resize(NUM_NETS);
        MSC::NUM_NODES_OFST = 0;
        for(int i = 0; i < NUM_NETS; i++) {
            net_ids[i] = i;
        }
        REM_NETS = NUM_NETS;
    }
    __global__ void get_gpu_net_edge_cnt(int bofst, int blen) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= blen) {
            return;
        }
        int deg = d_net_deg_data[d_npart_ids[bofst + tid]];
        MSC::td_net_nodes_cnt[bofst + tid] = deg;
        MSC::td_net_block_len[bofst + tid] = deg * deg;
        MSC::td_net_block_len1[bofst + tid] = ((deg * (deg - 1)) >> 1);
    }
    void get_bbx_data() {
        bx_arr = (bbx_t *)malloc(sizeof(bbx_t) * NUM_NETS);
        for(int i = 0; i < NUM_NETS; i++) {
            bx_arr[i] = nets[i].bx_data;
            bbx_t cbx_val = bx_arr[i];
            assert(cbx_val.minx <= cbx_val.maxx); assert(cbx_val.miny <= cbx_val.maxy); assert(cbx_val.minz <= cbx_val.maxz);
            assert(cbx_val.minx != MAXX); assert(cbx_val.miny != MAXY); assert(cbx_val.minz != MAXZ);
            assert(cbx_val.maxx != -1); assert(cbx_val.maxy != -1); assert(cbx_val.maxz != -1);
        }
    }
    void part_ydir_vec(vector<int> &pot_ids) {
        sort(pot_ids.begin(), pot_ids.end(), [](int a, int b) {
            if(bx_arr[a].miny == bx_arr[b].miny) {
                return bx_arr[a].maxy < bx_arr[b].maxy;
            }
            return bx_arr[a].miny < bx_arr[b].miny;
        });
        gcd_dt cmaxy = -1;
        for(int i = 0; i < pot_ids.size(); i++) {
            int cnet_id = pot_ids[i];
            if(bx_arr[cnet_id].miny > cmaxy || bx_arr[cnet_id].minx > CMAXX) {
                SPDLOG_ASSERT(vis_data[cnet_id] == 0, "");
                vis_data[cnet_id] = 1;
                cmaxy = max(cmaxy, bx_arr[cnet_id].maxy);
                RU_BIT->range_update(bx_arr[cnet_id].minx, bx_arr[cnet_id].maxx, 1);
                npart_len[PART_ID]++; npart_ids[GP_ID++] = cnet_id; TKN++;
            }
        }
    }
    void part_ydir(int st, int ed) {
        sort(net_ids + st, net_ids + ed, [](int a, int b) {
            if(bx_arr[a].miny == bx_arr[b].miny) {
                return bx_arr[a].maxy < bx_arr[b].maxy;
            }
            return bx_arr[a].miny < bx_arr[b].miny;
        });
        gcd_dt cmaxy = -1;
        for(int i = st; i < ed; i++) {
            int cnet_id = net_ids[i];
            if(bx_arr[cnet_id].miny > cmaxy || bx_arr[cnet_id].minx > CMAXX) {
                SPDLOG_ASSERT(vis_data[cnet_id] == 0, "");
                vis_data[cnet_id] = 1;
                cmaxy = max(cmaxy, bx_arr[cnet_id].maxy); CMAXX = max(CMAXX, bx_arr[cnet_id].maxx);
                npart_len[PART_ID]++; npart_ids[GP_ID++] = cnet_id; TKN++;
            }
        }
        sort(net_ids + st, net_ids + ed, [](int a, int b) {
            if(bx_arr[a].minx == bx_arr[b].minx) {
                return bx_arr[a].maxx < bx_arr[b].maxx;
            }
            return bx_arr[a].minx < bx_arr[b].minx;
        });
    }
    void add_net_edge_data(int net_id) {
        net_t cnet = nets[net_id];
        vector<pin_t> cpins = cnet.pins; EDGE_ID = 0;
        if(cpins.size() > MAX_EDGE_DATA_PINS) {
            spdlog::warn("skipping_net_add_edge net_id {} net_name {} num_pins {}", net_id, cnet.net_name, cpins.size());
            return;
        }
        for(int i = 0; i < cpins.size(); i++) {
            pin_t pin1 = cpins[i];
            genc_dt enc1 = UTIL::get_enc<WIRE_ENC>(pin1.pcd.z, pin1.pcd.x, pin1.pcd.y);
            for(int j = i + 1; j < cpins.size(); j++) {
                pin_t pin2 = cpins[j];
                genc_dt enc2 = UTIL::get_enc<WIRE_ENC>(pin2.pcd.z, pin2.pcd.x, pin2.pcd.y);
                // DB_INIT::add_query_edge_data(EDGE_QUERY, net_id, EDGE_ID, enc1, enc2);
                if(EDGE_ID && EDGE_ID % EDGE_QUERY_EXEC_THRESHOLD == 0) {
                    EDGE_QUERY.back() = ';';
                    // DB_INIT::batch_insert_query(EDGE_QUERY); EDGE_QUERY = DB_INIT::edge_data_ref_query;
                }
                EDGE_ID++;
            }
        }
    }
    void add_net_meta_data(int &batch_id) {
        // SPDLOG_ASSERT(NUM_NETS == nets.size(), "NUM_NETS = {} nets.size() = {}", NUM_NETS, nets.size());
        // int version = 0;
        // string res_query = DBI::net_metadata_ref_query;
        // int bofst = npart_ofst[batch_id], blen = npart_len[batch_id];
        // for(int i = bofst; i < bofst + blen; i++) {
        //     net_t cnet = nets[npart_ids[i]];
        //     int deg = cnet.pins.size();
        //     SPDLOG_ASSERT(cnet.net_name.size() <= 300, "actual len is {} not {}", cnet.net_name.size(), 300);
        //     DB_INIT::add_query_net_metadata(res_query, version, batch_id, npart_ids[i], cnet.net_name, deg, cnet.bx_data.minx, cnet.bx_data.maxx, cnet.bx_data.miny, cnet.bx_data.maxy, cnet.bx_data.minz, cnet.bx_data.maxz);
        // }
        // SPDLOG_ASSERT(res_query.size() > 0, "total size is {}", res_query.size());
        // res_query.back() = ';';
        // DBI::batch_insert_query(res_query);
    }
    bool is_rbt_bbx_intrsct(box_t crect) {
        vector<box_t> res;
        res.clear();
        net_bbx.query(bgi::intersects(crect), back_inserter(res));
        if(res.size() > 0) {
            return true;
        }
        return false;
    }
    bool is_rbt_bbx_cintrsct(bxtree_t &net_bbx1, box_t crect) {
        vector<box_t> res;
        res.clear();
        net_bbx1.query(bgi::intersects(crect), back_inserter(res));
        if(res.size() > 0) {
            return true;
        }
        return false;
    }
    void part_xdim() {
        int targ_num = (REM_NETS + NUM_XDIM_PARTS - 1) / NUM_XDIM_PARTS;
        int xmax_id = -1;
        auto cmp1 = [](int a, int b) {
            return bx_arr[a].maxx < bx_arr[b].maxx;
        };
        auto cmp2 = [](int a, int b) {
            gcd_dt a1 = (bx_arr[a].maxx - bx_arr[a].minx + 1) * (bx_arr[a].maxy - bx_arr[a].miny + 1);
            gcd_dt a2 = (bx_arr[b].maxx - bx_arr[b].minx + 1) * (bx_arr[b].maxy - bx_arr[b].miny + 1);
            return a1 < a2;
        };
        auto cmp3 = [](int a, int b) {
            gcd_dt a1 = (bx_arr[a].maxx - bx_arr[a].minx + 1) + (bx_arr[a].maxy - bx_arr[a].miny + 1);
            gcd_dt a2 = (bx_arr[b].maxx - bx_arr[b].minx + 1) + (bx_arr[b].maxy - bx_arr[b].miny + 1);
            return a1 < a2;
        };
        auto cmp4 = [](int a, int b) {
            return nets[a].pins.size() < nets[b].pins.size();
        };
        sort(net_ids, net_ids + REM_NETS, [&](int a, int b) {
                return cmp1(a, b);
        });
        vector<int> xid_net(REM_NETS, -1), part_ids_vec(NUM_XDIM_PARTS, -1);
        int cmax_id = -1, cur_part_id = 0, num_pnets = 0;
        vector<int> part_len(NUM_XDIM_PARTS, 0);
        for(int i = 0; i < REM_NETS; i++) {
            int cnet_id = net_ids[i], minx = bx_arr[cnet_id].minx, maxx = bx_arr[cnet_id].maxx;
            if(cmax_id < minx) {
                xid_net[i] = cur_part_id; num_pnets++;
                part_len[cur_part_id]++;
            }
            if(num_pnets == targ_num) {
                part_ids_vec[cur_part_id] = i + 1;
                num_pnets = 0; cur_part_id++; cmax_id = maxx;
            }
        }
        SPDLOG_ASSERT(cur_part_id < NUM_XDIM_PARTS, "{} {}", cur_part_id, NUM_XDIM_PARTS);
        vector<vector<int>> pick_ids(cur_part_id + 10);
        for(int i = 0; i <= cur_part_id; i++) {
            pick_ids[i] = vector<int>(part_len[i], 0);
        }
        for(int i = 0; i < REM_NETS; i++) {
            if(xid_net[i] != -1) {
                int prt_id = xid_net[i], fid = part_len[prt_id]--;
                pick_ids[prt_id][fid - 1] = net_ids[i];
            }
        }
        #pragma omp parallel for num_threads(cur_part_id)
        for(int i = 0; i <= cur_part_id; i++) {
            sort(pick_ids[i].begin(), pick_ids[i].end(), [&](int a, int b) {
                return cmp3(a, b);
            });
            bxtree_t net_bbx1 = bxtree_t();
            net_bbx1.clear();
            for(int pnet_id = 0; pnet_id < pick_ids[i].size(); pnet_id++) {
                int nid = pick_ids[i][pnet_id];
                SPDLOG_ASSERT(vis_data[nid] == 0, "");
                box_t cbx = box_t(point_t(bx_arr[nid].minx, bx_arr[nid].miny), point_t(bx_arr[nid].maxx, bx_arr[nid].maxy));
                if(!is_rbt_bbx_cintrsct(net_bbx1, cbx)) {
                    net_bbx1.insert(cbx);
                    #pragma omp critical
                    {
                        npart_len[PART_ID]++; npart_ids[GP_ID++] = nid; TKN++;
                    }
                    vis_data[nid] = 1;
                }
            }
        }
    }
    void part_nets() {
        TIME_();
        spdlog::info("started net partition");
        init_data();
        get_bbx_data(); PART_ID = 0;
        #if ADD_PART_DATA_DB == 1
            long *time_arr = (long *)malloc(sizeof(long) * NUM_NETS);
        #endif
        auto cmp1 = [](int a, int b) {
            if(bx_arr[a].minx == bx_arr[b].minx) {
                return bx_arr[a].maxx < bx_arr[b].maxx;
            }
            return bx_arr[a].minx < bx_arr[b].minx;
        };
        auto cmp2 = [](int a, int b) {
            gcd_dt a1 = (bx_arr[a].maxx - bx_arr[a].minx + 1) * (bx_arr[a].maxy - bx_arr[a].miny + 1);
            gcd_dt a2 = (bx_arr[b].maxx - bx_arr[b].minx + 1) * (bx_arr[b].maxy - bx_arr[b].miny + 1);
            return a1 < a2;
        };
        if(RBT_PART || ORD_PART) {
            sort(net_ids, net_ids + REM_NETS, [&](int a, int b) {
                return cmp2(a, b);
            });
        } else {
            sort(net_ids, net_ids + REM_NETS, [&](int a, int b) {
                return cmp1(a, b);
            });
        }
        GP_ID = 0; PART_OFST = 0;
        int prev_cnt_ofst = 0;
        if(ORD_PART) {
            RU_BIT = make_unique<ru_bit_tree>(X);
        }
        while(REM_NETS > 0) {
            #if ADD_PART_DATA_DB == 1
                auto start_time = chrono::high_resolution_clock::now();
            #endif
            // spdlog::info("REMAINING NETS: {}", REM_NETS);
            //sort them based on bbx x coordinate data
            int pnet_id = 0;
            TKN = 0; CMAXX = -1;
            if(RBT_PART) {
                net_bbx = bxtree_t();
                net_bbx.clear();
                while(pnet_id < REM_NETS) {
                    int nid = net_ids[pnet_id];
                    SPDLOG_ASSERT(vis_data[nid] == 0, "");
                    box_t cbx = box_t(point_t(bx_arr[nid].minx, bx_arr[nid].miny), point_t(bx_arr[nid].maxx, bx_arr[nid].maxy));
                    if(!is_rbt_bbx_intrsct(cbx)) {
                        // net_bbx.insert(box_t(point_t(bx_arr[nid].minx, bx_arr[nid].miny), point_t(bx_arr[nid].maxx, bx_arr[nid].maxy)));
                        net_bbx.insert(cbx);
                        npart_len[PART_ID]++; npart_ids[GP_ID++] = nid; TKN++;
                        vis_data[nid] = 1;
                    }
                    pnet_id++;
                }
            } else if(ORD_PART) {
                vector<int> pot_ids;
                RU_BIT->reset();
                while(pnet_id < REM_NETS) {
                    pot_ids.clear();
                    while(pnet_id < REM_NETS) {
                        int nid = net_ids[pnet_id];
                        SPDLOG_ASSERT(vis_data[nid] == 0, "");
                        int net_cnt = RU_BIT->range_query(bx_arr[nid].minx, bx_arr[nid].maxx);
                        if(net_cnt == 0 && pot_ids.size() < XDIM_NETS) {
                            pot_ids.push_back(nid);
                        }
                        pnet_id++;
                    }
                    part_ydir_vec(pot_ids);
                }
            } else if(BATCH_PART) {
                part_xdim();
            } else {
                while(pnet_id < REM_NETS) {
                    while(pnet_id < REM_NETS) {
                        int nid = net_ids[pnet_id];
                        // cmaxx = max(cmaxx, bx_arr[])
                        SPDLOG_ASSERT(vis_data[nid] == 0, "");
                        if(bx_arr[nid].minx > CMAXX) {
                            break;
                        }
                        pnet_id++;
                    }
                    part_ydir(pnet_id, min(pnet_id + XDIM_NETS, REM_NETS));
                    pnet_id = min(pnet_id + XDIM_NETS, REM_NETS);
                }
            }
            //sort based on vis_data
            int id = 0, fid = 0;
            while(id < REM_NETS) {
                if(vis_data[net_ids[fid]] == 0) {
                    fid++;
                }
                if(id > fid && vis_data[net_ids[id]] == 0) {
                    SPDLOG_ASSERT(vis_data[net_ids[fid]] == 1, "vis data is not as expected");
                    swap(net_ids[fid], net_ids[id]);
                    fid++;
                }
                id++;
            }
            #if ADD_PART_DATA_DB == 1
            auto end_time = chrono::high_resolution_clock::now();
            long total_time_taken = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
            time_arr[PART_ID] = total_time_taken;
            #endif
            REM_NETS -= TKN;
            // spdlog::info("GOT {} NETS", TKN);
            PART_OFST += prev_cnt_ofst; prev_cnt_ofst = npart_len[PART_ID]; npart_ofst[PART_ID] = PART_OFST;
            SPDLOG_ASSERT(prev_cnt_ofst <= MAX_BATCH_SIZE, "{} {}", prev_cnt_ofst, MAX_BATCH_SIZE);
            sort(npart_ids + PART_OFST, npart_ids + PART_OFST + prev_cnt_ofst, [](int a, int b) {
               return  nets[a].pins.size() < nets[b].pins.size();
            });
            CUDAMEMCPYHD_1D(d_npart_ids + PART_OFST, npart_ids + PART_OFST, prev_cnt_ofst, int);
            int num_blks = CEIL_DIV(prev_cnt_ofst, NUM_GPU_THREADS);
            get_gpu_net_edge_cnt<<<num_blks, NUM_GPU_THREADS>>>(PART_OFST, prev_cnt_ofst);
            cudaDeviceSynchronize();
            thrust::device_ptr<long> d_edge_cnt_ptr = thrust::device_pointer_cast(MSC::td_net_block_len);
            thrust::device_ptr<long> d_edge_cnt_optr = thrust::device_pointer_cast(MSC::td_net_block_pfx);
            thrust::device_ptr<long> d_edge_cnt_ptr1 = thrust::device_pointer_cast(MSC::td_net_block_len1);
            thrust::device_ptr<long> d_edge_cnt_optr1 = thrust::device_pointer_cast(MSC::td_net_block_pfx1);
            thrust::device_ptr<int> d_net_nodes_cnt_ptr = thrust::device_pointer_cast(MSC::td_net_nodes_cnt);
            thrust::device_ptr<int> d_net_nodes_cnt_optr = thrust::device_pointer_cast(MSC::td_net_nodes_cnt_pfx);
            thrust::exclusive_scan(d_edge_cnt_ptr + PART_OFST, d_edge_cnt_ptr + PART_OFST + prev_cnt_ofst, d_edge_cnt_optr + PART_OFST);
            thrust::exclusive_scan(d_edge_cnt_ptr1 + PART_OFST, d_edge_cnt_ptr1 + PART_OFST + prev_cnt_ofst, d_edge_cnt_optr1 + PART_OFST);
            thrust::exclusive_scan(d_net_nodes_cnt_ptr + PART_OFST, d_net_nodes_cnt_ptr + PART_OFST + prev_cnt_ofst, d_net_nodes_cnt_optr + PART_OFST);
            cudaDeviceSynchronize();
            CUDAMEMCPYDH_1D(MSC::net_nodes_cnt.data() + PART_OFST, MSC::td_net_nodes_cnt + PART_OFST, prev_cnt_ofst, int);
            CUDAMEMCPYDH_1D(MSC::net_nodes_cnt_pfx.data() + PART_OFST, MSC::td_net_nodes_cnt_pfx + PART_OFST, prev_cnt_ofst, int);
            CUDAMEMCPYDH_1D(MSC::net_block_len.data() + PART_OFST, MSC::td_net_block_len + PART_OFST, prev_cnt_ofst, long);
            CUDAMEMCPYDH_1D(MSC::net_block_len1.data() + PART_OFST, MSC::td_net_block_len1 + PART_OFST, prev_cnt_ofst, long);
            CUDAMEMCPYDH_1D(MSC::net_block_pfx.data() + PART_OFST, MSC::td_net_block_pfx + PART_OFST, prev_cnt_ofst, long);
            CUDAMEMCPYDH_1D(MSC::net_block_pfx1.data() + PART_OFST, MSC::td_net_block_pfx1 + PART_OFST, prev_cnt_ofst, long);
            #if TEST_HANAN_STPTS_GEN == 1
            // spdlog::info("completed batch_id {}", PART_ID);
            HNS::get_stpts(PART_OFST, prev_cnt_ofst);
            #endif
            #if ADD_NET_METADATA_DB == 1
            // spdlog::info("adding net metadata to db");
            add_net_meta_data(PART_ID);
            #endif
            PART_ID++;
            net_gen_bsem.release();
        }
        spdlog::info("number of partitions {}", PART_ID);
        #if NET_PART_CHECK == 1
        spdlog::warn("testing net partition");
        if(XDIM_NETS % TEST_DNM == 0) PART_UTIL::test_part(bx_arr, PART_ID);
        #endif
        net_gen_bsem.release();
        #if ADD_PART_DATA_DB == 1
        for(int i = 0; i < PART_ID; i++) {
            int version = 0;
            long num_prs;
            UTIL::get_num_pts_batch(num_prs, i);
            DB_INIT::insert_part_data(version, i, npart_len[i], num_prs, time_arr[i]);
        }
        free(time_arr);
        #endif
        #if ADD_EDGE_DATA == 1
            spdlog::warn("adding edge_data");
            EDGE_ID = 0; EDGE_QUERY.clear(); EDGE_QUERY = DB_INIT::edge_data_ref_query;
            for(int i = 0; i < PART_ID; i++) {
                long ofst = npart_ofst[i], len = npart_len[i];
                for(int j = ofst; j < ofst + len; j++) {
                    int cnet_id = net_ids[j];
                    add_net_edge_data(cnet_id);
                }
            }
            if(EDGE_QUERY.back() == ',') EDGE_QUERY.back() = ';';
            DB_INIT::batch_insert_query(EDGE_QUERY);
            spdlog::warn("done adding edge_data");
            EDGE_QUERY.clear();
        #endif
        free(vis_data); 
        // proc_bsem.acquire();
        // GEN_STATE = false;
    }
}
