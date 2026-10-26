#pragma once
#include "../common.h"
#include "./lpattern.hpp"
#include "./lpattern_util.hpp"
#include "./lpattern_common.h"
#include "../util.hpp"
#include "../mst_files/mst_common.h"
#include "../mst_files/mst_build.hpp"
#include "../brvsk_files/brvsk_common.h"
#include "../brvsk_files/brvsk_build.hpp"
#include "../bit_util.hpp"

using namespace std;
using namespace COMMON_DATA;
namespace LPC = LPAT_COMMON;
namespace LPU = LPAT_UTIL;
namespace MSB = MST_BUILD;
namespace BSB = BRVSK_BUILD;
namespace BULS = BIT_UTIL_SPACE;

namespace TEST_LPAT {
    void test_gpu_get_pfx() {
        wcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE); 
        vcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE); 
        hd_wcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE);
        CUDA_MALLOC_1D(d_wcost_cube, CSIZE, gcb_dt);
        CUDA_MALLOC_1D(d_vcost_cube, CSIZE, gcb_dt);
        UTIL::get_cost_data();
        CUDAMEMCPYDH_1D(wcost_cube, d_wcost_cube, CSIZE, gcb_dt);
        CUDAMEMCPYDH_1D(vcost_cube, d_vcost_cube, CSIZE, gcb_dt);
        // UTIL::get_ofcost();
        UTIL::pfx_cube_cost();
        string pfx_res = "";
        auto get_thrust_res = [&pfx_res]() {
            auto start_time = chrono::high_resolution_clock::now();
            thrust::device_ptr<gcb_dt> pfx_ptr = thrust::device_pointer_cast(d_vcost_cube);
            thrust::exclusive_scan(pfx_ptr, pfx_ptr + CSIZE, pfx_ptr);
            cudaDeviceSynchronize();
            auto end_time = chrono::high_resolution_clock::now();
            auto total_time_taken = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
            pfx_res += "(thrust," + to_string(total_time_taken) + "us),";
            spdlog::info("total time taken by thrust is {}us", total_time_taken);
        };
        get_thrust_res();
        auto get_num_blocks = []() {
            int res = 0;
            for(int i = 0; i < Z; i++) {
                res += (i & 1) ? X : Y;
            }
            return res;
        };
        auto check_arr_similar = []() {
            for(int i = 0; i < CSIZE; i++) {
                if(abs(hd_wcost_cube[i] - wcost_cube[i]) > diff_tolerance) {
                    gcd_dt x1, y1;
                    gvia_dt z1;
                    UTIL::get_dec_all<WIRE_ENC>(i, z1, x1, y1);
                    spdlog::info("failed for [z.x.y] {} {} {} enc {}", z1, x1, y1, i);
                    SPDLOG_ASSERT(abs(hd_wcost_cube[i] - wcost_cube[i]) <= diff_tolerance, "mismatch at index {} vals {} {}", i, hd_wcost_cube[i], wcost_cube[i]);
                }
            }
        };
        int num_blocks = get_num_blocks();
        for(int max_mem = 500; max_mem <= 2500; max_mem += 100) {
            cudaFuncSetAttribute(UTIL::d_get_pfx_cost_cube, cudaFuncAttributeMaxDynamicSharedMemorySize, max_mem * sizeof(gcb_dt));
            #if TEST_GPU_PFX_LOGS == 1
            spdlog::info("max_mem {}", max_mem);
            #endif
            string temp_res = "";
            for(int num_ths = 8; num_ths <= NUM_GPU_THREADS; num_ths += 8) {
                #if TEST_GPU_PFX_LOGS_DETAIL == 1
                spdlog::info("[max_mem, num_ths] [{},{}]", max_mem, num_ths);
                #endif
                auto start_time = chrono::high_resolution_clock::now();
                // spdlog::info("max_mem {} num_ths {}", max_mem, num_ths);
                UTIL::get_cost_data();
                UTIL::d_get_epfx_cost_cube_inplace<<<num_blocks, NUM_GPU_THREADS, sizeof(gcb_dt) * max_mem>>>(d_wcost_cube, num_ths, max_mem);
                int nnm_blocks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
                cudaDeviceSynchronize();
                auto end_time = chrono::high_resolution_clock::now();
                auto total_time_taken = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
                // spdlog::info("total time taken is {}us", total_time_taken);
                CUDAMEMCPYDH_1D(hd_wcost_cube, d_wcost_cube, CSIZE, gcb_dt);
                temp_res += "(" + to_string(max_mem) + "," + to_string(num_ths) + "," + to_string(total_time_taken) + "us),";
                check_arr_similar();
            }
            pfx_res += temp_res;
        }
        cudaFree(d_wcost_cube); cudaFree(d_vcost_cube); 
        free(hd_wcost_cube); free(wcost_cube); free(vcost_cube);
        spdlog::info("pfx_time_data {}", pfx_res);
    }
    void test_gpu_get_epfx_inplace() {
        wcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE);
        vcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE);
        hd_wcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE); 
        CUDA_MALLOC_1D(d_wcost_cube, CSIZE, gcb_dt);
        CUDA_MALLOC_1D(d_vcost_cube, CSIZE, gcb_dt);
        int nm_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
        UTIL::get_cost_data();
        cudaDeviceSynchronize();
        CUDAMEMCPYDH_1D(wcost_cube, d_wcost_cube, CSIZE, gcb_dt);
        // UTIL::get_ofcost();
        UTIL::pfx_cube_cost();        
        auto get_num_blocks = []() {
            int res = 0;
            for(int i = 0; i < Z; i++) {
                res += (i & 1) ? X : Y;
            }
            return res;
        };
        auto check_arr_similar = []() {
            for(int i = 0; i < CSIZE; i++) {
                if(abs(hd_wcost_cube[i] - wcost_cube[i]) > diff_tolerance) {
                    gcd_dt x1, y1;
                    gvia_dt z1;
                    UTIL::get_dec_all<WIRE_ENC>(i, z1, x1, y1);
                    spdlog::info("failed for [z.x.y] {} {} {} enc {}", z1, x1, y1, i);
                    SPDLOG_ASSERT(abs(hd_wcost_cube[i] - wcost_cube[i]) <= diff_tolerance, "mismatch at index {} vals {} {}", i, hd_wcost_cube[i], wcost_cube[i]);
                }
            }
        };
        int num_blocks = get_num_blocks();
        int max_mem = max(X, Y), num_ths = 128;
        // cudaFuncSetAttribute(LPU::d_get_epfx_cost_cube_inplace, cudaFuncAttributeMaxDynamicSharedMemorySize, max_mem * sizeof(gcb_dt));
        UTIL::d_get_epfx_cost_cube_inplace<<<num_blocks, NUM_GPU_THREADS, sizeof(gcb_dt) * max_mem>>>(d_wcost_cube, num_ths, max_mem);
        CUDAMEMCPYDH_1D(hd_wcost_cube, d_wcost_cube, CSIZE, gcb_dt);
        check_arr_similar();
        cudaFree(d_wcost_cube); cudaFree(d_vcost_cube);
        free(hd_wcost_cube); free(wcost_cube); free(vcost_cube);
    }
    // void test_lpat() {
        // SPDLOG_ASSERT(VIA_WT == WIRE_WT, "VIA_WT {} WIRE_WT {}", VIA_WT, WIRE_WT);
        // int batch_id = 0;
        // wcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE); 
        // vcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE);
        // LPC::lpat_dist = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE);
        // LPC::lpat_par = (int *)malloc(sizeof(int) * CSIZE);
        // DJC::djkst_oval.clear(); DJC::djkst_par.clear();
        // DJC::djkst_oval = vector<vector<gcb_dt>>(CSIZE, vector<gcb_dt>(MAX_DJKST_BENDS + 1, MAX_OFCOST));
        // DJC::maze_val.clear(); 
        // DJC::maze_val.resize(CSIZE);
        // DJC::maze_par.clear();
        // DJC::maze_par.resize(CSIZE);
        // UTIL::get_cost_data();
        // CUDAMEMCPYDH_1D(wcost_cube, d_wcost_cube, CSIZE, gcb_dt);
        // CUDAMEMCPYDH_1D(vcost_cube, d_vcost_cube, CSIZE, gcb_dt);
        // UTIL::pfx_cube_cost();
        // net_gen_bsem.acquire();
        // while(batch_id < PART_ID) {
        //     long num_prs;
        //     #if RANDOM_GEN == 0
        //     LPAT_NUM_PTS = npart_len[batch_id];
        //     int ofst = npart_ofst[batch_id];
        //     UTIL::get_num_pts_batch(num_prs, batch_id);
        //     spdlog::info("number of nets, edges in batch {} is {}, {}", batch_id, LPAT_NUM_PTS, num_prs);
        //     rd_enc = (pt_prs *)malloc(sizeof(pt_prs) * num_prs);
        //     int pts_id = 0;
        //     for(int i = ofst; i < ofst + LPAT_NUM_PTS; i++) {
        //         int net_id = npart_ids[i];
        //         vector<pin_t> pins = nets[net_id].pins;
        //         for(int j = 0; j < pins.size(); j++) {
        //             coord_t ccd = pins[j].pcd;
        //             gcd_dt x1 = ccd.x, y1 = ccd.y;
        //             gvia_dt z1 = ccd.z;
        //             genc_dt enc1 = UTIL::get_enc<CUBE_ENC>(z1, x1, y1);
        //             for(int k = j + 1; k < pins.size(); k++) {
        //                 coord_t ncd = pins[k].pcd;
        //                 gcd_dt x2 = ncd.x, y2 = ncd.y;
        //                 gvia_dt z2 = ncd.z;
        //                 genc_dt enc2 = UTIL::get_enc<CUBE_ENC>(z2, x2, y2);
        //                 rd_enc[pts_id++] = {enc1:enc1, enc2:enc2};
        //             }
        //         }
        //     }
        //     SPDLOG_ASSERT(pts_id == num_prs, "pts_id {} num_prs {} are not matching", pts_id, num_prs);
        //     #else
        //     LPAT_NUM_PTS = RGEN_SPACE::rlen[batch_id];
        //     int ofst = RGEN_SPACE::rofst[batch_id];
        //     num_prs = long(LPAT_NUM_PTS);
        //     spdlog::info("number of edges in batch {} is {}", batch_id, num_prs);
        //     rd_enc = (pt_prs *)malloc(sizeof(pt_prs) * num_prs);
        //     int pts_id = 0;
        //     for(int i = ofst; i < ofst + num_prs; i++) {
        //         rd_enc[pts_id++] = RGEN_SPACE::r_pair[i];
        //     }
        //     SPDLOG_ASSERT(pts_id == num_prs, "pts_id {} num_prs {} are not matching", pts_id, num_prs);
        //     #endif
        //     for(int i = 0; i < num_prs; i++) {
        //         // if(i % 50000 == 0 && i) {
        //         //     spdlog::info("completed {} edges", i);
        //         // }
        //         pt_prs cpr = rd_enc[i];
        //         genc_dt enc1 = cpr.enc1, enc2 = cpr.enc2;
        //         gcd_dt x1, y1, x2, y2;
        //         gvia_dt z1, z2;
        //         // enc1 = UTIL::get_enc(2,0,427); enc2 = UTIL::get_enc(0,1,435);
        //         UTIL::get_dec_all<CUBE_ENC>(enc1, z1, x1, y1);
        //         UTIL::get_dec_all<CUBE_ENC>(enc2, z2, x2, y2);
        //         if(z1 == z2 && y1 == y2 && x1 == x2) {
        //             continue;
        //         }
        //         if(abs(x1 - x2) * abs(y1 - y2) > 5000) {
        //             continue;
        //         }
        //         gcb_dt res1, res2;
        //         lpat_dt pat;
        //         // spdlog::info("checking for [z.x.y] ({},{},{}) ({},{},{})", z1, x1, y1, z2, x2, y2);
        //         LPATTERN::get_lpat_dist(enc1, enc2, res1, pat);
        //         // DJ::get_djkst_val(enc1, enc2, res1);
        //         DJ::get_maze_val(enc1, enc2, res2);
        //         gvia_dt z3, z4;
        //         LPU::get_zdec(pat, z3, z4);
        //         if(abs(res2 - res1) > diff_tolerance) {
        //             DJU::get_maze_path(enc1, enc2);
        //             // LPATTERN::get_lpat_djkst(enc2, enc1, res2);
        //             // if(abs(res1 - res2) > diff_tolerance) {
        //             //     LPAT_UTIL::get_djkst_path(enc2, enc1);
        //             //     LPAT_UTIL::print_lpat_path(enc1, enc2, pat);
        //             // }
        //             SPDLOG_ASSERT(abs(res2 - res1) <= diff_tolerance, "index {} distance for points [z,x,y] ({},{},{}) ({},{},{}) {} {} {} {}", i, z1, x1, y1, z2, x2, y2, res1, res2, z3, z4);
        //         }                
        //     }
        //     free(rd_enc);
        //     batch_id++;
        //     net_gen_bsem.acquire();
        // }
        // DJC::djkst_oval.clear(); DJC::djkst_par.clear();
        // spdlog::info("total number of batches are {}", batch_id);
    // }
    void test_gpu_lpat() {
        // SPDLOG_ASSERT(VIA_WT == WIRE_WT, "VIA_WT {} WIRE_WT {}", VIA_WT, WIRE_WT);
        UTIL::print_gpu_free_mem();
        int batch_id = 0, version_id = 1;
        int nm_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
        CUDA_MALLOC_1D(d_wcost_cube, CSIZE, gcb_dt);
        CUDA_MALLOC_1D(d_vcost_cube, CSIZE, gcb_dt);
        UTIL::get_cost_data();
        UTIL::init_bit_cost();
        if(BIT_PMODEL) {
            BULS::init_bit_data();
        }
        cudaDeviceSynchronize();
        #if LCPU_GPU_COMPARISON == 1 || TEST_LPATTERN_GEN == 1
        wcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE); 
        vcost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE); 
        CUDAMEMCPYDH_1D(wcost_cube, d_wcost_cube, CSIZE, gcb_dt);
        CUDAMEMCPYDH_1D(vcost_cube, d_vcost_cube, CSIZE, gcb_dt);
        UTIL::pfx_cube_cost();
        #endif
        #if VIA_MINOR == 1
            via_cost_cube = (gcb_dt *)malloc(sizeof(gcb_dt) * CSIZE);
            CUDA_MALLOC_1D(d_via_cost_cube, CSIZE, gcb_dt);
        #endif
        #if VIA_MINOR == 1
        CUDAMEMCPYHD_1D(d_via_cost_cube, via_cost_cube, CSIZE, gcb_dt);
        #endif
        #if TEST_LPAT_ADD_EDGE == 1
        CUDA_MALLOC_1D(LPC::d_pth_mat, CSIZE, pth_dt);
        #if LCPU_GPU_EDGEADD_COMPARISON == 1
        LPC::pth_mat = (in_t *)malloc(sizeof(in_t) * CSIZE);
        LPC::hd_pth_mat = (pth_dt *)malloc(sizeof(pth_dt) * CSIZE);
        #endif
        #endif
        auto get_num_pfx_blocks = []() {
            int res = 0;
            for(int i = 0; i < Z; i++) {
                res += (i & 1) ? X : Y;
            }
            return res;
        };
        int memsize = max(X, Y), nth = 128, num_pfx_blocks = get_num_pfx_blocks();
        UTIL::d_get_epfx_cost_cube_inplace<<<num_pfx_blocks, NUM_GPU_THREADS, sizeof(gcb_dt) * memsize>>>(d_wcost_cube, nth, memsize);
        cudaDeviceSynchronize();
        net_gen_bsem.acquire();
        #if ADD_LPATTERN_BATCH_DIST == 1
        LPC::LPAT_BATCH_QUERY.clear();
        LPC::LPAT_BATCH_QUERY = DB_INIT::lpat_dist_ref_query;
        #endif
        while(batch_id < PART_ID) {
            if(!BIT_PMODEL) {
                UTIL::get_cost_data();
                cudaDeviceSynchronize();
            }
            int memsize = max(X, Y), nth = 128, num_pfx_blocks = get_num_pfx_blocks();
            if(!BIT_PMODEL) {
                UTIL::d_get_epfx_cost_cube_inplace<<<num_pfx_blocks, NUM_GPU_THREADS, sizeof(gcb_dt) * memsize>>>(d_wcost_cube, nth, memsize);
                cudaDeviceSynchronize();
            }            
            #if ADD_LPATTERN_BATCH_DIST == 1
            auto start_time_batch_perf = chrono::high_resolution_clock::now();
            #endif
            long num_prs;
            #if RANDOM_GEN == 0
            int num_nets = npart_len[batch_id], ofst = npart_ofst[batch_id];
            UTIL::get_num_pts_batch(num_prs, batch_id);
            #else
            int num_nets = RGEN_SPACE::rlen[batch_id], ofst = RGEN_SPACE::rofst[batch_id];
            num_prs = long(num_nets);
            #endif
            #if LPAT_LOGS == 1 
            spdlog::info("number of edges in batch {} is {}", batch_id, num_prs);
            #endif
            LPC::hd_ptprs = (genc_dt *)malloc(sizeof(genc_dt) * num_prs * 2);
            rd_enc = (pt_prs *)malloc(sizeof(pt_prs) * num_prs);
            CUDA_MALLOC_1D(LPC::d_ptprs, num_prs * 2, genc_dt);
            CUDA_MALLOC_1D(LPC::d_lpat_res, num_prs, gcb_dt);
            CUDA_MALLOC_1D(LPC::d_lpat_pat, num_prs, lpat_dt);
            // lpat_pat = (lpat_dt *)malloc(sizeof(lpat_dt) * num_prs);
            // TOT_LPAT_PRS = num_prs;
            #if RANDOM_GEN == 0
            int pts_id = 0;
            for(int i = ofst; i < ofst + num_nets; i++) {
                int net_id = npart_ids[i];
                vector<pin_t> pins = nets[net_id].pins;
                for(int j = 0; j < pins.size(); j++) {
                    coord_t ccd = pins[j].pcd;
                    gcd_dt x1 = ccd.x, y1 = ccd.y;
                    gvia_dt z1 = ccd.z;
                    genc_dt enc1 = UTIL::get_enc<CUBE_ENC>(z1, x1, y1);
                    #if NET_ID_ORDER_REV == 1
                    for(int k = 0; k < j; k++) {
                    #else
                    for(int k = j + 1; k < pins.size(); k++) {
                    #endif
                        coord_t ncd = pins[k].pcd;
                        gcd_dt x2 = ncd.x, y2 = ncd.y;
                        gvia_dt z2 = ncd.z;
                        genc_dt enc2 = UTIL::get_enc<CUBE_ENC>(z2, x2, y2);
                        rd_enc[pts_id++] = {enc1:enc1, enc2:enc2};
                    }
                }
            }
            SPDLOG_ASSERT(pts_id == num_prs, "pts_id {} num_prs {} are not matching", pts_id, num_prs);
            #else
            int pts_id = 0;
            for(int i = ofst; i < ofst + num_nets; i++) {
                rd_enc[pts_id++] = RGEN_SPACE::r_pair[i];
            }
            SPDLOG_ASSERT(pts_id == num_prs, "pts_id {} num_prs {} are not matching", pts_id, num_prs);
            #endif
            for(int i = 0; i < num_prs; i++) {
                pt_prs cpt = rd_enc[i];
                LPC::hd_ptprs[2 * i] = genc_dt(cpt.enc1); LPC::hd_ptprs[2 * i + 1] = genc_dt(cpt.enc2);
            }
            CUDAMEMCPYHD_1D(LPC::d_ptprs, LPC::hd_ptprs, num_prs * 2, genc_dt);
            #if LPAT_DETAIL_LOGS == 1
            spdlog::info("getting gpu lpatterns");
            #endif
            int num_pts_per_block = NUM_GPU_THREADS / Z1Z2, num_blocks = CEIL_DIV(num_prs, num_pts_per_block);
            #if LPAT_DETAIL_LOGS == 1
            spdlog::info("num_blocks for gpu L pattern {}", num_blocks);
            #endif
            SPDLOG_ASSERT(num_pts_per_block == 5, "num_pts_per_block {}", num_pts_per_block);
            auto start_time = chrono::high_resolution_clock::now();
            #if LPAT_DETAIL_LOGS == 1
            spdlog::info("invoking gpu now");
            #endif
            if(BIT_PMODEL) {
                BULS::get_gpu_lpat_dist<5><<<num_blocks, NUM_GPU_THREADS>>>(LPC::d_lpat_pat, LPC::d_ptprs, LPC::d_lpat_res, num_prs);
            } else {
                LPATTERN::get_gpu_lpat_dist<5><<<num_blocks, NUM_GPU_THREADS>>>(LPC::d_lpat_pat, LPC::d_ptprs, LPC::d_lpat_res, d_wcost_cube, d_vcost_cube, num_prs);
            }
            cudaDeviceSynchronize();
            auto end_time = chrono::high_resolution_clock::now();
            auto total_time_taken = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
            #if LPAT_DETAIL_LOGS == 1
            spdlog::info("num_pts {} total_time_taken_by_gpulpat {}us", num_prs, total_time_taken);
            #endif
            #if ADD_LPAT_PERF_DATA_DB == 1
            int num_ptss = num_prs;
            DB_INIT::insert_lpat_perf(version_id, batch_id, num_ptss, total_time_taken);
            #endif
            #if ADD_RANDOM_LPAT_PERF_DB == 1
            int num_ptss = num_prs;
            DB_INIT::insert_random_lpat_perf(batch_id, version_id, total_time_taken, num_ptss);
            #endif
            #if TEST_LPAT_ADD_EDGE == 1
            #if MST_TREE == 1
            MSB::get_mst_tree<0>(batch_id, num_prs);
            #endif
            #if BRVSK_TREE == 1
            BSB::get_mst_tree(batch_id, num_prs);
            #endif
            #if LCPU_GPU_EDGEADD_COMPARISON == 1
            // LPC::hd_lpat_pat = (lpat_dt *)malloc(sizeof(lpat_dt) * num_prs);
            CUDAMEMCPYDH_1D(LPC::hd_pth_mat, LPC::d_pth_mat, CSIZE, pth_dt);
            memset(LPC::pth_mat, 0, sizeof(in_t) * CSIZE);
            LPU::add_edge_update(batch_id);
            for(int i = 0; i < CSIZE; i++) {
                SPDLOG_ASSERT(LPC::hd_pth_mat[i] == LPC::pth_mat[i], "idx {} [hd_pth_mat[i],pth_mat[i]] [{},{}]", i, LPC::hd_pth_mat[i], LPC::pth_mat[i]);
            }                
            #endif
            free(LPC::hd_lpat_pat);
            #endif
            #if LCPU_GPU_COMPARISON == 1
            hd_lpat_res = (gcb_dt *)malloc(sizeof(gcb_dt) * num_prs);
            CUDAMEMCPYDH_1D(hd_lpat_res, LPC::d_lpat_res, num_prs, gcb_dt);
            start_time = chrono::high_resolution_clock::now();
            for(int i = 0; i < num_prs; i++) {
                pt_prs cpr = rd_enc[i];
                gcb_dt res1;
                lpat_dt pat;
                genc_dt enc1 = cpr.enc1, enc2 = cpr.enc2;            
                LPATTERN::get_lpat_dist(enc1, enc2, res1, pat);
                gcb_dt res2 = gcb_dt(hd_lpat_res[i]);
                gcd_dt x1, x2, y1, y2;
                gvia_dt z1, z2;
                if(abs(res2 - res1) > diff_tolerance) {
                    UTIL::get_dec_all<VIA_ENC>(enc1, z1, x1, y1);
                    UTIL::get_dec_all<VIA_ENC>(enc2, z2, x2, y2);
                    SPDLOG_ASSERT(abs(res2 - res1) <= diff_tolerance, "index {} distance for points [z,x,y] ({},{},{}) ({},{},{}) {} {}", i, z1, x1, y1, z2, x2, y2, res1, res2);
                }
            }
            end_time = chrono::high_resolution_clock::now();
            total_time_taken = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count();
            spdlog::info("num_pts {} total_time_taken_by_cpulpat {}us", num_prs, total_time_taken);
            free(LPC::hd_lpat_res);
            #endif
            #if ADD_LPATTERN_BATCH_DIST == 1
                auto end_time_batch_perf = chrono::high_resolution_clock::now();
                total_time_taken = chrono::duration_cast<chrono::microseconds>(end_time_batch_perf - start_time_batch_perf).count();
                CUDAMEMCPYDH_1D(hd_lpat_res, LPC::d_lpat_res, num_prs, gcb_dt);
                gcb_dt res_sum = 0;
                for(long i = 0; i < num_prs; i++) {
                    res_sum += hd_lpat_res[i];
                }
                DB_INIT::add_query_lpat_data(LPC::LPAT_BATCH_QUERY, version_id, batch_id, res_sum, total_time_taken);
            #endif
            #if TEST_LPATTERN_GEN == 1
            hd_lpat_res = (gcb_dt *)malloc(sizeof(gcb_dt) * num_prs);
            LPC::hd_lpat_pat = (lpat_dt *)malloc(sizeof(lpat_dt) * num_prs);
            CUDAMEMCPYDH_1D(LPC::hd_lpat_pat, LPC::d_lpat_pat, num_prs, lpat_dt);
            CUDAMEMCPYDH_1D(hd_lpat_res, LPC::d_lpat_res, num_prs, gcb_dt);
            for(int i = 0; i < num_prs; i++) {
                pt_prs cpr = rd_enc[i];
                in_t enc1 = cpr.enc1, enc2 = cpr.enc2, x1, y1, z1, x2, y2, z2;
                UTIL::get_dec_all(enc1, z1, x1, y1); UTIL::get_dec_all(enc2, z2, x2, y2);
                lpat_dt cpat = LPC::hd_lpat_pat[i];
                gcb_dt ref_res = hd_lpat_res[i], res;
                LUL::get_lpat_res1(enc1, enc2, cpat, res);
                SPDLOG_ASSERT(abs(res - ref_res) <= diff_tolerance, "idx {} pat {} res {} ref_res {} [zxy] [{},{},{}] [{},{},{}]", i, cpat, res, ref_res, z1, x1, y1, z2, x2, y2);
            }
            free(LPC::hd_lpat_pat); free(hd_lpat_res);
            #endif
            free(rd_enc); free(LPC::hd_ptprs);
            cudaFree(LPC::d_ptprs); cudaFree(LPC::d_lpat_res); cudaFree(LPC::d_lpat_pat);
            batch_id++;
            net_gen_bsem.acquire();
        }
        write_route_bsem.release();
        #if TEST_LPAT_ADD_EDGE == 1
        cudaFree(LPC::d_pth_mat);
            #if LCPU_GPU_EDGEADD_COMPARISON == 1
            free(LPC::pth_mat); free(LPC::hd_pth_mat);
            #endif
        #endif
        #if ADD_LPATTERN_BATCH_DIST == 1
        LPC::LPAT_BATCH_QUERY.back() = ';';
        DB_INIT::batch_insert_query(LPC::LPAT_BATCH_QUERY);
        LPC::LPAT_BATCH_QUERY.clear();
        #endif
        spdlog::info("total number of batches are {}", batch_id);
        #if LCPU_GPU_COMPARISON == 1 || TEST_LPATTERN_GEN == 1
        free(wcost_cube); free(vcost_cube);
        #endif
        cudaFree(d_wcost_cube); cudaFree(d_vcost_cube); 
        if(BIT_PMODEL) {
            cudaFree(d_cube_data_bit); cudaFree(d_cost_data_bit); cudaFree(d_bit_vis_data);

        }
        #if VIA_MINOR == 1
        free(via_cost_cube);
        cudaFree(d_via_cost_cube);
        #endif
    }
}
