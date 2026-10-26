#include "common.h"
#include "read_data.hpp"
#include "./part_net/net_part.hpp"
#include "l_pattern_files/test_lpat.hpp"
#include "l_pattern_files/lpattern_util.hpp"
#include "mst_files/mst_common.h"
#include "writer_files/writer.hpp"

using namespace std;
namespace LPU = LPAT_UTIL;
namespace MSC = MST_COMMON;
namespace RW = ROUTE_WRITER;
namespace RWC = ROUTE_WRITER_COMMON;
// namespace ZPU = ZPAT_UTIL;

void parse_args(int argc, char* argv[]) {
  int i = 0;
  while (i < argc) {
    auto s = std::string(argv[i]);
    if (s == "-cap") {
      assert(i + 1 < argc);
      COMMON_DATA::CAP_FILE = string(argv[++i]);
    } else if (s == "-net") {
      assert(i + 1 < argc);
      COMMON_DATA::NET_FILE = string(argv[++i]);
    } else if (s == "-gr_out") {
      assert(i + 1 < argc);
      COMMON_DATA::OUT_FILE_PATH = string(argv[++i]);
    } else if(s == "-chip_name") {
      assert(i + 1 < argc);
      COMMON_DATA::CHIP_NAME = string(argv[++i]);
    } else if(s == "-new_fmt") {
        COMMON_DATA::NEW_FMT = true;
    } else if(s == "-mst_cpu_threads") {
        MSC::MST_THREADS = stoi(argv[++i]);
    } else if(s == "-bit_pmodel") {
      COMMON_DATA::BIT_PMODEL = true;
    } else if(s == "-rbt_part") {
      COMMON_DATA::RBT_PART = true;
    } else if(s == "-rseg_fmt") {
      COMMON_DATA::RSEG_FMT = true;
    } else if(s == "-ord_part") {
      COMMON_DATA::ORD_PART = true;
    } else if(s == "-batch_part") {
      COMMON_DATA::BATCH_PART = true;
    } else if(s == "-num_xdim_parts") {
      COMMON_DATA::NUM_XDIM_PARTS = stoi(argv[++i]);
    } else if(s == "-ispd24_new_nf_path") {
      COMMON_DATA::ISPD24_NEW_NF_PATH = argv[++i];
    } else if(s == "-upd_ispd24_net_file") {
      COMMON_DATA::UPD_ISPD24_NET_FILe = true;
    } else if(s == "-writer_threads") {
      RWC::NUM_THREADS = stoi(argv[++i]);
    } else if(s == "-wrt_mrg") {
      COMMON_DATA::WRT_MRG = true;
    } else if(s == "-oprd_ppatch") {
      COMMON_DATA::OPRD_PPATCH = true;
    }
    ++i;
  }
}

int main(int argc, char *argv[]) {
    TIME_();
    parse_args(argc, argv);
    bool truncate_file = true;
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(OUT_FILE_PATH, truncate_file);
    gr_file_logger = spdlog::basic_logger_mt("gr_file_mt_logger", OUT_FILE_PATH);
    gr_file_logger->set_pattern("%v");
    if(UPD_ISPD24_NET_FILe) {
      READ_DATA::read_write_ispd24_nets(COMMON_DATA::NET_FILE);
      return;
    }
    READ_DATA::read_cap_file(COMMON_DATA::CAP_FILE);
    READ_DATA::read_net_file(COMMON_DATA::NET_FILE);
    #if ADD_MST_PERF_ALL_BATCH_DATA == 1
    auto start_time = chrono::high_resolution_clock::now();
    #endif
    #if RANDOM_GEN == 0
    COMMON_DATA::XDIM_NETS = 50;
    thread t1(NET_PART::part_nets);
    #else
    thread t1(RGEN_SPACE::get_rand_pts);
    #endif
    thread t2(TEST_LPAT::test_gpu_lpat);
    #if WRITE_ROUTE == 1
    thread t3(RW::write_route);
    #endif
    t1.join();
    t2.join();
    #if WRITE_ROUTE == 1
    t3.join();
    #endif
    #if ADD_MST_PERF_ALL_BATCH_DATA == 1
    auto end_time = chrono::high_resolution_clock::now();
    long total_time_taken_ms = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
    #if MST_TREE_GPU == 1 && MST_BRVSK == 1
    DBI::insert_mst_perf_data(BRVSK_VERSION, -1, MSC::BLEN, total_time_taken_ms);
    #elif MST_TREE_GPU == 0
    DBI::insert_mst_perf_data(KRSK_CPU_VERSION, -1, MSC::BLEN, total_time_taken_ms);
    #endif
    #endif
    #if EVAL_OFCOST == 1
    vector<int> eval_wire_cnt; eval_wire_cnt.clear(); eval_wire_cnt.resize(CSIZE);
    vector<int> eval_via_cnt; eval_via_cnt.clear(); eval_via_cnt.resize(CSIZE);
    CUDAMEMCPYDH_1D(eval_via_cnt.data(), d_eval_via_cnt, CSIZE, int);
    CUDAMEMCPYDH_1D(eval_wire_cnt.data(), d_eval_wire_cnt, CSIZE, int);
    eval_dt total_ofcost = 0.0, total_via_cost = 0, total_wire_cost = 0;
    for(int z = 1; z < Z; z++) {
        eval_dt lwire_cnt = 0.0, lvia_cnt = 0.0, lofcost = 0.0;
        eval_dt max_layer_wire_cnt = 0, max_layer_via_cnt = 0, max_layer_cap = 0.0, max_overflow = 0.0;
        for(int x = 0; x < X; x++) {
            for(int y = 0; y < Y; y++) {
                eval_dt cwire_cnt = eval_wire_cnt[UTIL::get_enc<WIRE_ENC>(z, x, y)];
                eval_dt cvia_cnt = eval_via_cnt[UTIL::get_enc<WIRE_ENC>(z, x, y)];
                max_layer_wire_cnt = max(max_layer_wire_cnt, cwire_cnt);
                max_layer_via_cnt = max(max_layer_via_cnt, cvia_cnt);
                max_layer_cap = max(max_layer_cap, cube_data[UTIL::get_enc<WIRE_ENC>(z, x, y)]);
                eval_dt cdemand = 2 * cwire_cnt + cvia_cnt;
                eval_dt of_val = cdemand  - 2 * cube_data[UTIL::get_enc<WIRE_ENC>(z, x, y)];
                max_overflow = max(max_overflow, of_val / 2);
                eval_dt cofcost = expf(of_val * 0.25);
                lofcost += cofcost; lwire_cnt += cwire_cnt; lvia_cnt += cvia_cnt;
            }
        }
        total_ofcost += lofcost; total_wire_cost += lwire_cnt; total_via_cost += lvia_cnt;
        spdlog::info("[layer, max_wire_cnt, max_via_cnt, max_layer_cap, max_overflow] [{},{},{},{},{}]", z, max_layer_wire_cnt, max_layer_via_cnt / 2, max_layer_cap, max_overflow);
        // spdlog::info("layer {} wire_cnt {}", z, lwire_cnt);
        // spdlog::info("layer {} via_cnt {}", z, lvia_cnt);
        // spdlog::info("layer {} ofcost {}", z, lofcost);
    }
    spdlog::info("overall wire_cnt {}", total_wire_cost);
    spdlog::info("overall via_cnt {}", total_via_cost);
    spdlog::info("overall ofcost {}", total_ofcost);
    #endif
    return 0;
}
