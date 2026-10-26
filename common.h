#pragma once
#include <queue>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <cassert>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/box.hpp>
#include "cuda.h"
#include "cuda_runtime.h"
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <semaphore>
#include <thread>
#include <algorithm>
#include <spdlog/fmt/fmt.h>
// #include "cuda_api.h"

using namespace std;
typedef long in_t;
typedef long zp_dt;
typedef double fr_t;
typedef float gcb_dt;
typedef int genc_dt;
typedef int gpfx_dt;
typedef int gvia_dt;
typedef int gcd_dt;
typedef int djdir_dt;
typedef int pth_dt;
typedef double eval_dt ;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

//EVAL PARAMS
#define EVAL_OFCOST 0
#define EVAL_CONN_CHECK 0
#define EVAL_CONN_VIS_BITS 26
#define EVAL_VIA_OFCOST 0

//WRITE_PARAMS
#define WRITE_ROUTE 1

//COST TYPE
#define LIVE_COST 1

//DB PARAMS
#define WT_FLOAT 1
#define ADD_PART_DATA_DB 0
#define ADD_NET_METADATA_DB 0
#define QUERY_EXEC_THRESHOLD 500000
#define EDGE_QUERY_EXEC_THRESHOLD 1000000
#define ADD_LPAT_PERF_DATA_DB 0
#define ADD_RANDOM_LPAT_PERF_DB 0
#define ADD_EDGE_DATA 0

#define FWRITE_ASSERT 0

#define ZPAT_ENCBITS 32
#define ZPAT_ENCMASK ((1L << ZPAT_ENCBITS) - 1L)
#define ZPAT_PATBITS 4
#define ZPAT_PATMASK ((1 << ZPAT_PATBITS) - 1)
#define DJKST_QUEUE 1
#define RANDOM_GEN 0
#define VIA_WT_PDT 1.0
#define WIRE_WT_PDT 1
#define CD_SYSTEM 1
#define VIA_MINOR 0
#define MAXX 100000
#define MAXY 100000
#define MAXZ 100
#define CEIL_DIV(a, b) (((a) + (b) - 1) / (b))
#define ZDIR 0
#define XDIR 1
#define YDIR 2
#define NUM_GPU_THREADS 256
#define MAX_NUM_PINS 20000
#define MAX_EDGE_DATA_PINS 2000
#define MAX_OFGPU_COST gcb_dt(1.0 * 1e195)
#define GPU_INF_COST gcb_dt(1.0 * 1e150)
#define NET_PART_CHECK 0
#define OFCOST_FACTOR 1
#define MAX_BATCH_SIZE 2000000
#define DOUBLE_ACC_VAL (1.0 * 0.0000000001)
#define FLOAT_ACC_VAL (1.0 * 0.001)
#define NET_ID_ORDER_REV 1

#define WIRE_ENC 0
#define VIA_ENC 1
#define WV_ENC WIRE_ENC
#define CUBE_ENC WIRE_ENC
#define INC_OF_FORMULA 1
#define HUGE_NET_DEG 8000
// #define 

#define NUM_PIN_CDS 3
#define NUM_BBX_CDS 6
#define CPU_ASSERT 0
#define ZCD 0
#define XCD 1
#define YCD 2
#define ZMIN_CD 0
#define ZMAX_CD 1
#define XMIN_CD 2
#define XMAX_CD 3
#define YMIN_CD 4
#define YMAX_CD 5

typedef bg::model::point<gcd_dt, 2, bg::cs::cartesian> point_t;
typedef bg::model::box<point_t> box_t;
typedef std::pair<box_t, int> rtree_t;
typedef bg::index::rtree<box_t, bg::index::rstar<16>> bxtree_t;


#define CUDA_MALLOC_1D(ptr, count, type) \
    cudaMalloc((void**)&(ptr), (count) * sizeof(type))

#define CUDAMEMCPYHD_1D(ptr_dst, ptr_src, count, type) \
    cudaMemcpy(ptr_dst, ptr_src, (count) * sizeof(type), cudaMemcpyHostToDevice)

#define CUDAMEMCPYDD_1D(ptr_dst, ptr_src, count, type) \
    cudaMemcpy(ptr_dst, ptr_src, (count) * sizeof(type), cudaMemcpyDeviceToDevice)

#define CUDAMEMCPYDH_1D(ptr_dst, ptr_src, count, type) \
    cudaMemcpy(ptr_dst, ptr_src, (count) * sizeof(type), cudaMemcpyDeviceToHost)

#define CUDAMEMCPY_SYMBOL(ptr_dst, ptr_src, count, type) \
    cudaMemcpyToSymbol(ptr_dst, ptr_src, (count) * (sizeof(type)))

#define CUDA_MEMSET(ptr, value, count, type) \
    cudaMemset(ptr, value, ((count) * sizeof(type)))

#define SPDLOG_ASSERT(cond, ...) \
    do { \
        if (!(cond)) { \
            spdlog::critical(__VA_ARGS__); \
            assert(cond); \
        } \
    } while (0)

class MeasureTime { // measures run time of a function call using high resolution CPU clock
    private:
        const std::string _name;
        std::chrono::nanoseconds* _rt;
        std::chrono::high_resolution_clock::time_point _begin;
    public:
        MeasureTime(const std::string& name, std::chrono::nanoseconds* rt = nullptr) : _name(name), _rt(rt)
        {
        _begin = std::chrono::high_resolution_clock::now();
        }
        ~MeasureTime()
        {
        auto difft = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - _begin);
        if (_rt) {
            (*_rt) += difft;
        } else {
            spdlog::info("{} runtime : {}(s)", _name, difft.count()/1.e9);
        }
        }
};
  
  
  
#define TIME_A(X) MeasureTime __FUNC__##t(__PRETTY_FUNCTION__, X)
#define TIME_()  MeasureTime __FUNC__##t(__PRETTY_FUNCTION__)
#define TIME_S(str)  MeasureTime __FUNC__##t(str)

namespace COMMON_DATA {
    class bit_tree {
        public:
        vector<gcd_dt> bit_arr;
        int sz;
        bit_tree() {
            bit_arr.clear();
        }
        bit_tree(int csz) {
            sz = csz + 1;
            bit_arr.clear(); bit_arr.resize(sz + 1, 0);
        }
        bit_tree(vector<gcd_dt> &in_vec) {
            bit_arr.clear(); bit_arr.resize(in_vec.size() + 1, 0);
            sz = bit_arr.size();
            for(int i = 0; i < in_vec.size(); i++) {
                update_bit(i, in_vec[i]);
            }
        }
        ~bit_tree() {
            bit_arr.clear();
        }
        void reset() {
            bit_arr.clear(); bit_arr.resize(sz + 1, 0);
        }
        void update_bit(int id, gcd_dt val) {
            id++;
            while(id < sz) {
                bit_arr[id] += val;
                id = id + (id & (-id));
            }
        }
        gcd_dt query_bit(int id) {
            gcb_dt res = 0;
            id++;
            while(id > 0) {
                res += bit_arr[id];
                id = id - (id & (-id));
            }
            return res;
        }
        gcb_dt range_query(int l, int r) {
            if(l > r) {
                return gcb_dt(0);
            }
            l++; r++;
            return query_bit(r) - query_bit(l - 1);
        }
    };
    class ru_bit_tree {
        public:
        unique_ptr<bit_tree> bit1, bit2;
        int sz;
        ru_bit_tree() {
            bit1 = make_unique<bit_tree>();
            bit2 = make_unique<bit_tree>();
        }
        ru_bit_tree(int sz) {
            bit1 = make_unique<bit_tree>(sz);
            bit2 = make_unique<bit_tree>(sz);
        }
        ru_bit_tree(vector<gcd_dt> &in_vec) {
            sz = in_vec.size();
            bit1 = make_unique<bit_tree>(sz);
            bit2 = make_unique<bit_tree>(sz);
            for(int i = 0; i < sz; i++) {
                range_update(i, i, in_vec[i]);
            }
        }
        void reset() {
            bit1->reset(); bit2->reset();
        }
        void range_update(int l, int r, gcd_dt val) {
            if(l > r || val == 0) {
                return;
            }
            bit1->update_bit(l, val);
            bit1->update_bit(r + 1, -val);
            bit2->update_bit(l, (l - 1) * val);
            bit2->update_bit(r + 1, -r * val);
        }
        gcd_dt query_bit(int id) {
            return bit1->query_bit(id) * id - bit2->query_bit(id);
        }
        gcd_dt range_query(int l, int r) {
            if (l > r) {
                return gcd_dt(0);
            }
            return query_bit(r) - query_bit(l - 1);
        }
    };
    class minterval_t {
        public:
        genc_dt enc1, enc2;
        minterval_t(genc_dt cenc1, genc_dt cenc2) {
            enc1 = cenc1; enc2 = cenc2;
        }
    };
    struct pt_prs {
        genc_dt enc1, enc2;
    };
    struct ldata_t {
        string layer_name;
        bool dir;
        gcb_dt min_cost_len;
    };
    struct bbx_t {
        gcd_dt minx, miny, maxx, maxy;
        gvia_dt minz, maxz;
    };
    struct coord_t {
        gcd_dt x, y, z;
    };
    struct pin_t {
        gcb_dt pslack;
        coord_t pcd;
        string pname;
    };
    struct net_t {
        bbx_t bx_data;
        string net_name;
        vector<pin_t> pins;
    };
    class net_id_t {
        public:
            long st_id, ed_id;
            net_id_t(long st_id_0, long ed_id_0) {
                st_id = st_id_0; ed_id = ed_id_0;
            }
    };
    unique_ptr<ru_bit_tree> RU_BIT;
    pt_prs *rd_enc;
    gcd_dt X, Y, Z, CSIZE;
    gvia_dt Z1, Z2, Z1Z2;
    __constant__ int DX, DY, DZ, DCSIZE, DZ1, DZ2, DZ1Z2, D_LPAT_PTS_PER_BLK;
    gcb_dt ULEN_COST, UVIA_COST, *short_cost_data, *cube_data, *min_cost_len;
    int NUM_NETS, ODEG_NUM_NETS = 0;
    gcd_dt *hlen_data, *vlen_data;
    int *npart_ofst, *npart_len, *npart_ids;
    __managed__ gcb_dt *d_wcost_cube, *d_vcost_cube;
    __managed__ gcb_dt *d_cube_data, *d_cube_data_bit, *d_cost_data_bit;
    __managed__ bool BIT_PMODEL = false;
    __managed__ bool *d_bit_vis_data;
    gcb_dt *wcost_cube, *vcost_cube, *hd_wcost_cube, *hd_lpat_res;
    bool *layer_dirs, NEW_FMT, RSEG_FMT = false;
    long TOT_LPAT_PRS;
    vector<ldata_t> layer_data;
    vector<genc_dt> net_pins; 
    bxtree_t net_bbx;
    bool RBT_PART = false, ORD_PART = false, BATCH_PART = false;
    vector<net_t> nets, one_deg_nets;
    int XDIM_NETS = 50;
    binary_semaphore net_gen_bsem(0), write_route_bsem(0);
    int PART_ID, PART_OFST;
    bool GEN_STATE;
    int MAX_EXP = 17;
    gcb_dt MAX_OFCOST = 1e75;
    gcb_dt diff_tolerance = 1e-9;
    gcd_dt LPAT_NUM_PTS = 1e6;
    string NVOL_FILE="vol.txt", CHIP_NAME, EDGE_QUERY;
    int MAX_DJKST_BENDS = 2, NUM_NODES_LEN = 0;
    long EDGE_ID = 0;
    string CAP_FILE, NET_FILE;
    vector<net_id_t> net_id_vec;
    gcd_dt VIA_WT = VIA_WT_PDT, WIRE_WT = WIRE_WT_PDT;
    int NUM_XDIM_PARTS = 4, WR_THREADS = 4;
    bool WRT_MRG = false, OPRD_PPATCH = false;
    FILE *OUT_FILE;
    string OUT_FILE_PATH, ISPD24_NEW_NF_PATH = "";
    bool UPD_ISPD24_NET_FILe = false;
    shared_ptr<spdlog::logger> gr_file_logger;
    __managed__ genc_dt *d_net_pins;
    __managed__ long *d_net_pins_pfx, *d_net_pins_len;
    __managed__ int *d_net_deg_data, *d_pre_demand_data, *d_npart_ids;
    __managed__ int *d_eval_wire_cnt, *d_eval_via_cnt, *d_eval_vis_data;
    //init test data structures here
}
