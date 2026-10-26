#pragma once
#include "common.h"
#include "util.hpp"
#include "writer_files/writer_util.hpp"

using namespace std;
using namespace COMMON_DATA;
namespace CDT = COMMON_DATA;
namespace RWU = ROUTE_WRITER_UTIL;

namespace READ_DATA {
    char *buf;
    long buf_id, buf_size;
    bool is_digit(char &ch) {
        return ch >= '0' && ch <= '9';
    }
    void get_long(long &val) {
        val = 0;
        while(buf[buf_id] < '0' || buf[buf_id] > '9') {buf_id++;}
        while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
            val = val * 10 + buf[buf_id++] - '0';
        }
    }
    void get_gcd_dt(gcd_dt &val) {
        val = 0;
        while(buf[buf_id] < '0' || buf[buf_id] > '9') {buf_id++;}
        while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
            val = val * 10 + buf[buf_id++] - '0';
        }
    }
    void get_string(string &res) {
        while(buf_id < buf_size && (buf[buf_id] == ' ' || buf[buf_id] == '\n' || buf[buf_id] == ',')) {
            buf_id++;
        }
        res = "";
        while(buf_id < buf_size && (buf[buf_id] != ' ' && buf[buf_id] != '\n' && buf[buf_id] != ',')) {
            res += buf[buf_id];
            buf_id++;
        }
    }
    void get_float(gcb_dt &val) {
        val = 0.0;
        while(buf[buf_id] != '-' && (buf[buf_id] < '0' || buf[buf_id] > '9')) {buf_id++;}
        bool neg = false;
        if(buf[buf_id] == '-') {neg = true; buf_id++;}
        while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
            val = val * 10.0 + buf[buf_id++] - '0';
        }
        if(buf[buf_id] == '.') {
            float mul = 0.1; buf_id++;
            while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
                val = val + float(buf[buf_id++] - '0') * mul;
                mul *= 0.1;
            }
        }
        if(buf[buf_id] == 'e') {
            buf_id++;long mul = 1.0;bool mul_neg = false;int pw = 0;
            if(buf[buf_id] == '-') {mul_neg = true; buf_id++;}
            if(buf[buf_id] == '+') {buf_id++;}
            while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
                pw = pw * 10 + buf[buf_id] - '0'; buf_id++;
            }
            for(int i = 0; i < pw; i++) {val = (mul_neg) ? (val * (1.0)/(10.0)) : (val * 10.0);}
        }
        val = (neg) ? -val : val;

    }
    void get_double(gcb_dt &val) { 
        val = 0.0;
        while(buf[buf_id] != '-' && (buf[buf_id] < '0' || buf[buf_id] > '9')) {buf_id++;}
        bool neg = false;
        if(buf[buf_id] == '-') {neg = true; buf_id++;}
        while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
            val = val * 10.0 + buf[buf_id++] - '0';
        }
        if(buf[buf_id] == '.') {
            gcb_dt mul = 0.1; buf_id++;
            while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
                val = val + gcb_dt(buf[buf_id++] - '0') * mul;
                mul *= 0.1;
            }
        }
        if(buf[buf_id] == 'e') {
            buf_id++;long mul = 1.0;bool mul_neg = false;int pw = 0;
            if(buf[buf_id] == '-') {mul_neg = true; buf_id++;}
            if(buf[buf_id] == '+') {buf_id++;}
            while(buf[buf_id] >= '0' && buf[buf_id] <= '9') {
                pw = pw * 10 + buf[buf_id] - '0'; buf_id++;
            }
            for(int i = 0; i < pw; i++) {val = (mul_neg) ? (val * (1.0)/(10.0)) : (val * 10.0);}
        }
        val = (neg) ? -val : val;
    }
    void read_net_metadata() {
        get_gcd_dt(Z);
        get_gcd_dt(X);
        get_gcd_dt(Y);
        CSIZE = X * Y * Z;
        Z1 = Z; Z2 = CEIL_DIV(Z1, 2); Z1Z2 = Z1 * Z2;
        int XSZ = X, YSZ = Y, ZSZ = Z, CCSIZE = CSIZE, lpat_pts_per_block = NUM_GPU_THREADS / Z1Z2;
        SPDLOG_ASSERT(Z1Z2 < NUM_GPU_THREADS, "");
        spdlog::info("X {} Y {} Z {} CSIZE {}", XSZ, YSZ, ZSZ, CCSIZE);
        CUDAMEMCPY_SYMBOL(DX, &XSZ, 1, int);
        // cudaMemcpyToSymbol(DY, &YSZ, sizeof(int));
        CUDAMEMCPY_SYMBOL(DY, &YSZ, 1, int);
        CUDAMEMCPY_SYMBOL(DZ, &ZSZ, 1, int);
        CUDAMEMCPY_SYMBOL(DZ1, &Z1, 1, int);
        CUDAMEMCPY_SYMBOL(DZ2, &Z2, 1, int);
        CUDAMEMCPY_SYMBOL(DZ1Z2, &Z1Z2, 1, int);
        CUDAMEMCPY_SYMBOL(D_LPAT_PTS_PER_BLK, &lpat_pts_per_block, 1, int);
        CUDAMEMCPY_SYMBOL(DCSIZE, &CCSIZE, 1, int);
        get_double(ULEN_COST);
        get_double(UVIA_COST);
        short_cost_data = (gcb_dt *)malloc(sizeof(gcb_dt) * Z);
        for(int i = 0; i < Z; i++) {
            short_cost_data[i] = 0;
            gcb_dt val=0.0;
            get_double(val);
            short_cost_data[i] = val;
        }
        gcb_dt tmp_val;
        if(NEW_FMT) {get_double(tmp_val); get_double(tmp_val); get_double(tmp_val); get_double(tmp_val);}
    }
    void read_hlen_data() {
        hlen_data = (gcd_dt *)malloc(sizeof(gcd_dt) * (X - 1));
        for(int i = 0; i < X - 1; i++) {
            long val = 0;
            get_long(val);
            hlen_data[i] = val;
        }
    }
    void read_vlen_data() {
        vlen_data = (gcd_dt *)malloc(sizeof(gcd_dt) * (Y - 1));
        for(int i = 0; i < Y - 1; i++) {
            long val = 0;
            get_long(val);
            vlen_data[i] = val;
        }
    }
    void read_layer_data(int layer_id) {
        int cube_id = layer_id * X * Y;
        for(int i = 0; i < Y; i++) {
            for(int j = 0; j < X; j++) {
                genc_dt idx = UTIL::get_enc<CUBE_ENC>(layer_id, j, i);
                gcb_dt val = 0.0; get_float(val); cube_data[idx] = gcb_dt(val);
            }
        }
    }
    void read_cube_data() {
        layer_dirs = (bool *)malloc(sizeof(bool) * Z); min_cost_len = (gcb_dt *)malloc(sizeof(gcb_dt) * Z);
        cube_data = (gcb_dt *)malloc(sizeof(gcb_dt) * X * Y * Z);
        layer_data.resize(Z);
        for(int i = 0; i < Z; i++) {
            get_string(layer_data[i].layer_name);
            long val = 0;
            get_long(val);
            layer_data[i].dir = val; layer_dirs[i] = (val & 1); val = 0;
            get_long(val);
            min_cost_len[i] = val; layer_data[i].min_cost_len = val;
            read_layer_data(i);
        }
        #if EVAL_OFCOST == 1
        CUDA_MALLOC_1D(d_eval_via_cnt, CSIZE, int); CUDA_MEMSET(d_eval_via_cnt, 0, CSIZE, int);
        CUDA_MALLOC_1D(d_eval_wire_cnt, CSIZE, int); CUDA_MEMSET(d_eval_wire_cnt, 0, CSIZE, int);
        #if EVAL_CONN_CHECK == 1
        CUDA_MALLOC_1D(d_eval_vis_data, CSIZE, int);
        #endif
        #endif
    }
    void read_data(string fpath) {
        FILE *f = fopen(fpath.c_str(), "rb");
        if(f == NULL) {
            spdlog::critical("error in reading file {}", fpath);
            exit(0);
        }
        fseek(f, 0, SEEK_END);
        buf_size = ftell(f);
        spdlog::info("reading file {}", fpath);
        spdlog::info("read {} bytes", buf_size);
        fseek(f, 0, SEEK_SET);
        buf = (char *)malloc(1LL * sizeof(char) * (buf_size + 1)); buf[buf_size] = 0;
        fread(buf, buf_size, 1, f); fclose(f); buf_id = 0; buf[buf_size] = 0;
    }
    void read_cap_file(string fpath) {
        TIME_();
        spdlog::info("reading cap file");
        read_data(fpath); read_net_metadata(); read_hlen_data(); read_vlen_data(); read_cube_data();
        buf_id = 0;
        // for(int i= 0; i < Z; i++) {
        //     ldata_t cldata = layer_data[i];
        //     spdlog::info("layer_name {} layer_dir {} layer_min_cost_len {}", cldata.layer_name, cldata.dir, cldata.min_cost_len);
        // }
        // free(buf);        
        spdlog::info("read cap file");
    }
    void get_num_nets() {
        buf_id = 0; NUM_NETS = 0;
        int n_bckts = 0;
        while(buf_id < buf_size) {
            if(buf[buf_id] == '(') {
                n_bckts++;
            } else if(buf[buf_id] == ')') {
                n_bckts--;
                if(n_bckts == 0) {
                    NUM_NETS++;
                }
            }
            buf_id++;
        }
        buf_id = 0;
    }
    
    void cprnt_net_data(int net_id) {
        spdlog::info("{}", nets[net_id].net_name);
        for(pin_t cpin : nets[net_id].pins) {
            spdlog::info("[zxy] [{} {} {}]", cpin.pcd.z, cpin.pcd.x, cpin.pcd.y);
        }
    }

    void get_num_nets_temp() {
        buf_id = 0; NUM_NETS = 0;
        int n_bckts = 0;
        net_id_vec.clear();
        while(buf_id < buf_size) {
            if(buf[buf_id] == '(') {
                if(n_bckts == 0) {
                    net_id_t cnet_id(buf_id, -1);
                    net_id_vec.push_back(cnet_id);
                }
                n_bckts++;
            } else if(buf[buf_id] == ')') {
                n_bckts--;
                if(n_bckts == 0) {
                    net_id_vec.back().ed_id = buf_id;
                    NUM_NETS++;
                }
            }
            buf_id++;
        }
        buf_id = 0;
    }
    void read_comma(long &bid, long &ref_id) {
        while(bid < ref_id && buf[bid] != ',') {
            bid++;
        }
        bid++;
    }
    void read_obrckt(long &bid, long &ref_id) {
        while(bid < ref_id && buf[bid] != '(') {
            bid++;
        }
        bid++;
    }
    int read_cd_till_comma(long &bid) {
        string snum = "";
        while(bid < buf_size && buf[bid] != ',') {
            if(is_digit(buf[bid])) {
                snum += buf[bid];
            }
            bid++;
        }
        bid++;
        // spdlog::info("{}", snum);
        return stoi(snum);
    }
    int read_cd_till_cbrck(long &bid) {
        string snum = "";
        while(bid < buf_size && buf[bid] != ')') {
            if(is_digit(buf[bid])) {
                snum += buf[bid];
            }
            bid++;
        }
        bid++;
        // spdlog::info("{}", snum);
        return stoi(snum);
    }
    void read_nt(int &net_id) {
        long st_id, ed_id;
        if(net_id == 0) {
            st_id = 0; ed_id = net_id_vec[net_id].st_id;
        } else {
            st_id = net_id_vec[net_id - 1].ed_id + 1; ed_id = net_id_vec[net_id].st_id;
        }
        string cnet_name(buf + st_id, buf + ed_id);
        nets[net_id].net_name = move(cnet_name);
        nets[net_id].pins.clear();
        ed_id = net_id_vec[net_id].ed_id;
        while(true) {
            read_comma(st_id, ed_id);
            read_comma(st_id, ed_id);
            if(st_id >= ed_id)
                break;
            pin_t cpin;
            cpin.pcd.z = read_cd_till_comma(st_id);
            cpin.pcd.x = read_cd_till_comma(st_id);
            cpin.pcd.y = read_cd_till_cbrck(st_id);
            nets[net_id].pins.push_back(cpin);
        }
    }
    void read_ispd24_net(int &net_id) {
        long st_id, ed_id;
        if(net_id == 0) {
            st_id = 0; ed_id = net_id_vec[net_id].st_id;
        } else {
            st_id = net_id_vec[net_id - 1].ed_id + 1; ed_id = net_id_vec[net_id].st_id;
        }
        string cnet_name(buf + st_id, buf + ed_id);
        nets[net_id].net_name = move(cnet_name);
        nets[net_id].pins.clear();
        ed_id = net_id_vec[net_id].ed_id;
        while(true) {
            read_obrckt(st_id, ed_id);
            if(st_id + 5 >= ed_id)
                break;
            vector<pin_t> pin_vec;
            pin_vec.clear();
            while(true) {
                pin_t cpin;
                cpin.pcd.z = read_cd_till_comma(st_id);
                cpin.pcd.x = read_cd_till_comma(st_id);
                cpin.pcd.y = read_cd_till_cbrck(st_id);
                pin_vec.push_back(cpin);
                SPDLOG_ASSERT(buf[st_id - 1] == ')', "{}", buf[st_id]);
                if(buf[st_id] == ']') break;
                st_id++;
            }
            SPDLOG_ASSERT(pin_vec.size() > 0, "{}", pin_vec.size());
            nets[net_id].pins.push_back(pin_vec[0]);
        }
    }
    void read_ispd24_net_single(int &net_id) {
        long st_id, ed_id;
        if(net_id == 0) {
            st_id = 0; ed_id = net_id_vec[net_id].st_id;
        } else {
            st_id = net_id_vec[net_id - 1].ed_id + 1; ed_id = net_id_vec[net_id].st_id;
        }
        string cnet_name(buf + st_id, buf + ed_id);
        nets[net_id].net_name = move(cnet_name);
        nets[net_id].pins.clear();
        ed_id = net_id_vec[net_id].ed_id;
        while(true) {
            read_obrckt(st_id, ed_id);
            if(st_id + 5 >= ed_id)
                break;
            pin_t cpin;
            cpin.pcd.z = read_cd_till_comma(st_id);
            cpin.pcd.x = read_cd_till_comma(st_id);
            cpin.pcd.y = read_cd_till_cbrck(st_id);
            SPDLOG_ASSERT(buf[st_id] == ']', "{}", buf[st_id]);
            nets[net_id].pins.push_back(cpin);
        }
    }
    void write_ispd24_nets() {
        unique_ptr<ofstream> GROUT = make_unique<ofstream>(ISPD24_NEW_NF_PATH);
        for(net_t &cnet : nets) {
            string res_str = fmt::format("{}(\n", cnet.net_name);
            for(pin_t &cpin : cnet.pins) {
                res_str += fmt::format("[({}, {}, {})]\n", cpin.pcd.z, cpin.pcd.x, cpin.pcd.y);
            }
            res_str += fmt::format(")");
            (*GROUT) << res_str;
        }
    }
    void remove_single_deg_nets() {
        int fr_id = NUM_NETS - 1;
        NUM_NODES_LEN = 0; ODEG_NUM_NETS = 0;
        for(int i = 0; i < nets.size(); i++) {
            if(nets[i].pins.size() == 1) {
                ODEG_NUM_NETS++;
            } else {
                NUM_NODES_LEN += nets[i].pins.size();
            }
        }
        vector<net_t> tmp_nets;
        one_deg_nets.clear(); one_deg_nets.resize(ODEG_NUM_NETS);
        NUM_NETS -= ODEG_NUM_NETS;
        tmp_nets.resize(NUM_NETS);
        int oid = 0, nid = 0;
        for(int i = 0; i < nets.size(); i++) {
            if(nets[i].pins.size() == 1) {
                one_deg_nets[oid++] = nets[i];
            } else {
                tmp_nets[nid++] = nets[i];
            }
        }
        nets = move(tmp_nets);
        SPDLOG_ASSERT(nid == NUM_NETS, "{} {}", nid, NUM_NETS);
        SPDLOG_ASSERT(oid == ODEG_NUM_NETS, "{} {}", oid, ODEG_NUM_NETS);
        // for(int i = 0; i < fr_id; i++) {
        //     while(fr_id > i && nets[fr_id].pins.size() == 1) {
        //         fr_id--;
        //     }
        //     if(nets[i].pins.size() == 1) {
        //         swap(nets[i], nets[fr_id]);
        //     } else {
        //         NUM_NODES_LEN += nets[i].pins.size();
        //     }
        // }
        // ODEG_NUM_NETS = NUM_NETS - fr_id - 1;
        // spdlog::info("one_deg_nets {}", ODEG_NUM_NETS);
        
        // for(int i = 0; i < ODEG_NUM_NETS; i++) {
        //     one_deg_nets[i] = nets.back();
        //     nets.pop_back();
        // }
        // NUM_NETS = nets.size();
        // for(int i = 0; i < NUM_NETS; i++) {
        //     SPDLOG_ASSERT(nets[i].pins.size() > 1, "unacceptable");
        // }
    }
    void remove_redundant_pins(net_t &cnet) {
        map<in_t, pin_t> mp;
        mp.clear();
        for(pin_t cpin : cnet.pins) {
            genc_dt enc = UTIL::get_enc<CUBE_ENC>(cpin.pcd.z, cpin.pcd.x, cpin.pcd.y);
            if(mp.count(enc) == 0) {
                mp[enc] = cpin;
            } else {
                pin_t opin = mp[enc];
                if(opin.pslack > cpin.pslack) {
                    mp[enc] = cpin;
                }
            }
        }
        cnet.pins.clear(); cnet.pins.resize(mp.size());
        int id = 0;
        for(auto  mp_nd : mp) {
            cnet.pins[id++] = mp_nd.second;
        }
    }
    void route_one_deg_nets() {
        vector<net_t> temp_nets;
        temp_nets.clear();
        for(net_t cnet : nets) {            
            if(cnet.pins.size() == 1) {
                string net_gr = "";
                pin_t cpin = cnet.pins[0];
                #if FWRITE_ASSERT == 1
                SPDLOG_ASSERT(cpin.pcd.z == 0, "z cd {}", cpin.pcd.z);
                #endif
                RWU::write_net_name(net_gr, cnet.net_name);
                // UTIL::write_net_name(cnet.net_name);
                if(RSEG_FMT) {
                    RWU::add_vias_cpu(net_gr, RWU::get_gran_cd(cpin.pcd.x), RWU::get_gran_cd(cpin.pcd.y), cpin.pcd.z, cpin.pcd.z+1);
                } else {
                    RWU::add_vias_cpu(net_gr, (cpin.pcd.x), (cpin.pcd.y), cpin.pcd.z, cpin.pcd.z+1);
                }                
                // UTIL::add_vias_cpu(UTIL::get_gran_cd(cpin.pcd.x), UTIL::get_gran_cd(cpin.pcd.y), cpin.pcd.z, cpin.pcd.z+1);
                // cube_data[UTIL::get_enc<CUBE_ENC>(cpin.pcd.z + 1, cpin.pcd.x, cpin.pcd.y)] -= WIRE_WT_PDT;
                cube_data[UTIL::get_enc<CUBE_ENC>(cpin.pcd.z, cpin.pcd.x, cpin.pcd.y)] -= WIRE_WT_PDT;
                RWU::write_closed_bracket(net_gr);
                gr_file_logger->info(net_gr);
            } else {
                temp_nets.push_back(cnet);
            }
        }
        nets.clear();
        NUM_NODES_LEN = 0;
        for(net_t cnet : temp_nets) {
            nets.push_back(cnet);
            NUM_NODES_LEN += cnet.pins.size();
        }
        NUM_NETS = nets.size();
        spdlog::info("{} NUM_NETS have degree > 1", NUM_NETS);
    }
    void remove_large_nets() {
        vector<net_t> tmp_nets;
        int id = 0;
        for(net_t cnet : nets) {
            if(cnet.pins.size() <= MAX_NUM_PINS) {
                remove_redundant_pins(cnet);
                tmp_nets.push_back(cnet);
            } else {
                spdlog::info("ignoring_net degree {} name {} id {}", cnet.pins.size(), cnet.net_name, id);
            }
            id++;
        }
        nets.clear();
        for(net_t cnet : tmp_nets) {
            nets.push_back(cnet);
        }
        NUM_NETS = nets.size();
    }
void add_net_gpu_data() {
        TIME_();
        //multi_thread this later
        vector<int> deg_data; deg_data.clear(); deg_data.resize(NUM_NETS);
        for(int i = 0; i < NUM_NETS; i++) {
            deg_data[i] = nets[i].pins.size();
        }
        CUDA_MALLOC_1D(d_net_deg_data, NUM_NETS, int);
        CUDAMEMCPYHD_1D(d_net_deg_data, deg_data.data(), NUM_NETS, int);
        CUDA_MALLOC_1D(d_npart_ids, NUM_NETS, int);
        CUDA_MALLOC_1D(d_net_pins_len, NUM_NETS, long);
        CUDA_MALLOC_1D(d_net_pins_pfx, NUM_NETS, long);
        int num_blks = CEIL_DIV(NUM_NETS, NUM_GPU_THREADS);
        UTIL::copy_deg_net_len<<<num_blks, NUM_GPU_THREADS>>>(NUM_NETS);
        cudaDeviceSynchronize();
        thrust::device_ptr<long> d_nets_len_ptr = thrust::device_pointer_cast(d_net_pins_len);
        thrust::device_ptr<long> d_nets_len_optr = thrust::device_pointer_cast(d_net_pins_pfx);
        thrust::exclusive_scan(d_nets_len_ptr, d_nets_len_ptr + NUM_NETS, d_nets_len_optr);
        long max_ofst, max_len;
        CUDAMEMCPYDH_1D(&max_ofst, d_net_pins_pfx + NUM_NETS - 1, 1, long);
        CUDAMEMCPYDH_1D(&max_len, d_net_pins_len + NUM_NETS - 1, 1, long);
        long max_pins = max_len + max_ofst;
        CUDA_MALLOC_1D(d_net_pins, max_pins, int);
        net_pins.clear(); net_pins.resize(max_pins);
        long fid = 0;
        for(int i = 0; i < NUM_NETS; i++) {
            net_t cnet = nets[i];
            for(pin_t cpin : cnet.pins) {
                gcd_dt x = cpin.pcd.x, y = cpin.pcd.y;
                gvia_dt z = cpin.pcd.z;
                genc_dt enc = UTIL::get_enc<CUBE_ENC>(z, x, y);
                net_pins[fid++] = enc;
            }
        }
        #if CPU_ASSERT == 1
        SPDLOG_ASSERT(fid == max_pins, "{} {}", fid, max_pins);
        #endif
        CUDAMEMCPYHD_1D(d_net_pins, net_pins.data(), max_pins, genc_dt);
    }
    void get_net_bbx_data(int &net_id) {
        bbx_t cbx_data; cbx_data.minx = MAXX; cbx_data.miny = MAXY; cbx_data.minz = MAXZ; cbx_data.maxx = -1; cbx_data.maxy = -1; cbx_data.maxz = -1;
        for(pin_t &cpin_data : nets[net_id].pins) {
            cbx_data.minx = min(cbx_data.minx, cpin_data.pcd.x); cbx_data.maxx = max(cbx_data.maxx, cpin_data.pcd.x);
            cbx_data.miny = min(cbx_data.miny, cpin_data.pcd.y); cbx_data.maxy = max(cbx_data.maxy, cpin_data.pcd.y);
            cbx_data.minz = min(cbx_data.minz, cpin_data.pcd.z); cbx_data.maxz = max(cbx_data.maxz, cpin_data.pcd.z);
        }
        nets[net_id].bx_data = cbx_data;
    }
    void read_net_file(string fpath) {
        TIME_();
        spdlog::info("reading net file");
        read_data(fpath);
        get_num_nets_temp();
        spdlog::info("NUM_NETS : {} {}", NUM_NETS, net_id_vec.size());
        nets.resize(NUM_NETS);
        int net_id = 0, num_bckts = 0;
        #pragma omp parallel for num_threads(7)
        for(int i = 0; i < NUM_NETS; i++) {
            if(NEW_FMT) {
                read_nt(i);
            } else {
                read_ispd24_net_single(i);
            }
            remove_redundant_pins(nets[i]);
            get_net_bbx_data(i);
        }
        free(buf);
        remove_single_deg_nets();
        CUDA_MALLOC_1D(d_cube_data, CSIZE, gcb_dt);
        CUDA_MALLOC_1D(d_pre_demand_data, CSIZE, int);
        CUDAMEMCPYHD_1D(d_cube_data, cube_data, CSIZE, gcb_dt);
        if(BIT_PMODEL) {
            CUDA_MALLOC_1D(d_cube_data_bit, CSIZE, gcb_dt);
            CUDA_MALLOC_1D(d_cost_data_bit, CSIZE, gcb_dt);
            CUDA_MALLOC_1D(d_bit_vis_data, CSIZE, bool);
            CUDAMEMCPYHD_1D(d_cube_data_bit, cube_data, CSIZE, gcb_dt);
        }
        add_net_gpu_data();
        spdlog::info("read net file");
    }
    void read_write_ispd24_nets(string fpath) {
        TIME_();
        spdlog::info("reading net file");
        read_data(fpath);
        get_num_nets_temp();
        spdlog::info("NUM_NETS : {} {}", NUM_NETS, net_id_vec.size());
        nets.resize(NUM_NETS);
        int net_id = 0, num_bckts = 0;
        #pragma omp parallel for num_threads(7)
        for(int i = 0; i < NUM_NETS; i++) {
            read_ispd24_net(i);
        }
        free(buf);
        SPDLOG_ASSERT(ISPD24_NEW_NF_PATH != "", "{}", ISPD24_NEW_NF_PATH);
        write_ispd24_nets();
    }
}
