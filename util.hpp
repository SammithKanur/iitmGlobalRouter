#pragma once
#include "common.h"

using namespace std;
using namespace COMMON_DATA;
namespace CDT = COMMON_DATA;

namespace UTIL {
    template<typename T>
    void get_epfx(T *r_arr, T *arr, int len) {
        T pfx = 0;
        for(int i = 0; i < len; i++) {
            r_arr[i] = pfx;
            pfx += arr[i];
        }
    }
    template<typename T>
    void get_vec_epfx(vector<T> &r_arr, vector<T> &arr, int len) {
        T pfx = 0;
        for(int i = 0; i < len; i++) {
            r_arr[i] = pfx;
            pfx += arr[i];
        }
    }    
    void write_bbx() {
        string res = "";
        for(int i = 0; i < NUM_NETS; i++) {
            net_t cnet = nets[i];
            in_t xd = cnet.bx_data.maxx - cnet.bx_data.minx + 1, yd = cnet.bx_data.maxy - cnet.bx_data.miny + 1, zd = cnet.bx_data.maxz - cnet.bx_data.minz + 1;
            in_t vol = xd*yd*zd;
            res += to_string(xd*10*yd) + "\n";
        }
        FILE *f = fopen(NVOL_FILE.c_str(), "w");
        fprintf(f, "%s", res.c_str());
        fclose(f);
    }
    template<int tp>
    inline genc_dt get_enc(gvia_dt z, gcd_dt x, gcd_dt y) {
        if(tp == WIRE_ENC) {
            return (z&1) ? z * X * Y + x * Y + y : z * X * Y + y * X + x;
        } else if(tp == VIA_ENC) {
            return (z&1) ? x * Z * Y + y * Z + z : y * Z * X + x * Z + z; 
        }
    }
    template<int tp>
    __inline__ __device__ genc_dt get_gpu_enc(gvia_dt z, gcd_dt x, gcd_dt y) {
        if(tp == WIRE_ENC) {
            return ((z&1) ? z * DX * DY + x * DY + y : z * DX * DY + y * DX + x);
        } else if(tp == VIA_ENC) {
            return ((z&1) ? x * DZ * DY + y * DZ + z : y * DZ * DX + x * DZ + z);
        }
    }
    template<int tp>
    inline void get_dec_all(genc_dt val, gvia_dt &z, gcd_dt &x, gcd_dt &y) {
        if(tp == WIRE_ENC) {
            z = val/(X*Y); x=(z&1)?val/Y%X:val%X; y=(z&1)?val%Y:val/X%Y;
        } else if(tp == VIA_ENC) {
            z = val % Z; x = (z&1)?val/(Z*Y):(val/Z)%X; y = (z&1)?(val/Z)%Y:val/(X*Z);
        }
    }
    template<int tp>
    __inline__ __device__ void get_gpu_dec_all(genc_dt &val, gvia_dt &z, gcd_dt &x, gcd_dt &y) {
        if(tp == WIRE_ENC) {
            z = (val)/(DX*DY); x=((z&1)?val/DY%DX:val%DX); y=((z&1)?val%DY:val/DX%DY);
        } else if(tp == VIA_ENC) {
            z = val % DZ; x = (z&1)?val/(DZ*DY):(val/DZ)%DX; y = (z&1)?(val/DZ)%DY:val/(DX*DZ);
        }
    }
    inline void get_dec_z(genc_dt &val, gvia_dt &z) {
        z=val/X/Y;
    }
    inline void get_dec_x(genc_dt &val, gcd_dt &x) {
        x=((val/X/Y)&1)?val/Y%X:val%X;
    }
    inline void get_dec_y(genc_dt &val, gcd_dt &y) {
        y=((val/X/Y)&1)?val%Y:val/X%Y;
    }
    inline void get_bbx_size(bbx_t &bx, int &res) {
        res = (bx.maxx - bx.minx + 1) * (bx.maxy - bx.miny + 1) * 10;
    }
    inline int get_bbx_sq_size(bbx_t &bx) {
        return (bx.maxx - bx.minx + 1) * (bx.maxy - bx.miny + 1);
    }
    void get_num_pts_batch(long &res, int batch_id) {
        int ofst = npart_ofst[batch_id], len = npart_len[batch_id];
        res = 0;
        for(int i = ofst; i < ofst + len; i++) {
            int nid = npart_ids[i];            
            long num_pins = nets[nid].pins.size();
            res += ((num_pins) * (num_pins - 1)) >> 1;
        }
    }
    void get_batches(long &MAX_RESOURCES, long *itms, int &num_items, long *bofst, long *blen, long &num_batches) {
        long dmd = 0, cblen = 0, onum_batches = 0;
        for(int i = 0; i < num_items; i++) {
            long cresource = itms[i];
            SPDLOG_ASSERT(cresource <= MAX_RESOURCES, "single item demand too high cur_demand {} max_resource {}", cresource, MAX_RESOURCES);
            if(dmd > MAX_RESOURCES) {
                onum_batches++; dmd = 0;
            }
            dmd += cresource;
        }
        blen = (long *)malloc(sizeof(long) * onum_batches);
        bofst = (long *)malloc(sizeof(long) * onum_batches);
        onum_batches++; num_batches = 0; cblen = 0;
        for(int i = 0; i < num_items; i++) {
            long cresource = itms[i];
            SPDLOG_ASSERT(cresource <= MAX_RESOURCES, "single item demand too high cur_demand {} max_resource {}", cresource, MAX_RESOURCES);
            if(dmd > MAX_RESOURCES) {
                blen[num_batches++] = cblen; cblen = 0; dmd = 0;
            }
            dmd += cresource;
            cblen++;
        }
        num_batches++;
        SPDLOG_ASSERT(onum_batches == num_batches, "onumbatches not matching num_batches {} {}", onum_batches, num_batches);
        get_epfx<long>(bofst, blen, num_batches);
    }
    __inline__ __device__ gcb_dt get_gpu_of_cost(gcb_dt cap) {
        // if(cap < 0.1) {
        //     return GPU_INF_COST;
        // }
        return (expf(-cap * 0.5));
        // cap *= 0.5;
        // return -(cap * cap * cap * cap * cap);
    }
    __inline__ __device__ gcb_dt get_gpu_of_cost_inc(gcb_dt &cap, float inc_val) {
        return get_gpu_of_cost(cap - inc_val) - get_gpu_of_cost(cap);
    }
    __inline__ gcb_dt get_cpu_of_cost(gcb_dt &cap) {
        return (expf(-cap * 0.5) * OFCOST_FACTOR);
    }
    __global__ void get_gpu_cost_cube(gcb_dt *wcost_cube, gcb_dt *vcost_cube, gcb_dt *cube_data, int mx_items) {
        genc_dt tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= mx_items) {
            return;
        }
        gcb_dt cur_cap = cube_data[tid];
        #if INC_OF_FORMULA == 1
        wcost_cube[tid] = get_gpu_of_cost_inc(cur_cap, 1);
        #else
        cost_cube[tid] = get_gpu_of_cost(cur_cap);
        #endif
        gvia_dt z;
        gcd_dt x, y;
        get_gpu_dec_all<WIRE_ENC>(tid, z, x, y);
        // if(z == 0) return;
        genc_dt vtid = get_gpu_enc<VIA_ENC>(z, x, y);
        if(z & 1) {
            if(y > 0 && y < DY - 1) {
                #if INC_OF_FORMULA == 1
                vcost_cube[vtid - 1] = get_gpu_of_cost_inc(cube_data[tid-1], 0.5);
                vcost_cube[vtid] = get_gpu_of_cost_inc(cube_data[tid], 0.5);
                #else
                vcost_cube[vtid - 1] = get_gpu_of_cost(cube_data[tid]);
                vcost_cube[vtid] = get_gpu_of_cost(cube_data[tid]);
                #endif
            } else if(y == 0) {
                #if INC_OF_FORMULA == 1
                vcost_cube[vtid] = get_gpu_of_cost_inc(cube_data[tid], 1);
                #else
                vcost_cube[vtid] = get_gpu_of_cost(cube_data[tid]);
                #endif
            } else {
                #if INC_OF_FORMULA == 1
                vcost_cube[vtid-1] = get_gpu_of_cost_inc(cube_data[tid-1], 1);
                #else
                vcost_cube[vtid-1] = get_gpu_of_cost(cube_data[tid-1]);
                #endif
            }   
        } else {
            if(x > 0 && x < DX - 1) {
                #if INC_OF_FORMULA == 1
                vcost_cube[vtid - 1] = get_gpu_of_cost_inc(cube_data[tid-1], 0.5);
                vcost_cube[vtid] = get_gpu_of_cost_inc(cube_data[tid], 0.5);
                #else
                vcost_cube[vtid - 1] = get_gpu_of_cost(cube_data[tid]);
                vcost_cube[vtid] = get_gpu_of_cost(cube_data[tid]);
                #endif
            } else if(x == 0) {
                #if INC_OF_FORMULA == 1
                vcost_cube[vtid] = get_gpu_of_cost_inc(cube_data[tid], 1);
                #else
                vcost_cube[vtid] = get_gpu_of_cost(cube_data[tid]);
                #endif
            } else {
                #if INC_OF_FORMULA == 1
                vcost_cube[vtid-1] = get_gpu_of_cost_inc(cube_data[tid-1], 1);
                #else
                vcost_cube[vtid-1] = get_gpu_of_cost(cube_data[tid-1]);
                #endif
            }   
        }
    }    
    void get_gpu_pfx_sum_blks(int &res) {
        res = 0;
        for(int i = 0; i < Z; i++) {
            res += (i & 1) ? X : Y;
        }
    }
    __global__ void d_get_pfx_cost_cube(gcb_dt *cost_cube, gcb_dt *pfx_cost_cube, int nth, int mem_size) {
        extern __shared__ gcb_dt smarr[];
        __shared__ gpfx_dt clyr;
        __shared__ int ctrk, enc_ofst, rem_val, thload;
        __shared__ gcb_dt pfx_ofst;
        if(threadIdx.x >= nth) {
            return;
        }
        if(threadIdx.x == 0) {
            clyr = 0; ctrk = 0; pfx_ofst = 0;
            for(gpfx_dt i = 0; i < DZ && ctrk <= blockIdx.x; i++) {
                ctrk += (i & 1) ? DX : DY;
                clyr++;
            }
            clyr--;
            // assert(clyr >= 0 && clyr < DZ);
            ctrk -= (clyr&1) ? DX:DY; assert(ctrk <= blockIdx.x);
            ctrk = blockIdx.x - ctrk; enc_ofst = DX * DY * clyr + (clyr & 1) ? DX * ctrk : DY * ctrk; rem_val = (clyr & 1) ? DY : DX;
            enc_ofst = DX * DY * clyr + ((clyr & 1) ? DY * ctrk : DX * ctrk);
            thload = CEIL_DIV(mem_size, nth);
            enc_ofst = 0;
        }
        __syncthreads();
        while(rem_val > 0) {
            __syncthreads();
            //load values to smarr
            gvia_dt zcd = clyr;
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                smarr[i] = cost_cube[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //get the prefix sum for each partition parallely
            gcb_dt pfx=smarr[max(0, min(threadIdx.x * thload - 1, mem_size-1))], prev;
            if(threadIdx.x == 0) pfx = 0;
            __syncthreads();
            for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload, mem_size); i++) {
                prev=smarr[i]; smarr[i]=pfx; pfx+=prev;
            }
            __syncthreads();
            //propagate the pfx end value of each partition to other partitions
            if(threadIdx.x == 0) {
                pfx = smarr[thload - 1];
                for(int i = thload + thload - 1; i < mem_size - 1; i += thload) {
                    pfx += smarr[i]; smarr[i] = pfx;
                }
                if(thload - 1 != mem_size - 1) {
                    pfx += smarr[mem_size - 1]; smarr[mem_size - 1] = pfx;
                }                
            }
            __syncthreads();
            //spread the propagated value to each partition
            if(threadIdx.x > 0) {
                pfx = smarr[min(threadIdx.x * thload - 1, mem_size-1)];
                for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload-1, mem_size-1); i++) {
                    smarr[i] += pfx;
                }
            }
            __syncthreads();
            //write back values to dram
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {                
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                pfx_cost_cube[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)] = smarr[i] + pfx_ofst;
            }
            __syncthreads();
            //update enc_ofst and rem_val
            if(threadIdx.x == 0) {
                gcd_dt xcd = ((zcd&1)?ctrk:enc_ofst + min(rem_val, mem_size) - 1), ycd = ((zcd&1)?enc_ofst + min(rem_val, mem_size) - 1:ctrk);
                pfx_ofst += pfx + cost_cube[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)];
                // pfx_ofst += smarr[min(mem_size, rem_val) - 1];
                enc_ofst += mem_size;
                rem_val -= mem_size;
            }
            __syncthreads();
        }
    }
    __global__ void d_get_epfx_cost_cube_inplace(gcb_dt *cost_cube, int nth, int mem_size) {
        extern __shared__ gcb_dt smarr[];
        __shared__ gpfx_dt clyr;
        __shared__ int ctrk, enc_ofst, rem_val, thload;
        __shared__ gcb_dt pfx_ofst, tpfx_ofst;
        if(threadIdx.x >= nth) {
            return;
        }
        if(threadIdx.x == 0) {
            clyr = 0; ctrk = 0; pfx_ofst = 0; tpfx_ofst = 0;
            for(gpfx_dt i = 0; i < DZ && ctrk <= blockIdx.x; i++) {
                ctrk += (i & 1) ? DX : DY;
                clyr++;
            }
            clyr--;
            // assert(clyr >= 0 && clyr < DZ);
            ctrk -= (clyr&1) ? DX:DY; assert(ctrk <= blockIdx.x);
            ctrk = blockIdx.x - ctrk; enc_ofst = DX * DY * clyr + (clyr & 1) ? DX * ctrk : DY * ctrk; rem_val = (clyr & 1) ? DY : DX;
            enc_ofst = DX * DY * clyr + ((clyr & 1) ? DY * ctrk : DX * ctrk);
            thload = CEIL_DIV(mem_size, nth);
            enc_ofst = 0;
        }
        __syncthreads();
        while(rem_val > 0) {
            __syncthreads();
            //load values to smarr
            gvia_dt zcd = clyr;
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                smarr[i] = cost_cube[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //get the prefix sum for each partition parallely
            gcb_dt pfx=smarr[max(0, min(threadIdx.x * thload - 1, mem_size-1))], prev;
            if(threadIdx.x == 0) pfx = 0;
            __syncthreads();
            for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload, mem_size); i++) {
                prev=smarr[i]; smarr[i]=pfx; pfx+=prev;
            }
            __syncthreads();
            //propagate the pfx end value of each partition to other partitions
            if(threadIdx.x == 0) {
                pfx = smarr[thload - 1];
                for(int i = thload + thload - 1; i < mem_size - 1; i += thload) {
                    pfx += smarr[i]; smarr[i] = pfx;
                }
                if(thload - 1 != mem_size - 1) {
                    pfx += smarr[mem_size - 1]; smarr[mem_size - 1] = pfx;
                }
                gcd_dt xcd = ((zcd&1)?ctrk:enc_ofst + min(rem_val, mem_size) - 1), ycd = ((zcd&1)?enc_ofst + min(rem_val, mem_size) - 1:ctrk);
                tpfx_ofst = cost_cube[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)];
            }
            __syncthreads();
            //spread the propagated value to each partition
            if(threadIdx.x > 0) {
                pfx = smarr[min(threadIdx.x * thload - 1, mem_size-1)];
                for(int i = threadIdx.x * thload; i < min(threadIdx.x * thload + thload-1, mem_size-1); i++) {
                    smarr[i] += pfx;
                }
            }
            __syncthreads();
            //write back values to dram
            for(int i = threadIdx.x; i < min(rem_val, mem_size); i += nth) {
                gcd_dt xcd = ((zcd&1)?ctrk:i+enc_ofst), ycd = ((zcd&1)?i+enc_ofst:ctrk);
                cost_cube[UTIL::get_gpu_enc<CUBE_ENC>(zcd, xcd, ycd)] = smarr[i] + pfx_ofst;
            }
            __syncthreads();
            //update enc_ofst and rem_val
            if(threadIdx.x == 0) {
                pfx_ofst += pfx + tpfx_ofst;
                // pfx_ofst += smarr[min(mem_size, rem_val) - 1];
                enc_ofst += mem_size;
                rem_val -= mem_size;
            }
            __syncthreads();
        }
    }
    __global__ void d_init_bit_cost(int max_size) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_size) {
            return;
        }
        d_cost_data_bit[tid] = get_gpu_of_cost(d_cube_data_bit[tid]);
    }
    void write_integer(int num) {
        if(num == 0) {
            putc_unlocked('0', OUT_FILE);
        } else {
            static char temp[20];
            int len = 0;
            while(num) temp[len++] = num % 10, num /= 10;
            for(int i = len - 1; i >= 0; i--) putc_unlocked('0' + temp[i], OUT_FILE);
        }
    }
    void add_vias_cpu(long x, long y, long z1, long z2) {
      if(z1 > z2) {swap(z1, z2);}
      for(int i = z1;i < z2; i++) {
        write_integer(x);
        putc_unlocked(' ', OUT_FILE);
        write_integer(y);
        putc_unlocked(' ', OUT_FILE);
        putc_unlocked('m', OUT_FILE);
        putc_unlocked('e', OUT_FILE);
        putc_unlocked('t', OUT_FILE);
        putc_unlocked('a', OUT_FILE);
        putc_unlocked('l', OUT_FILE);
        write_integer(i + 1);
        putc_unlocked(' ', OUT_FILE);
        write_integer(x);
        putc_unlocked(' ', OUT_FILE);
        write_integer(y);
        putc_unlocked(' ', OUT_FILE);
        putc_unlocked('m', OUT_FILE);
        putc_unlocked('e', OUT_FILE);
        putc_unlocked('t', OUT_FILE);
        putc_unlocked('a', OUT_FILE);
        putc_unlocked('l', OUT_FILE);
        write_integer(i + 2);
        putc_unlocked('\n', OUT_FILE);
      }
    }
    void add_horizontal_cpu(long z, long y, long x1, long x2) {    
        if (x1 != x2) {
            write_integer(x1);
            putc_unlocked(' ', OUT_FILE);
            write_integer(y);
            putc_unlocked(' ', OUT_FILE);
            putc_unlocked('m', OUT_FILE);
            putc_unlocked('e', OUT_FILE);
            putc_unlocked('t', OUT_FILE);
            putc_unlocked('a', OUT_FILE);
            putc_unlocked('l', OUT_FILE);
            write_integer(z + 1);
            putc_unlocked(' ', OUT_FILE);
            write_integer(x2);
            putc_unlocked(' ', OUT_FILE);
            write_integer(y);
            putc_unlocked(' ', OUT_FILE);
            putc_unlocked('m', OUT_FILE);
            putc_unlocked('e', OUT_FILE);
            putc_unlocked('t', OUT_FILE);
            putc_unlocked('a', OUT_FILE);
            putc_unlocked('l', OUT_FILE);
            write_integer(z + 1);
            putc_unlocked('\n', OUT_FILE);
        }
    }
    void add_vertical_cpu(long z, long x, long y1, long y2) {
      if (y1 != y2) {
        write_integer(x);
        putc_unlocked(' ', OUT_FILE);
        write_integer(y1);
        putc_unlocked(' ', OUT_FILE);
        putc_unlocked('m', OUT_FILE);
        putc_unlocked('e', OUT_FILE);
        putc_unlocked('t', OUT_FILE);
        putc_unlocked('a', OUT_FILE);
        putc_unlocked('l', OUT_FILE);
        write_integer(z + 1);
        putc_unlocked(' ', OUT_FILE);
        write_integer(x);
        putc_unlocked(' ', OUT_FILE);
        write_integer(y2);
        putc_unlocked(' ', OUT_FILE);
        putc_unlocked('m', OUT_FILE);
        putc_unlocked('e', OUT_FILE);
        putc_unlocked('t', OUT_FILE);
        putc_unlocked('a', OUT_FILE);
        putc_unlocked('l', OUT_FILE);
        write_integer(z + 1);
        putc_unlocked('\n', OUT_FILE);
      }
    }
    __inline__ void write_net_name(string &nname) {
        for(int i = 0; i < nname.size(); i++) {
            putc_unlocked(nname[i], OUT_FILE);
        }
        putc_unlocked('\n', OUT_FILE);
        putc_unlocked('(', OUT_FILE);
        putc_unlocked('\n', OUT_FILE);
    }
    __inline__ void write_closed_bracket() {
        putc_unlocked(')', OUT_FILE);
        putc_unlocked('\n', OUT_FILE);
    }
    __inline__ long get_gran_cd(long cd) {
        return 4200L * cd + 2100L;
    }
    void pfx_cube_cost() {
        TIME_();
        int enc_id = 0;
        for(gvia_dt z = 0; z < Z; z++) {
            if(z & 1) {
                for(gcd_dt x = 0; x < X; x++) {
                    gcb_dt pfx = 0, prv = 0;
                    for(gcd_dt y = 0; y < Y; y++) {
                        enc_id = UTIL::get_enc<WIRE_ENC>(z, x, y);
                        prv=wcost_cube[enc_id]; wcost_cube[enc_id] = pfx; enc_id++;
                        pfx += prv;
                    }
                }
            } else {
                for(int y = 0; y < Y; y++) {
                    gcb_dt pfx = 0, prv = 0;
                    for(gcd_dt x = 0; x < X; x++) {
                        enc_id = UTIL::get_enc<WIRE_ENC>(z, x, y);
                        prv = wcost_cube[enc_id]; wcost_cube[enc_id] = pfx; enc_id++;
                        pfx += prv;
                    }
                }
            }
        }
    }     
    void get_cost_data() {
        int nm_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
        UTIL::get_gpu_cost_cube<<<nm_blks, NUM_GPU_THREADS>>>(d_wcost_cube, d_vcost_cube, d_cube_data, CSIZE);
        cudaDeviceSynchronize();
    }
    void init_bit_cost() {
        int nm_blks = CEIL_DIV(CSIZE, NUM_GPU_THREADS);
        if(BIT_PMODEL) {
            spdlog::info("initialising bit cost");
            d_init_bit_cost<<<nm_blks, NUM_GPU_THREADS>>>(CSIZE);
        }
        cudaDeviceSynchronize();
    }
    __global__ void copy_deg_net_len(int max_nets) {
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        if(tid >= max_nets) {
            return;
        }
        d_net_pins_len[tid] = d_net_deg_data[tid];
    }
    void print_gpu_free_mem() {
        size_t free_byte, total_byte;
        cudaError_t cuda_status = cudaMemGetInfo(&free_byte, &total_byte);
        SPDLOG_ASSERT(cuda_status == cudaSuccess, "{} {}", cuda_status, cudaSuccess);
        double free_gb = static_cast<double>(free_byte) / (double(1L << 30));
        double total_gb = static_cast<double>(total_byte) / (double(1L << 30));
        double used_gb = total_gb - free_gb;
        spdlog::info("[free_gb,used_gb,total_gb] [{},{},{}]", free_gb, used_gb, total_gb);
    }
}
