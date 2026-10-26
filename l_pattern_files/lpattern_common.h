#pragma once
#include "../common.h"

using namespace std;
using namespace COMMON_DATA;

typedef int16_t lpat_dt;
typedef long ledid_dt;

#define ADD_LPATTERN_BATCH_DIST 0

#define LCPU_GPU_COMPARISON 0
#define LTHREAD_PER_EDGE 50
#define TEST_LPATTERN_GEN 0

#define LGPU_ASSERT 1
#define LCPU_ASSERT 1

#define LPAT_LOGS 0
#define LPAT_DETAIL_LOGS 0

#define TEST_GPU_PFX_LOGS 0
#define TEST_GPU_PFX_LOGS_DETAIL 0
#define TEST_LPAT_ADD_EDGE 1
#define LCPU_GPU_EDGEADD_COMPARISON 0
#define LWIRE_FACTOR 1
#define LVIA_FACTOR 1
#define COMMIT_LPAT 0

namespace LPAT_COMMON {
    struct dcd_t {
        gcd_dt x, y;
        gvia_dt z;
        gcb_dt dist;
        int cdir, ndir_chng;
    };
    int nx[6] = {-1, 0, 1, 0, 0, 0};
    int ny[6] = {0, 1, 0, -1, 0, 0};
    int nz[6] = {0, 0, 0, 0, 1, -1};
    gcb_dt *lpat_dist;
    int *lpat_par;
    dcd_t *djkst_oval;
    string LPAT_BATCH_QUERY = "";
    __managed__ lpat_dt *d_lpat_pat;
    lpat_dt *hd_lpat_pat, *lpat_pat;
    gcb_dt *lpat_res, *hd_lpat_res;
    genc_dt *d_ptprs, *hd_ptprs;
    gcb_dt *d_lpat_res;
    ledid_dt *d_led_id;
    pth_dt *d_pth_mat, *hd_pth_mat;
    gcd_dt *pth_mat;
}
