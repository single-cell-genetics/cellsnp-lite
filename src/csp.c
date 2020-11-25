/* Utils
 * Author: Xianjie Huang <hxj5@hku.hk>
 */

#include <stdlib.h>
#include "htslib/sam.h"
#include "csp.h"

/*
* File API
 */

static inline csp_bam_fs* csp_bam_fs_init(void) { return (csp_bam_fs*) calloc(1, sizeof(csp_bam_fs)); }

static inline void csp_bam_fs_destroy(csp_bam_fs *p) {
    if (p) {
        if (p->idx) { hts_idx_destroy(p->idx); }
        if (p->hdr) { sam_hdr_destroy(p->hdr); }
        if (p->fp)  { hts_close(p->fp); }
        free(p);
    }
}

static inline csp_bam_fs* csp_bam_fs_build(const char *fn, int *ret) {
    csp_bam_fs *p;
    if (NULL == fn) { *ret = -1; return NULL; }
    if (NULL == (p = csp_bam_fs_init())) { *ret = -1; return NULL; }
    if (NULL == (p->fp = hts_open(fn, "rb"))) { *ret = -2; goto fail; }
    if (NULL == (p->hdr = sam_hdr_read(p->fp))) { *ret = -2; goto fail; }
    if (NULL == (p->idx = sam_index_load(p->fp, fn))) { *ret = -2; goto fail; }
    *ret = 0;
    return p;
  fail:
    csp_bam_fs_destroy(p);
    return NULL;		
}

static inline int csp_bam_fs_reset(csp_bam_fs *p, const char *fn) {
    if (NULL == p) { return -1; }
    if (p->idx) { hts_idx_destroy(p->idx); }
    if (p->hdr) { sam_hdr_destroy(p->hdr); }
    if (p->fp)  { hts_close(p->fp); }		
    if (NULL == (p->fp = hts_open(fn, "rb"))) { return -1; }
    if (NULL == (p->hdr = sam_hdr_read(p->fp))) { return -1; }
    if (NULL == (p->idx = sam_index_load(p->fp, fn))) { return -1; }
    return 0;
}

