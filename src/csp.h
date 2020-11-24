/* Utils
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_CSP_H
#define CSP_CSP_H

#include <stdlib.h>
#include "htslib/sam.h"

/*
* File API
 */

/*@abstract  Packing the common bam file related pointers into a structure. */
typedef struct {
    htsFile *fp;
    sam_hdr_t *hdr;   // hdr is needed by sam_read1().
    hts_idx_t *idx;
} csp_bam_fs;

/*@abstract  Create a csp_bam_fs structure.
@return  Pointer to the csp_bam_fs structure if success, NULL otherwise.

@note    The pointer returned successfully by csp_bam_fs_init() should be freed
         by csp_bam_fs_destroy() when no longer used.
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

/*@abstract  Build the csp_bam_fs structure.  
@param fn    Filename of the bam/sam/cram.
@param ret   Pointer to the state.
             0 if success, -1 for memory error, -2 for open/parse sam file error.
@return      The pointer to the csp_bam_fs structure if success, NULL otherwise.

@note        The pointer returned successfully by csp_bam_fs_build() should be freed
             by csp_bam_fs_destroy() when no longer used.
*/
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

#endif
