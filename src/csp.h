/* Utils
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_CSP_H
#define CSP_CSP_H

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
static inline csp_bam_fs* csp_bam_fs_init(void);

static inline void csp_bam_fs_destroy(csp_bam_fs *p);

/*@abstract  Build the csp_bam_fs structure.  
@param fn    Filename of the bam/sam/cram.
@param ret   Pointer to the state.
             0 if success, -1 for memory error, -2 for open/parse sam file error.
@return      The pointer to the csp_bam_fs structure if success, NULL otherwise.

@note        The pointer returned successfully by csp_bam_fs_build() should be freed
             by csp_bam_fs_destroy() when no longer used.
*/
static inline csp_bam_fs* csp_bam_fs_build(const char *fn, int *ret);

static inline int csp_bam_fs_reset(csp_bam_fs *p, const char *fn);

#endif
