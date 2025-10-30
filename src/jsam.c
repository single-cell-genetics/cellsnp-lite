// jsam.c - sequence alignment file operations.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "jsam.h"

/* 
* BAM/SAM/CRAM File API 
*/
const char csp_nt5_str[] = "ACGTN";

/*@note If the translation failed, this function will try "<name> - chr" if name starts with "chr", try
        "chr + <name>" otherwise. */
int csp_sam_hdr_name2id(sam_hdr_t *hdr, const char *name, kstring_t *s) {
    int tid;
    if ((tid = sam_hdr_name2tid(hdr, name)) < 0) {
        if (tid < -1) { return tid; }
        if (0 == strncmp(name, "chr", 3)) { return sam_hdr_name2tid(hdr, name + 3); }
        else {
            kputs("chr", s); kputs(name, s);
            return sam_hdr_name2tid(hdr, ks_str(s));
        }
    } else { return tid; }
}

//@note No need to free the returned char* pointer when success.
const char* csp_fmt_chr_name(const char *name, sam_hdr_t *hdr, kstring_t *s) {
    int tid;
    if ((tid = csp_sam_hdr_name2id(hdr, name, s)) < 0) { return NULL; }
    else { return sam_hdr_tid2name(hdr, tid); }
}

/*@note 1. To speed up, the caller should guarantee parameters b and tag are valid. 
        2. The data of the pointer returned by this function is part of bam1_t, so do not double free!
 */
char* get_bam_aux_str(bam1_t *b, const char tag[2]) {
    uint8_t *data;
    if (NULL == (data = bam_aux_get(b, tag))) { return NULL; }
    return bam_aux2Z(data);
}