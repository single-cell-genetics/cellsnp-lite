/* SNP operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_SNP_H
#define CSP_SNP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "kvec.h"

/* 
* SNP List API
*/
/*@abstruct    The basic data structure to store SNP info.
@param chr     Name of chr.
@param pos     0-based coordinate in the reference sequence.
@param ref     Ref base (a letter). 0 means no ref in the input SNP file for the pos.
@param alt     Alt base (a letter). 0 means no alt in the input SNP file for the pos.
 */
typedef struct {
    char *chr;   
    hts_pos_t pos; 
    int8_t ref, alt;
} csp_snp_t;

/*@abstract  Initilize the csp_snp_t structure.
@return      Pointer to the csp_snp_t structure if success, NULL if failure.
@note        The pointer returned successfully by csp_snp_init() should be freed
             by csp_snp_destroy() when no longer used.
 */
static inline csp_snp_t* csp_snp_init(void) { return (csp_snp_t*) calloc(1, sizeof(csp_snp_t)); }

static inline void csp_snp_destroy(csp_snp_t *p) { 
    if (p) { free(p->chr); free(p); } 
}

static inline void csp_snp_reset(csp_snp_t *p) {
    if (p) { free(p->chr); memset(p, 0, sizeof(csp_snp_t)); }
}

/*@abstract  A list containing of several pointers to the csp_snp_t structure.
@param a     Array of csp_snp_t* pointers.
@param n     The next pos of the unused element of the array.
@param m     Size of the whole array.

@note        The csp_snplist_t structure should be freed by csp_snplist_destroy() when no longer used.
             csp_snplist_init() function should be called immediately after the structure was created.

@example (A simple example from kvec.h)
    kvec_t(int) array;
    kv_init(array);
    kv_push(int, array, 10); // append
    kv_A(array, 20) = 4;
    kv_destroy(array);
 */
typedef kvec_t(csp_snp_t*) csp_snplist_t;   /* kvec_t from kvec.h */
#define csp_snplist_init(v) kv_init(v)
#define csp_snplist_resize(v, size) kv_resize(csp_snp_t*, v, size)
#define csp_snplist_push(v, x) kv_push(csp_snp_t*, v, x)
#define csp_snplist_A(v, i) kv_A(v, i)
#define csp_snplist_size(v) kv_size(v)
#define csp_snplist_max(v) kv_max(v)
#define csp_snplist_destroy(v) {								\
    size_t __j;											\
    for (__j = 0; __j < csp_snplist_size(v); __j++) csp_snp_destroy(csp_snplist_A(v, __j));	\
    kv_destroy(v);										\
}

/*@abstract    Extract SNP info from bcf/vcf file.
@param fn      Filename of bcf/vcf.
@param pl      Pointer to client data used to store the extracted SNP info.
@param ret     Pointer to store running state. 0 if success, -1 otherwise.
@param print_skip  If print the skipped SNPs. 1, yes, 0, no.
@return        Num of elements successfully added into the snplist.

@note          If length of Ref or Alt is larger than 1, then the SNP would be skipped.
               If length of Ref or Alt is 0, then their values would be infered during pileup.
 */
static size_t get_snplist_from_vcf(const char *fn, csp_snplist_t *pl, int *ret, int print_skip) {
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    csp_snp_t *ip = NULL;
    size_t l, m, n = 0;  /* use n instead of kv_size(*pl) in case that pl is not empty at the beginning. */
    int r;
    *ret = -1;
    if (NULL == fn || NULL == pl) { return 0; }
    if (NULL == (fp = hts_open(fn, "rb"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, fn); return 0; }
    if (NULL == (hdr = bcf_hdr_read(fp))) { fprintf(stderr, "[E::%s] could not read header for '%s'\n", __func__, fn); goto fail; }
    if (NULL == (rec = bcf_init1())) { fprintf(stderr, "[E::%s] could not initialize the bcf structure.\n", __func__); goto fail; }
    for (m = 1; (r = bcf_read1(fp, hdr, rec)) >= 0; m++) {
        if (NULL == (ip = csp_snp_init())) { 
            fprintf(stderr, "[E::%s] could not initialize the csp_snp_t structure.\n", __func__); 
            goto fail; 
        }
        ip->chr = safe_strdup(bcf_hdr_id2name(hdr, rec->rid));
        if (NULL == ip->chr) {
            if (print_skip) { fprintf(stderr, "[W::%s] skip No.%ld SNP: could not get chr name.\n", __func__, m); }
            csp_snp_destroy(ip);
            continue;
        } ip->pos = rec->pos;
        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele > 0) {
            if (1 == (l = strlen(rec->d.allele[0]))) { ip->ref = rec->d.allele[0][0]; }
            else if (l > 1) { 
                if (print_skip) { fprintf(stderr, "[W::%s] skip No.%ld SNP: ref_len > 1.\n", __func__, m); }
                csp_snp_destroy(ip); 
                continue; 
            } // else: do nothing. keep ref = 0, 0 is its init value.               
            if (2 == rec->n_allele) {
                if (1 == (l = strlen(rec->d.allele[1]))) { ip->alt = rec->d.allele[1][0]; }
                else if (l > 1) {
                    if (print_skip) { fprintf(stderr, "[W::%s] skip No.%ld SNP: alt_len > 1.\n", __func__, m); }
                    csp_snp_destroy(ip);
                    continue; 					
                } // else: do nothing. keep alt = 0, 0 is its init value.
            } else if (rec->n_allele > 2) {
                if (print_skip) { fprintf(stderr, "[W::%s] skip No.%ld SNP: n_allele > 2.\n", __func__, m); }
                csp_snp_destroy(ip);
                continue;                 
            } // else: keep alt = 0.
        } // else: do nothing. keep ref = alt = 0, 0 is their init values.
        csp_snplist_push(*pl, ip);
        n++;
    }
    if (-1 == r) { // end of bcf file.
        csp_snplist_resize(*pl, csp_snplist_size(*pl));
    } else { 
        fprintf(stderr, "[E::%s] error when parsing '%s'\n", __func__, fn); 
        goto fail; 
    }
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    *ret = 0;
    return n;
  fail:
    if (rec) bcf_destroy1(rec);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp) hts_close(fp);
    return n;
}
#define get_snplist(fn, pl, ret, print_skip) get_snplist_from_vcf(fn, pl, ret, print_skip)

#endif
