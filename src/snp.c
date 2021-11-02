/* snp.c - SNP operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "kvec.h"
#include "jstring.h"
#include "snp.h"

/* 
* SNP List API
*/

snp_t* snp_init(void) { return (snp_t*) calloc(1, sizeof(snp_t)); }

void snp_destroy(snp_t *p) { 
    if (p) { free(p->chr); free(p); } 
}

void snp_reset(snp_t *p) {
    if (p) { free(p->chr); memset(p, 0, sizeof(snp_t)); }
}

size_t get_snplist_from_vcf(const char *fn, snplist_t *pl, int *ret, int print_skip) {
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    snp_t *ip = NULL;
    size_t l, m, n = 0;  /* use n instead of kv_size(*pl) in case that pl is not empty at the beginning. */
    int r;
    *ret = -1;
    if (NULL == fn || NULL == pl) 
        return 0;
    if (NULL == (fp = hts_open(fn, "rb"))) {
        fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, fn);
        return 0;
    }
    if (NULL == (hdr = bcf_hdr_read(fp))) {
        fprintf(stderr, "[E::%s] could not read header for '%s'\n", __func__, fn);
        goto fail;
    }
    if (NULL == (rec = bcf_init1())) {
        fprintf(stderr, "[E::%s] could not initialize the bcf structure.\n", __func__);
        goto fail;
    }
    for (m = 1; (r = bcf_read1(fp, hdr, rec)) >= 0; m++) {
        if (NULL == (ip = snp_init())) { 
            fprintf(stderr, "[E::%s] could not initialize the snp_t structure.\n", __func__); 
            goto fail; 
        }
        ip->chr = safe_strdup(bcf_hdr_id2name(hdr, rec->rid));
        if (NULL == ip->chr) {
            if (print_skip)
                fprintf(stderr, "[W::%s] skip No.%ld SNP: could not get chr name.\n", __func__, m);
            snp_destroy(ip);
            continue;
        }
        ip->pos = rec->pos;
        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele > 0) {
            if (1 == (l = strlen(rec->d.allele[0]))) 
                ip->ref = rec->d.allele[0][0];
            else if (l > 1) { 
                if (print_skip) 
                    fprintf(stderr, "[W::%s] skip No.%ld SNP: ref_len > 1.\n", __func__, m);
                snp_destroy(ip); 
                continue; 
            } // else: do nothing. keep ref = 0, 0 is its init value.               
            if (2 == rec->n_allele) {
                if (1 == (l = strlen(rec->d.allele[1]))) 
                    ip->alt = rec->d.allele[1][0];
                else if (l > 1) {
                    if (print_skip) 
                        fprintf(stderr, "[W::%s] skip No.%ld SNP: alt_len > 1.\n", __func__, m);
                    snp_destroy(ip);
                    continue; 					
                } // else: do nothing. keep alt = 0, 0 is its init value.
            } else if (rec->n_allele > 2) {
                if (print_skip)
                    fprintf(stderr, "[W::%s] skip No.%ld SNP: n_allele > 2.\n", __func__, m);
                snp_destroy(ip);
                continue;                 
            } // else: keep alt = 0.
        } // else: do nothing. keep ref = alt = 0, 0 is their init values.
        kv_push(snp_t*, *pl, ip);
        n++;
    }
    if (-1 == r) { // end of bcf file.
        kv_resize(snp_t*, *pl, kv_size(*pl));
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

/*
 * Bi-Allele API
 */
biallele_t* biallele_init(void) { return (biallele_t*) calloc(1, sizeof(biallele_t)); }
void biallele_destroy(biallele_t *p) { if (p) { free(p); } }
void biallele_reset(biallele_t *p) {}

