// mplp.c - pileup and mpileup operations.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "kvec.h"
#include "jnumeric.h"
#include "jfile.h"
#include "jmempool.h"
#include "config.h"
#include "mplp.h"

/*
* Pileup and MPileup API
 */

csp_pileup_t* csp_pileup_init(void) {
    csp_pileup_t *p = (csp_pileup_t*) malloc(sizeof(csp_pileup_t));
    if (p) {
        if (NULL == (p->b = bam_init1())) { free(p); return NULL; }
    }
    return p;
}

void csp_pileup_destroy(csp_pileup_t *p) { 
    if (p) {
        if (p->b) bam_destroy1(p->b);	
        free(p);
    } 
}

int csp_pileup_reset(csp_pileup_t *p) {
    if (p) {
        if (p->b) { bam_destroy1(p->b); }
        memset(p, 0, sizeof(csp_pileup_t)); 
        if (NULL == (p->b = bam_init1())) { return -1; }
    } // else: do nothing.
    return 0;
}

void csp_pileup_reset_(csp_pileup_t *p) { }

void csp_pileup_print(FILE *fp, csp_pileup_t *p) {
    fprintf(fp, "qpos = %d\n", p->qpos);
    fprintf(fp, "base = %c, qual = %d\n", p->base, p->qual);
    fprintf(fp, "is_refskip = %d, is_del = %d\n", p->is_refskip, p->is_del);
    fprintf(fp, "umi = %s, cb = %s\n", p->umi, p->cb);
    fprintf(fp, "len_aln = %d\n", p->laln);
}

umi_unit_t* umi_unit_init(void) {
    umi_unit_t *p = (umi_unit_t*) calloc(1, sizeof(umi_unit_t));
    return p;   /* will set values just after this function is called so no need to set init values here. */
}
void umi_unit_destroy(umi_unit_t *p) { free(p); }

list_uu_t* list_uu_init(void) {
    list_uu_t *v = (list_uu_t*) malloc(sizeof(list_uu_t));
    if (v) { kv_init(*v); }
    return v;
}

JNUMERIC_INIT(cu_d, double)
JNUMERIC_INIT(cu_s, size_t)

int get_qual_vector(double qual, double cap_bq, double min_bq, double *rv) {
    double bq = max2(min2(cap_bq, qual), min_bq);
    double p = pow(0.1, bq / 10);
    rv[0] = log(1 - p);              rv[1] = log(0.75 - 2.0 / 3 * p);
    rv[2] = log(0.5 - 1.0 / 3 * p);  rv[3] = log(p);
    return 0;
}

int qual_matrix_to_geno(double qm[][4], size_t *bc, int8_t ref_idx, int8_t alt_idx, int db, double *gl, int *n) {
    int8_t other_idx[4], noth, i;
    size_t ref_read, alt_read;
    double *ref_qual, *alt_qual, oth_qual, tmp_qual;
    for (i = 0, noth = 0; i < 5; i++) {
        if (i != ref_idx && i != alt_idx) other_idx[noth++] = i;
    }
    ref_read = bc[ref_idx];
    alt_read = bc[alt_idx];
    ref_qual = qm[ref_idx];
    alt_qual = qm[alt_idx];
    for (i = 0, oth_qual = 0, tmp_qual = 0; i < noth; i++) {
        oth_qual += qm[other_idx[i]][3];
        tmp_qual += bc[other_idx[i]];
    }
    oth_qual += log(2.0 / 3) * tmp_qual;
    gl[0] = oth_qual + ref_qual[0] + alt_qual[3] + log(1.0 / 3) * alt_read;
    gl[1] = oth_qual + ref_qual[2] + alt_qual[2];
    gl[2] = oth_qual + ref_qual[3] + alt_qual[0] + log(1.0 / 3) * ref_read;
    if (db) {
        gl[3] = oth_qual + ref_qual[1] + log(1.0 / 4) * alt_read;
        gl[4] = oth_qual + alt_qual[1] + log(1.0 / 4) * ref_read;
        *n = 5;
    } else { *n = 3; }
    return 0;
}

//@note It's usually called when the input pos has no ref or alt.
void infer_allele(size_t *bc, int8_t *ref_idx, int8_t *alt_idx) {
    int8_t i, k1, k2;
    size_t m1, m2;
    if (bc[0] < bc[1]) { m1 = bc[1]; m2 = bc[0]; k1 = 1; k2 = 0; }
    else { m1 = bc[0]; m2 = bc[1]; k1 = 0; k2 = 1; }
    for (i = 2; i < 5; i++) {
        if (bc[i] > m1) { m2 = m1; k2 = k1; m1 = bc[i]; k1 = i; }
        else if (bc[i] > m2) { m2 = bc[i]; k2 = i; }
    }
    *ref_idx = k1; *alt_idx = k2;
}

void infer_alt(size_t *bc, int8_t ref_idx, int8_t *alt_idx) {
    int8_t i, k;
    size_t m = 0;   
    for (i = 0, k = -1; i < 5; i++) {
        if (ref_idx == i) { continue; }
        if (bc[i] > m) { k = i; m = bc[i]; }
    }
    if (-1 == k) { *alt_idx = ref_idx == 0 ? 1 : 0; }
    else { *alt_idx = k; }
}

/* note that the @p qu is also initialized after calling calloc(). */
csp_plp_t* csp_plp_init(void) { return (csp_plp_t*) calloc(1, sizeof(csp_plp_t)); }

void csp_plp_destroy(csp_plp_t *p) { 
    if (p) { 
        int i;
        for (i = 0; i < 5; i++) { kv_destroy(p->qu[i]); }
        if (p->hug) { map_ug_destroy(p->hug); }
        free(p); 
    }
}

void csp_plp_reset(csp_plp_t *p) {
    if (p) {   // TODO: reset based on is_genotype.
        int i;
        memset(p->bc, 0, sizeof(p->bc));
        p->tc = p->ad = p->dp = p->oth = 0;
        for (i = 0; i < 5; i++) { list_qu_reset(p->qu[i]); }
        memset(p->qmat, 0, sizeof(p->qmat));
        p->ngl = 0;
        if (p->hug) { map_ug_reset(p->hug); }
    }
}

void csp_plp_print(FILE *fp, csp_plp_t *p, char *prefix) {
    int i, j;
    khiter_t u;
    fprintf(fp, "%stotal read count = %ld\n", prefix, p->tc);
    fprintf(fp, "%sbase count (A/C/G/T/N):", prefix);
    for (i = 0; i < 5; i++)
        fprintf(fp, " %ld", p->bc[i]);
    fputc('\n', fp);
    fprintf(fp, "%squal matrix 5x4:\n", prefix);
    for (i = 0; i < 5; i++) {
        fprintf(fp, "%s\t", prefix);
        for (j = 0; j < 4; j++)
            fprintf(fp, " %.2f", p->qmat[i][j]);
        fputc('\n', fp);
    }
    fprintf(fp, "%snum of geno likelihood = %d\n", prefix, p->ngl);
    if (p->ngl) {
        fprintf(fp, "%sgeno likelihood:", prefix);
        for (i = 0; i < p->ngl; i++)
            fprintf(fp, " %.2f", p->gl[i]);
        fputc('\n', fp);
    }
    if (p->hug) {
        int size = kh_size(p->hug);
        fprintf(fp, "%ssize of the map_ug = %d\n", prefix, size);
        if (size) {
            fprintf(fp, "%s", prefix);
            for (u = kh_begin(p->hug); u != kh_end(p->hug); u++) {
                if (kh_exist(p->hug, u))
                    fprintf(fp, " %s", kh_key(p->hug, u));
            }
            fputc('\n', fp);
        }
    }
}

int csp_plp_str_vcf(csp_plp_t *p, kstring_t *s) {
    if (p->tc <= 0) {
        kputs(".:.:.:.:.:.", s);
        return 0;
    }
    double gl[5];
    int i, m;
    double tmp = -10 / log(10);
    char *gt[] = {"0/0", "1/0", "1/1"};
    m = get_idx_of_max(cu_d, p->gl, 3);
    kputs(gt[m], s);
    ksprintf(s, ":%ld:%ld:%ld:", p->ad, p->dp, p->oth);
    for (i = 0; i < p->ngl; i++) 
        gl[i] = p->gl[i] * tmp;
    if (join_arr_to_str(cu_d, gl, p->ngl, ',', "%.0f", s) < p->ngl) 
        return -1;
    kputc_(':', s);
    if (join_arr_to_str(cu_s, p->bc, 5, ',', "%ld", s) < 5) 
        return -1;
    return 0;
}

int csp_plp_to_vcf(csp_plp_t *p, jfile_t *s) {
    if (p->tc <= 0) {
        jf_puts(".:.:.:.:.:.", s);
        return 0;
    }
    double gl[5];
    int i, m;
    double tmp = -10 / log(10);
    char *gt[] = {"0/0", "1/0", "1/1"};
    m = get_idx_of_max(cu_d, p->gl, 3);
    jf_puts(gt[m], s);
    jf_printf(s, ":%ld:%ld:%ld:", p->ad, p->dp, p->oth);
    for (i = 0; i < p->ngl; i++)
        gl[i] = p->gl[i] * tmp;
    if (join_arr_to_str(cu_d, gl, p->ngl, ',', "%.0f", s->buf) < p->ngl)
        return -1;  // TODO: use internal buf directly is not good.
    jf_putc_(':', s);
    if (join_arr_to_str(cu_s, p->bc, 5, ',', "%ld", s->buf) < 5)
        return -1;
    return 0;
}

csp_mplp_t* csp_mplp_init(void) { 
    csp_mplp_t *p = (csp_mplp_t*) calloc(1, sizeof(csp_mplp_t));
    return p;
}

void csp_mplp_destroy(csp_mplp_t *p) { 
    if (p) {
        if (p->hsg) { map_sg_destroy(p->hsg); }
        if (p->hsg_iter) { free(p->hsg_iter); }
        if (p->pu) { jmempool_destroy(umi_unit, p->pu); }
        if (p->pl) { jmempool_destroy(list_umiunit, p->pl); }
        if (p->su) { jmempool_destroy(str, p->su); }
        if (p->fai) { fai_destroy(p->fai); }
        free(p); 
    }
}

void csp_mplp_reset(csp_mplp_t *p) {
    if (p) {
        memset(p->bc, 0, sizeof(p->bc));
        p->tc = p->ad = p->dp = p->oth = 0;
        p->nr_ad = p->nr_dp = p->nr_oth = 0;
        if (p->hsg) { map_sg_reset_val(p->hsg); }
        if (p->pu) { jmempool_reset(umi_unit, p->pu); }
        if (p->pl) { jmempool_reset(list_umiunit, p->pl); }
        if (p->su) { jmempool_reset(str, p->su); }
        memset(p->qvec, 0, sizeof(p->qvec));
    }
}

void csp_mplp_print(FILE *fp, csp_mplp_t *p, char *prefix) {
    int i;
    csp_plp_t *plp;
    kstring_t ks = KS_INITIALIZE, *s = &ks;
    fprintf(fp, "%sref_idx = %d, alt_idx = %d\n", prefix, p->ref_idx, p->alt_idx);
    fprintf(fp, "%sinf_rid = %d, inf_aid = %d\n", prefix, p->inf_rid, p->inf_aid);
    fprintf(fp, "%stotal base count = %ld\n", prefix, p->tc);
    fprintf(fp, "%sbase count (A/C/G/T/N):", prefix);
    for (i = 0; i < 5; i++) 
        fprintf(fp, " %ld", p->bc[i]);
    fputc('\n', fp);
    fprintf(fp, "%snum of sample group = %d\n", prefix, p->nsg);
    if (p->nsg) {
        kputs(prefix, s); kputc('\t', s);
        for (i = 0; i < p->nsg; i++) {
            fprintf(fp, "%sSG-%d = %s:\n", prefix, i, kh_key(p->hsg, p->hsg_iter[i]));
            plp = kh_val(p->hsg, p->hsg_iter[i]);
            csp_plp_print(fp, plp, ks_str(s));
        }
    }
    ks_free(s);
}

void csp_mplp_print_(FILE *fp, csp_mplp_t *p, char *prefix) {
    int i;
    fprintf(fp, "%sref_idx = %d, alt_idx = %d\n", prefix, p->ref_idx, p->alt_idx);
    fprintf(fp, "%sinf_rid = %d, inf_aid = %d\n", prefix, p->inf_rid, p->inf_aid);
    fprintf(fp, "%stotal base count = %ld\n", prefix, p->tc);
    fprintf(fp, "%sbase count (A/C/G/T/N):", prefix);
    for (i = 0; i < 5; i++) 
        fprintf(fp, " %ld", p->bc[i]);
    fputc('\n', fp);
    fprintf(fp, "%snum of sample group = %d\n", prefix, p->nsg);
}

int csp_mplp_set_sg(csp_mplp_t *p, char **s, const int n) {
    if (NULL == p || NULL == s || 0 == n) 
        return -1;
    int i, r;
    khiter_t k;
    if (NULL == p->hsg && NULL == (p->hsg = kh_init(sample_group))) 
        return -1;
    if (NULL == p->hsg_iter && NULL == (p->hsg_iter = (khiter_t*) malloc(sizeof(khiter_t) * n))) 
        return -1;
    for (i = 0; i < n; i++) {
        if (s[i]) { 
            k = kh_put(sample_group, p->hsg, s[i], &r); 
            if (r > 0) 
                kh_val(p->hsg, k) = NULL;
            else if (r < 0) 
                return -1;
            else
                return -2;  /* r = 0 means repeatd sgnames. */
        } else {
            return -1;
        }
    }
    // Storing iter index for each sg (sample group) name must be done after all sg names have been 
    // pushed into the HashMap in case that the internal arrays of HashMap autoly shrink or some 
    // else modifications.
    for (i = 0; i < n; i++) {
        k = kh_get(sample_group, p->hsg, s[i]);
        if (k == kh_end(p->hsg))
            return -1;
        else
            p->hsg_iter[i] = k;
    }
    p->nsg = n;
    return 0;
}

int csp_mplp_str_vcf(csp_mplp_t *mplp, kstring_t *s) {
    int i;
    for (i = 0; i < mplp->nsg; i++) {
        kputc_('\t', s);
        if (csp_plp_str_vcf(kh_val(mplp->hsg, mplp->hsg_iter[i]), s) < 0) 
            return -1; 
    } //s->s[--(s->l)] = '\0';    /* s->l could not be negative unless no csp_plp_t(s) are printed to s->s. */
    return 0;
}

int csp_mplp_to_vcf(csp_mplp_t *mplp, jfile_t *s) {
    int i;
    for (i = 0; i < mplp->nsg; i++) {
        jf_putc_('\t', s);
        if (csp_plp_to_vcf(kh_val(mplp->hsg, mplp->hsg_iter[i]), s) < 0) 
            return -1;
    } //s->s[--(s->l)] = '\0';    /* s->l could not be negative unless no csp_plp_t(s) are printed to s->s. */
    return 0;
}

int csp_mplp_str_mtx(csp_mplp_t *mplp, kstring_t *ks_ad, kstring_t *ks_dp, kstring_t *ks_oth, size_t idx) {
    csp_plp_t *plp;
    int i;
    for (i = 1; i <= mplp->nsg; i++) {
        plp = kh_val(mplp->hsg, mplp->hsg_iter[i]);
        if (plp->ad) ksprintf(ks_ad, "%ld\t%d\t%ld\n", idx, i, plp->ad);
        if (plp->dp) ksprintf(ks_dp, "%ld\t%d\t%ld\n", idx, i, plp->dp);
        if (plp->oth) ksprintf(ks_oth, "%ld\t%d\t%ld\n", idx, i, plp->oth);        
    }
    return 0; 
}

int csp_mplp_str_mtx_tmp(csp_mplp_t *mplp, kstring_t *ks_ad, kstring_t *ks_dp, kstring_t *ks_oth) {
    csp_plp_t *plp;
    int i;
    for (i = 1; i <= mplp->nsg; i++) {
        plp = kh_val(mplp->hsg, mplp->hsg_iter[i]);
        if (plp->ad) ksprintf(ks_ad, "%d\t%ld\n", i, plp->ad);
        if (plp->dp) ksprintf(ks_dp, "%d\t%ld\n", i, plp->dp);
        if (plp->oth) ksprintf(ks_oth, "%d\t%ld\n", i, plp->oth);        
    }
    kputc('\n', ks_ad); kputc('\n', ks_dp); kputc('\n', ks_oth);
    return 0; 
}

int csp_mplp_to_mtx(csp_mplp_t *mplp, jfile_t *fs_ad, jfile_t *fs_dp, jfile_t *fs_oth, size_t idx) {
    csp_plp_t *plp;
    int i;
    for (i = 1; i <= mplp->nsg; i++) {
        plp = kh_val(mplp->hsg, mplp->hsg_iter[i - 1]);
        if (plp->ad)
            fs_ad->is_tmp ? jf_printf(fs_ad, "%d\t%ld\n", i, plp->ad) : jf_printf(fs_ad, "%ld\t%d\t%ld\n", idx, i, plp->ad);
        if (plp->dp)
            fs_dp->is_tmp ? jf_printf(fs_dp, "%d\t%ld\n", i, plp->dp) : jf_printf(fs_dp, "%ld\t%d\t%ld\n", idx, i, plp->dp);
        if (plp->oth)
            fs_oth->is_tmp ? jf_printf(fs_oth, "%d\t%ld\n", i, plp->oth) : jf_printf(fs_oth, "%ld\t%d\t%ld\n", idx, i, plp->oth);      
    }
    if (fs_ad->is_tmp) jf_putc('\n', fs_ad);
    if (fs_dp->is_tmp) jf_putc('\n', fs_dp);
    if (fs_oth->is_tmp) jf_putc('\n', fs_oth);
    return 0; 
}