/* csp_pileup.c - cellsnp pileup method
 * Author: Xianejie Huang <hxj5@hku.hk>
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include "thpool.h"
#include "htslib/sam.h"
#include "config.h"
#include "csp.h"
#include "jfile.h"
#include "jsam.h"
#include "jstring.h"
#include "mplp.h"
#include "snp.h"

#if CSP_FIT_MULTI_SMP
    #include <errno.h>
#endif

/* auxiliary data used by @func mp_func. */
typedef struct {
    htsFile *fp;
    hts_itr_t *itr;
    const char *chrom;
    global_settings *gs;
} mp_aux_t;

/*@return   Pointer to mp_aux_t structure if success, NULL otherwise. */
static inline mp_aux_t* mp_aux_init(void) {
    return (mp_aux_t*) calloc(1, sizeof(mp_aux_t));
}

/*@note  As all elements are from external sources so no need to be freed/reset. */
static inline void mp_aux_destroy(mp_aux_t *p) {
    if (p) { free(p); }
}

static inline void mp_aux_reset(mp_aux_t *p) { }

/*!@func
@abstract    bam_plp_auto_f used by bam_mplp_init to extract valid reads to be pushed into bam_mpileup stack.
@param data  Pointer to auxiliary data.
@param b     Pointer to bam1_t structure.
@return      0 on success, -1 on end, < -1 on non-recoverable errors. refer to htslib/sam.h @func bam_plp_init.
@note        This function refers to @func mplp_func in bcftools/mpileup.c.   
*/
static int mp_func(void *data, bam1_t *b) {
    int ret;
    mp_aux_t *dat = (mp_aux_t*) data;
    global_settings *gs = dat->gs;
    bam1_core_t *c;
    do {
        if ((ret = sam_itr_next(dat->fp, dat->itr, b)) < 0) 
            break;
        c = &(b->core);
        if (c->tid < 0 || c->flag & BAM_FUNMAP) 
            continue;
        if (c->qual < gs->min_mapq) 
            continue;
        if (gs->rflag_filter && gs->rflag_filter & c->flag ) 
            continue;
        if (gs->rflag_require && ! (gs->rflag_require & c->flag)) 
            continue;
        if (gs->no_orphan && c->flag & BAM_FPAIRED && ! (c->flag & BAM_FPROPER_PAIR)) 
            continue;
        if (use_target(gs)) {
            int beg = c->pos, end = bam_endpos(b) - 1;
            int overlap = regidx_overlap(gs->targets, dat->chrom, beg, end, NULL);
            if (! overlap) 
                continue;
        }
        break;
    } while (1);
    return ret;
}

/*!@func
@abstract    Pileup one read.
@param pos   Pos of the reference sequence. 0-based.
@param bp    Pointer of bam_pileup1_t containing pileup-ed results.
@param p     Pointer of csp_pileup_t structure coming from csp_pileup_init() or csp_pileup_reset().
@param gs    Pointer of global settings.
@return      0 if success, -1 if error, 1 if the reads extracted are not in proper format, 2 if not passing filters.
@note        1. This function is modified from cigar_resolve2() function in sam.c of htslib.
                Link: https://github.com/samtools/htslib/blob/develop/sam.c#L4526
             2. Reads filtering is also applied inside this function, including:
                   UMI and cell tags and length of bases within alignment.
 */
static int pileup_read(hts_pos_t pos, const bam_pileup1_t *bp, csp_pileup_t *p, global_settings *gs) {
    /* Filter reads in order. For example, filtering according to umi tag and cell tag would speed up in the case
       that do not use UMI or Cell-barcode at all. */
    p->b = bp->b;
    if (use_umi(gs) && NULL == (p->umi = get_bam_aux_str(p->b, gs->umi_tag))) 
        return 1;
    if (use_barcodes(gs) && NULL == (p->cb = get_bam_aux_str(p->b, gs->cell_tag))) 
        return 1;

    bam1_core_t *c = &(p->b->core);
    uint32_t *cigar = bam_get_cigar(p->b);
    int k, op, l;  
    uint32_t laln;
    assert(c->pos <= pos);   // otherwise a bug.
    if (bp->is_del) 
        return 2;
    if (bp->is_refskip) 
        return 2;

    /* processing cigar string to get number of mapped positions. */
    if (gs->min_len > 0) {
        for (k = 0, laln = 0; k < c->n_cigar; k++) {
            op = get_cigar_op(cigar[k]);
            l = get_cigar_len(cigar[k]);
            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
                laln += l;
        }
        if (laln < gs->min_len) 
            return 2;
        else 
            p->laln = laln;
    }
    p->qpos = bp->qpos; 
    p->is_del = bp->is_del; 
    p->is_refskip = bp->is_refskip;
    if (p->qpos < c->l_qseq) { 
        p->base = bam_seqi(bam_get_seq(p->b), p->qpos); 
        p->qual = bam_get_qual(p->b)[p->qpos]; 
    } else {
        p->base = seq_nt16_char2idx('N');
        p->qual = 0;
    }
    return 0;
}

/*!@func
@abstract      Pileup One SNP.
@param pos     Pos of pileup-ed snp.
@param mp_n    Pointer of array containing numbers of bam_pileup1_t* pileup-ed from each file.
@param mp_plp  Pointer of array containing bam_pileup1_t* pileup-ed from each file.
@param nfs     Size of @p mp_nplp and @p mp_plp.
@param pileup  Pointer of csp_pileup_t structure.
@param mplp    Pointer of csp_mplp_t structure.
@param gs      Pointer of global_settings structure.
@return        0 if success, -1 if error, 1 if pileup failure without error.
*/
static int pileup_snp(hts_pos_t pos, int *mp_n, const bam_pileup1_t **mp_plp, int nfs, csp_pileup_t *pileup, csp_mplp_t *mplp, global_settings *gs) 
{
    const bam_pileup1_t *bp = NULL;
    int i, j, r, ret, st, state = -1;
    size_t npush1;

  #if DEBUG
    size_t npileup = 0;
  #endif

    for (i = 0, npush1 = 0; i < nfs; i++, npush1 = 0) {
        for (j = 0; j < mp_n[i]; j++) {
            bp = mp_plp[i] + j;
          #if DEBUG
            npileup++;
          #endif
            if (npush1 >= gs->max_depth) {
                npush1 = 0;
                break;
            }
            if (0 == (st = pileup_read(pos, bp, pileup, gs))) { // no need to reset pileup as the values in it will be immediately overwritten.
                if (use_barcodes(gs)) {
                    r = csp_mplp_push(pileup, mplp, -1, gs);
                } else if (use_sid(gs)) {
                    r = csp_mplp_push(pileup, mplp, i, gs);
                } else {
                    state = -1; goto fail;
                }
                if (r < 0) {
                    state = -1; goto fail;
                } else if (r == 0) {
                    npush1++;
                } // else if r > 0: do nothing
            } else if (st < 0) {
                state = -1; goto fail;
            }
        }
    }

  #if DEBUG
    fprintf(stderr, "[D::%s] before mplp statistics: npileup = %ld; the mplp is:\n", __func__, npileup);
    csp_mplp_print_(stderr, mplp, "\t");
  #endif

    mplp->pos = pos;
    if ((ret = csp_mplp_stat(mplp, gs)) != 0) {
        state = (ret > 0) ? 1 : -1;
        goto fail;
    }

  #if DEBUG
    fprintf(stderr, "[D::%s] after mplp statistics: the mplp is:\n", __func__);
    csp_mplp_print_(stderr, mplp, "\t");
  #endif

    return 0;

  fail:
    return state;
}

/*!@func
@abstract    Pileup regions (several chromosomes).
@param args  Pointer to thread_data structure.
@return      Num of SNPs, including those filtered, that are processed.
@note  1. The internal variable "ret" in thread_data structure saves the running state of the function:
         0 if success, 1 if error of other threads.
         -1, undefined errno in this thread.
         -2, open error in this thread.
       2. This function could be used by Mode2 & Mode1a(-T) & Mode1b(-T).
       3. Refer to @func mpileup from bam_plcmd.c in @repo samtools, the @p mp_iter, @p mp_plp and
          @p mp_n do not need to be reset every time calling bam_mplp_auto(). Guess that there may be 
          memory pools inside bam_mplp_* for these structures.
       4. TODO: check consistency among headers. For now not sure the pileup-ed @p tid from @func 
            bam_mplp_auto() is corresponded to which file's header for different files may have different 
            target_names in their headers.
       5. TODO: Refer to @func mpileup from bam_plcmd.c in @repo samtools, the hts_itr_t* structure is
            reused directly without destroying-creating again. is it a good way? it may speed up if 
            donot repeat the create-destroy-create-... process.
 */
static int csp_pileup_core(void *args) {
    thread_data *d = (thread_data*) args;
    global_settings *gs = d->gs;
    char **a = gs->chroms + d->n;
    int n = 0;                   /* n is the num of chroms that are successfully processed. */

    csp_bam_fs **bam_fs = d->bfs;
    int nfs = d->nfs;
    htsFile **fp = NULL;
    int nfp = 0;

    csp_pileup_t *pileup = NULL;
    csp_mplp_t *mplp = NULL;
    bam_mplp_t mp_iter = NULL;
    const bam_pileup1_t **mp_plp = NULL;
    int *mp_n = NULL;

    mp_aux_t **data = NULL;
    int ndat = 0;                 // num of elements in array of mp_aux_t data.

    int tid;
    int pos;
    int i, r, ret;
    size_t msnp, nsnp, unit = 200000;
    regitr_t *itr = NULL;
    kstring_t ks = KS_INITIALIZE, *s = &ks;

  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto clean; }
  #endif

  #if DEBUG
    fprintf(stderr, "[D::%s][Thread-%d] thread options:\n", __func__, d->i);
    thdata_print(stderr, d);
  #endif

    assert(d->nfs == gs->nin);
    assert(d->niter == d->m);
    assert(d->nitr == gs->nin);
    d->ret = -1;
    d->ns = d->nr_ad = d->nr_dp = d->nr_oth = 0;

    /* prepare data and structures. 
    */
  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto clean; }
  #endif

    if (jf_open(d->out_mtx_ad, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp mtx AD file '%s'.\n", __func__, d->out_mtx_ad->fn);
        d->ret = -2; goto clean;
    }
    if (jf_open(d->out_mtx_dp, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp mtx DP file '%s'.\n", __func__, d->out_mtx_dp->fn);
        d->ret = -2; goto clean;
    }
    if (jf_open(d->out_mtx_oth, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp mtx OTH file '%s'.\n", __func__, d->out_mtx_oth->fn);
        d->ret = -2; goto clean;
    }
    if (jf_open(d->out_vcf_base, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp vcf BASE file '%s'.\n", __func__, d->out_vcf_base->fn);
        d->ret = -2; goto clean;
    }
    if (gs->is_genotype) {
        if (jf_open(d->out_vcf_cells, NULL) <= 0) { 
            fprintf(stderr, "[E::%s] failed to open tmp vcf CELLS file '%s'.\n", __func__, d->out_vcf_cells->fn);
            d->ret = -2; goto clean;
        }
    }

    /* open input files */ 
  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto clean; }
  #endif

    fp = (htsFile**) calloc(gs->nin, sizeof(htsFile*));
    if (NULL == fp) {
        fprintf(stderr, "[E::%s] failed to open input files\n", __func__);
        goto clean;
    }
    for (; nfp < gs->nin; ) {
        if (d->i == 0) {               // the caller has opened input files for Thread-0 
            fp[nfp] = bam_fs[nfp]->fp;
            nfp++;
        } else if (NULL == (fp[nfp] = hts_open(gs->in_fns[nfp], "rb"))) {
            fprintf(stderr, "[E::%s] failed to open %s.\n", __func__, gs->in_fns[nfp]);
            d->ret = -2; goto clean;
        } else {
            nfp++;
        }
    }

    /* prepare mplp for pileup. */
  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto clean; }
  #endif

    if (NULL == (mplp = csp_mplp_init())) {
        fprintf(stderr, "[E::%s] could not init csp_mplp_t structure.\n", __func__);
        goto clean;
    }
    if (csp_mplp_prepare(mplp, gs) < 0) {
        fprintf(stderr, "[E::%s] could not prepare csp_mplp_t structure.\n", __func__);
        goto clean;
    }
    if (NULL == (pileup = csp_pileup_init())) { 
        fprintf(stderr, "[E::%s] Out of memory allocating csp_pileup_t struct.\n", __func__); 
        goto clean; 
    }

    /* create bam_mplp_* & mp_aux_t data structures */
    if (NULL == (data = (mp_aux_t**) malloc(nfs * sizeof(mp_aux_t*)))) {
        fprintf(stderr, "[E::%s] failed to allocate space for mp_aux_t data.\n", __func__);
        goto clean;
    }
    for (; ndat < nfs; ndat++) {
        if (NULL == (data[ndat] = mp_aux_init())) {
            fprintf(stderr, "[E::%s] failed to allocate space for mp_aux_t.\n", __func__);
            goto clean;
        } else {
            data[ndat]->fp = fp[ndat];
            data[ndat]->gs = gs;
        }
    }
    if (NULL == (mp_plp = (const bam_pileup1_t**) calloc(nfs, sizeof(bam_pileup1_t*)))) {
        fprintf(stderr, "[E::%s] failed to allocate space for mp_plp.\n", __func__);
        goto clean;
    }
    if (NULL == (mp_n = (int*) calloc(nfs, sizeof(int)))) {
        fprintf(stderr, "[E::%s] failed to allocate space for mp_nplp.\n", __func__);
        goto clean;
    }

    // init reg itr
    if (use_target(gs) && NULL == (itr = regitr_init(gs->targets))) {
        fprintf(stderr, "[E::%s] failed to init regitr.\n", __func__);
        goto clean;
    }

    /* pileup each SNP. 
    */
    // init mpileup 
    if (gs->max_depth > (1 << 20) / (float) nfs) {
        fprintf(stderr, "[W::%s] Combined max depth is above 1M. Potential memory hog!\n", __func__);
    }
    for (msnp = nsnp = 0; n < d->m; n++, msnp = nsnp = 0) {
      #if CSP_FIT_MULTI_SMP
        if (gs->tp_errno) { d->ret = 1; goto clean; }
      #endif

      #if VERBOSE
        fprintf(stderr, "[I::%s][Thread-%d] processing chrom %s ...\n", __func__, d->i, a[n]);
      #endif

        /* create bam_mplp_* mpileup structure from htslib */
        for (i = 0; i < ndat; i++) {
            data[i]->itr = d->iter[n][i];
            data[i]->chrom = a[n];
        }
        if (NULL == (mp_iter = bam_mplp_init(nfs, mp_func, (void**) data))) {
            fprintf(stderr, "[E::%s] failed to create mp_iter for chrom %s.\n", __func__, a[n]);
            goto clean;
        }
        bam_mplp_set_maxcnt(mp_iter, gs->max_depth);

        // As each query region is a chrom, so no need to call bam_mplp_init_overlaps() here?
        /* begin mpileup */
        mplp->chrom = a[n];
        while ((ret = bam_mplp_auto(mp_iter, &tid, &pos, mp_n, mp_plp)) > 0) {
          #if CSP_FIT_MULTI_SMP
            if (gs->tp_errno) { d->ret = 1; goto clean; }
          #endif
            if (tid < 0) 
                break;
            if (use_target(gs)) {
                int overlap = regidx_overlap(gs->targets, a[n], pos, pos, itr);
                biallele_t *ale = NULL;
                if (! overlap)  // no need to reset mplp_t here 
                    continue; 
                while (regitr_overlap(itr)) {
                    ale = regitr_payload(itr, biallele_t*);
                    break;
                }
                mplp->ref_idx = csp_mplp_base2int(ale->ref);
                mplp->alt_idx = csp_mplp_base2int(ale->alt);
            } else {
                mplp->ref_idx = -1;
                mplp->alt_idx = -1;
            }
            if ((r = pileup_snp(pos, mp_n, mp_plp, nfs, pileup, mplp, gs)) != 0) {
                if (r < 0) {
                    fprintf(stderr, "[E::%s] failed to pileup snp for %s:%d\n", __func__, a[n], pos);
                    goto clean; 
                } else {
                    csp_mplp_reset(mplp);
                    continue;
                }
            } else {
                d->ns++;
            }
            d->nr_ad += mplp->nr_ad;
            d->nr_dp += mplp->nr_dp;
            d->nr_oth += mplp->nr_oth;

            /* output mplp to mtx and vcf. */
            csp_mplp_to_mtx(mplp, d->out_mtx_ad, d->out_mtx_dp, d->out_mtx_oth, d->ns);
            ksprintf(s, "%s\t%d\t.\t%c\t%c\t.\tPASS\tAD=%ld;DP=%ld;OTH=%ld", a[n], pos + 1, \
                    seq_nt16_int2char(mplp->ref_idx), seq_nt16_int2char(mplp->alt_idx), mplp->ad, mplp->dp, mplp->oth);
            jf_puts(ks_str(s), d->out_vcf_base); jf_putc('\n', d->out_vcf_base);
            if (gs->is_genotype) {
                jf_puts(ks_str(s), d->out_vcf_cells);
                jf_puts("\tGT:AD:DP:OTH:PL:ALL", d->out_vcf_cells);
                csp_mplp_to_vcf(mplp, d->out_vcf_cells);
                jf_putc('\n', d->out_vcf_cells);
            }
            csp_mplp_reset(mplp); ks_clear(s);

          #if VERBOSE
            if ((++nsnp) - msnp >= unit) {
                fprintf(stderr, "[I::%s][Thread-%d] has pileup-ed %.2fM SNPs for chrom %s\n", __func__, d->i, nsnp / 1000000.0, a[n]);
                msnp = nsnp;
            }
          #endif
        }

        if (ret < 0) {
            fprintf(stderr, "[E::%s] failed to pileup chrom %s\n", __func__, a[n]);
            goto clean;
        }
        for (i = 0; i < ndat; i++) 
            mp_aux_reset(data[i]);

      #if VERBOSE
        fprintf(stderr, "[I::%s][Thread-%d] has pileup-ed in total %ld SNPs for chrom %s\n", __func__, d->i, nsnp, a[n]);
      #endif

    }

    d->ret = 0;

  clean:
  #if CSP_FIT_MULTI_SMP
    if (-2 == d->ret && EMFILE == errno) { gs->tp_errno |= TP_EMFILE; }
    else if (d->ret < 0) { gs->tp_errno |= TP_EUNDEF; }
  #endif
    if (s) { ks_free(s); }
    if (jf_isopen(d->out_mtx_ad)) { jf_close(d->out_mtx_ad); }
    if (jf_isopen(d->out_mtx_dp)) { jf_close(d->out_mtx_dp); }
    if (jf_isopen(d->out_mtx_oth)) { jf_close(d->out_mtx_oth); }
    if (jf_isopen(d->out_vcf_base)) { jf_close(d->out_vcf_base); }
    if (gs->is_genotype && jf_isopen(d->out_vcf_cells)) { jf_close(d->out_vcf_cells); }
    if (data) {
        for (i = 0; i < ndat; i++) { mp_aux_destroy(data[i]); }
        free(data); 
    }
    if (fp) {
        if (d->i > 0) {
            for (i = 0; i < nfp; i++) { hts_close(fp[i]); }
        } free(fp);
    }
    if (mp_plp) free(mp_plp);
    if (mp_n) free(mp_n);
    // do not free mp_iter here, otherwise will lead to double free error!!!
    // seems bam_mplp_* will free the mp_iter by default.
    //bam_mplp_destroy(mp_iter);   
    if (pileup) csp_pileup_destroy(pileup);
    if (mplp) { csp_mplp_destroy(mplp); }
    if (itr) { regitr_destroy(itr); }
    return n;
}

#if CSP_FIT_MULTI_SMP
/*!@func
@abstract  Infer appropriate number of threads for multi samples to
           avoid the issue that too many open files
@param gs  Pointer to global_settings
@return    Infered nthread value, no more than 1.
*/
static inline int infer_nthread(global_settings *gs) {
    if (gs->tp_ntry == 0) { return gs->mthread; }  // the first time to try, just use the value user specified
    if (gs->tp_ntry == 1) { 
        int n0 = 3;      // FIXME!!! the initial files opened by this program. eg. the program itself and lz, lhts etc.
        int n = gs->nin + 4 + gs->is_genotype + (NULL != gs->refseq_file);      // 4 is the 4 output files: base.vcf, ad.mtx, dp.mtx and oth.mtx
        int m = (gs->tp_max_open - n0) / n;
        if (m >= gs->mthread) { m = gs->mthread - 1; }
        if (m < 1) { m = 1; }
        return m;
    } else { return gs->nthread > 1 ? gs->nthread - 1 : 1; }
}
#endif

/*!@func
@abstract  Run cellSNP Mode with method of pileuping.
@param gs  Pointer to the global_settings structure.
@return    0 if success, -1 otherwise.
 */
int csp_pileup(global_settings *gs) {
    /* check options (input) */
    if (NULL == gs || gs->nin <= 0 || gs->nchrom <= 0 || NULL == gs->out_dir) {
        fprintf(stderr, "[E::%s] error options for fetch modes.\n", __func__);
        return -1;
    }

    /* prepare running data & options for each thread based on the checked global parameters.*/ 
    if (gs->nthread > 1 && NULL == gs->tp && NULL == (gs->tp = thpool_init(gs->nthread))) {
        fprintf(stderr, "[E::%s] could not initialize the thread pool.\n", __func__);
        return -1;
    }

    /* core part. */
    int nsample = use_barcodes(gs) ? gs->nbarcode : gs->nin;
    thread_data **td = NULL, *d = NULL;
    int ntd = 0, mtd;            // ntd: num of thread-data structures that have been created. mtd: size of td array.

    csp_bam_fs **bam_fs = NULL;       /* use array instead of single element to compatible with multi-input-files. */
    csp_bam_fs *bs = NULL;
    int nfs = 0;
    hts_itr_t ****titer = NULL;
    hts_itr_t ***iter = NULL;
    hts_itr_t **itr = NULL;
    int ntiter = 0, niter = 0, nitr = 0;

    const char *ref = NULL;
    char **a = NULL;
    kstring_t ks = KS_INITIALIZE, *s = &ks;
    int i, j, k, ret, is_ok = 0;

    size_t ns, nr_ad, nr_dp, nr_oth, ns_merge, nr_merge;
    jfile_t **out_tmp_mtx_ad, **out_tmp_mtx_dp, **out_tmp_mtx_oth, **out_tmp_vcf_base, **out_tmp_vcf_cells;

    /* calc number of threads and number of chroms for each thread. */
    mtd = gs->nthread > 1 ? gs->nchrom : 1;

    /* create output tmp filenames. */
    out_tmp_mtx_ad = out_tmp_mtx_dp = out_tmp_mtx_oth = out_tmp_vcf_base = out_tmp_vcf_cells = NULL;
    if (NULL == (out_tmp_mtx_ad = create_tmp_files(gs->out_mtx_ad, mtd, CSP_TMP_ZIP))) {
        fprintf(stderr, "[E::%s] fail to create tmp files for mtx_AD.\n", __func__);
        goto clean;
    }
    if (NULL == (out_tmp_mtx_dp = create_tmp_files(gs->out_mtx_dp, mtd, CSP_TMP_ZIP))) {
        fprintf(stderr, "[E::%s] fail to create tmp files for mtx_DP.\n", __func__);
        goto clean;
    }
    if (NULL == (out_tmp_mtx_oth = create_tmp_files(gs->out_mtx_oth, mtd, CSP_TMP_ZIP))) {
        fprintf(stderr, "[E::%s] fail to create tmp files for mtx_OTH.\n", __func__);
        goto clean;
    }
    if (mtd > 1) {
        if (NULL == (out_tmp_vcf_base = create_tmp_files(gs->out_vcf_base, mtd, CSP_TMP_ZIP))) {
            fprintf(stderr, "[E::%s] fail to create tmp files for vcf_BASE.\n", __func__);
            goto clean;
        }
        if (gs->is_genotype && NULL == (out_tmp_vcf_cells = create_tmp_files(gs->out_vcf_cells, mtd, CSP_TMP_ZIP))) {
            fprintf(stderr, "[E::%s] fail to create tmp files for vcf_CELLS.\n", __func__);
            goto clean;
        }
    }

    /* create csp_bam_fs structures */
    // open input files and construct hdr for Thread-0 and 
    // other threads would use directly hdr of Thread-0 and by themselves open input files.
    bam_fs = (csp_bam_fs**) calloc(gs->nin, sizeof(csp_bam_fs*));
    if (NULL == bam_fs) {
        fprintf(stderr, "[E::%s] could not initialize csp_bam_fs* array.\n", __func__);
        goto clean;
    }
    for (nfs = 0; nfs < gs->nin; nfs++) {
        if (NULL == (bs = csp_bam_fs_init())) {
            fprintf(stderr, "[E::%s] failed to create csp_bam_fs.\n", __func__);
            goto clean;
        }
        if (NULL == (bs->fp = hts_open(gs->in_fns[nfs], "rb"))) {
            fprintf(stderr, "[E::%s] failed to open %s.\n", __func__, gs->in_fns[nfs]);
            goto clean;
        }
        if (NULL == (bs->hdr = sam_hdr_read(bs->fp))) {
            fprintf(stderr, "[E::%s] failed to read header for %s.\n", __func__, gs->in_fns[nfs]);
            goto clean;
        }
        if (NULL == (bs->idx = sam_index_load(bs->fp, gs->in_fns[nfs]))) {
            fprintf(stderr, "[E::%s] failed to load index for %s.\n", __func__, gs->in_fns[nfs]);
            goto clean;
        }
        bam_fs[nfs] = bs;
    }
    bs = NULL;

    /* prepare hts_itr_t */
    titer = (hts_itr_t****) calloc(mtd, sizeof(hts_itr_t***));
    if (NULL == titer) {
        fprintf(stderr, "[E::%s] could not initialize hts_itr_t*** array.\n", __func__);
        goto clean;
    }

    /* prepare data for thread pool. */
    td = (thread_data**) calloc(mtd, sizeof(thread_data*));
    if (NULL == td) {
        fprintf(stderr, "[E::%s] could not initialize the array of thread_data structure.\n", __func__);
        goto clean;
    }

    /* start mpileup */
    for (; ntd < mtd; ntd++) {
        if (NULL == (d = thdata_init())) {
            fprintf(stderr, "[E::%s] could not initialize the thread_data structure.\n", __func__); 
            goto clean; 
        }
        d->n = ntd; d->m = mtd > 1 ? 1 : gs->nchrom;
        a = gs->chroms + d->n;
        d->i = ntd; d->gs = gs;

        // construct csp_bam_fs
        d->bfs = bam_fs; d->nfs = nfs;

        // construct hts_itr_t
        // we choose to create all hts_itrs here instead of creating several required hts_itrs in each thread
        // because in this way we could destroy all idx & hdr before running each thread to save memory as
        // idx & hdr are not needed for bam_mplp_auto() in each thread.
        // But in csp_fetch, the idx & hdr cannot be removed before running each thread as they are required
        // by sam_itr_queryi() in each thread.
        iter = (hts_itr_t***) calloc(d->m, sizeof(hts_itr_t**));
        if (NULL == iter) {
            fprintf(stderr, "[E::%s] failed to allocate space for hts_itr_t***\n", __func__);
            goto clean;
        }
        for (niter = 0; niter < d->m; niter++) {
            itr = (hts_itr_t**) calloc(gs->nin, sizeof(hts_itr_t*));
            if (NULL == itr) {
                fprintf(stderr, "[E::%s] failed to allocate space for hts_itr_t**\n", __func__);
                goto clean;
            }
            for (nitr = 0; nitr < gs->nin; nitr++) {
                if (NULL == (ref = csp_fmt_chr_name(a[niter], bam_fs[nitr]->hdr, s))) {
                    fprintf(stderr, "[E::%s] could not parse name for chrom %s.\n", __func__, a[niter]);
                    goto clean;
                } else {
                    ks_clear(s);
                }
                if (NULL == (itr[nitr] = sam_itr_querys(bam_fs[nitr]->idx, bam_fs[nitr]->hdr, ref))) {
                    fprintf(stderr, "[E::%s] could not parse region for chrom %s.\n", __func__, a[niter]);
                    goto clean;
                }
            }
            iter[niter] = itr;
        }
        itr = NULL;
        titer[ntiter++] = iter;
        d->iter = iter; d->niter = niter; d->nitr = nitr; iter = NULL;

      #if DEBUG
        assert(niter == d->m);
        assert(nitr == gs->nin);
      #endif

        // construct thdata
        d->out_mtx_ad = out_tmp_mtx_ad[ntd];
        d->out_mtx_dp = out_tmp_mtx_dp[ntd];
        d->out_mtx_oth = out_tmp_mtx_oth[ntd];
        if (mtd > 1) {
            d->out_vcf_base = out_tmp_vcf_base[ntd];
            d->out_vcf_cells = gs->is_genotype ? out_tmp_vcf_cells[ntd] : NULL;
        } else {
            d->out_vcf_base = gs->out_vcf_base;
            d->out_vcf_cells = gs->is_genotype ? gs->out_vcf_cells : NULL;
        }
        td[ntd] = d;
    }
    d = NULL;

    // clean hdr and idx
    for (i = 0; i < nfs; i++) {
        sam_hdr_destroy(bam_fs[i]->hdr); bam_fs[i]->hdr = NULL; 
    }
    for (i = 0; i < nfs; i++) {
        hts_idx_destroy(bam_fs[i]->idx); bam_fs[i]->idx = NULL;
    }

    // run threads
    if (mtd > 1) {
        for (i = 0; i < mtd; i++) {
            if (thpool_add_work(gs->tp, (void*) csp_pileup_core, td[i]) < 0) {
                fprintf(stderr, "[E::%s] could not add thread work (No. %d)\n", __func__, i);
                goto clean;
            }
        }
        thpool_wait(gs->tp);
    } else {
        csp_pileup_core(td[0]);
    }

    /* check running status of threads. */
  #if DEBUG
    for (i = 0; i < mtd; i++)
        fprintf(stderr, "[D::%s] ret of thread-%d is %d\n", __func__, i, td[i]->ret);
  #endif
    for (i = 0; i < mtd; i++) 
        if (td[i]->ret < 0)
            goto clean;

    /* merge tmp files. */
    ns = nr_ad = nr_dp = nr_oth = 0;
    for (i = 0; i < mtd; i++) {
        nr_ad += td[i]->nr_ad;
        nr_dp += td[i]->nr_dp;
        nr_oth += td[i]->nr_oth;
        ns += td[i]->ns;
    }

    if (jf_open(gs->out_mtx_ad, NULL) < 0) {
        fprintf(stderr, "[E::%s] failed to open mtx AD.\n", __func__);
        goto clean;
    }
    jf_printf(gs->out_mtx_ad, "%ld\t%d\t%ld\n", ns, nsample, nr_ad);
    merge_mtx(gs->out_mtx_ad, out_tmp_mtx_ad, mtd, &ns_merge, &nr_merge, &ret);
    if (ret < 0 || ns_merge != ns || nr_merge != nr_ad) {
        fprintf(stderr, "[E::%s] failed to merge mtx AD.\n", __func__);
        goto clean;
    }
    jf_close(gs->out_mtx_ad);

    if (jf_open(gs->out_mtx_dp, NULL) < 0) {
        fprintf(stderr, "[E::%s] failed to open mtx DP.\n", __func__);
        goto clean;
    }
    jf_printf(gs->out_mtx_dp, "%ld\t%d\t%ld\n", ns, nsample, nr_dp);
    merge_mtx(gs->out_mtx_dp, out_tmp_mtx_dp, mtd, &ns_merge, &nr_merge, &ret);
    if (ret < 0 || ns_merge != ns || nr_merge != nr_dp) {
        fprintf(stderr, "[E::%s] failed to merge mtx DP.\n", __func__);
        goto clean;
    }
    jf_close(gs->out_mtx_dp);

    if (jf_open(gs->out_mtx_oth, NULL) < 0) {
        fprintf(stderr, "[E::%s] failed to open mtx OTH.\n", __func__);
        goto clean;
    }
    jf_printf(gs->out_mtx_oth, "%ld\t%d\t%ld\n", ns, nsample, nr_oth);
    merge_mtx(gs->out_mtx_oth, out_tmp_mtx_oth, mtd, &ns_merge, &nr_merge, &ret);
    if (ret < 0 || ns_merge != ns || nr_merge != nr_oth) {
        fprintf(stderr, "[E::%s] failed to merge mtx OTH.\n", __func__);
        goto clean;
    }
    jf_close(gs->out_mtx_oth);

    if (mtd > 1) {
        if (jf_open(gs->out_vcf_base, NULL) < 0) {
            fprintf(stderr, "[E::%s] failed to open vcf BASE.\n", __func__);
            goto clean;
        }
        merge_vcf(gs->out_vcf_base, out_tmp_vcf_base, mtd, &ret);
        if (ret < 0) {
            fprintf(stderr, "[E::%s] failed to merge vcf BASE.\n", __func__);
            goto clean;
        }
        jf_close(gs->out_vcf_base);
        if (gs->is_genotype) {
            if (jf_open(gs->out_vcf_cells, NULL) < 0) {
                fprintf(stderr, "[E::%s] failed to open vcf CELLS.\n", __func__);
                goto clean;
            }
            merge_vcf(gs->out_vcf_cells, out_tmp_vcf_cells, mtd, &ret);
            if (ret < 0) {
                fprintf(stderr, "[E::%s] failed to merge vcf CELLS.\n", __func__);
                goto clean;
            }
            jf_close(gs->out_vcf_cells);     
        }
    }

    is_ok = 1;

    // hdr of other thdata should be set to NULL before being destroyed
    // otherwise will cause double free error!

  clean:
    if (td) {
        for (i = 0; i < mtd; i++) { thdata_destroy(td[i]); }
        free(td);
    }
    if (s) { ks_free(s); }
    if (titer) {
        for (i = 0; i < ntiter; i++) {
            for (j = 0; j < niter; j++) {
                for (k = 0; k < nitr; k++) { hts_itr_destroy(titer[i][j][k]); }
                free(titer[i][j]);
            }
            free(titer[i]);
        }
        free(titer);
    }
    if (bam_fs) {
        for (j = 0; j < nfs; j++) { csp_bam_fs_destroy(bam_fs[j]); }
        free(bam_fs);
    }
    if (out_tmp_mtx_ad && destroy_tmp_files(out_tmp_mtx_ad, mtd) < 0) {
        fprintf(stderr, "[W::%s] failed to remove tmp mtx AD files.\n", __func__);
    }
    if (out_tmp_mtx_dp && destroy_tmp_files(out_tmp_mtx_dp, mtd) < 0) {
        fprintf(stderr, "[W::%s] failed to remove tmp mtx DP files.\n", __func__);
    }
    if (out_tmp_mtx_oth && destroy_tmp_files(out_tmp_mtx_oth, mtd) < 0) {
        fprintf(stderr, "[W::%s] failed to remove tmp mtx OTH files.\n", __func__);
    }
    if (mtd > 1) {
        if (out_tmp_vcf_base && destroy_tmp_files(out_tmp_vcf_base, mtd) < 0) {
            fprintf(stderr, "[W::%s] failed to remove tmp vcf BASE files.\n", __func__);
        }
        if (out_tmp_vcf_cells && destroy_tmp_files(out_tmp_vcf_cells, mtd) < 0) {
            fprintf(stderr, "[W::%s] failed to remove tmp vcf CELLS files.\n", __func__);
        }
    }

    if (is_ok)
        return 0;

  //fail:
    if (d) { thdata_destroy(d); }
    if (iter) {
        for (j = 0; j < niter; j++) {
            for (k = 0; k < nitr; k++) { hts_itr_destroy(iter[j][k]); }
            free(iter[j]);
        }
        free(iter);        
    }
    if (itr) {
        for (k = 0; k < nitr; k++) { hts_itr_destroy(itr[k]); }
        free(itr);
    }
    if (bs) { csp_bam_fs_destroy(bs); }
    if (jf_isopen(gs->out_mtx_ad)) { jf_close(gs->out_mtx_ad); }
    if (jf_isopen(gs->out_mtx_dp)) { jf_close(gs->out_mtx_dp); }
    if (jf_isopen(gs->out_mtx_oth)) { jf_close(gs->out_mtx_oth); }
    if (jf_isopen(gs->out_vcf_base)) { jf_close(gs->out_vcf_base); }
    if (gs->is_genotype && jf_isopen(gs->out_vcf_cells)) { jf_close(gs->out_vcf_cells); }

  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno & TP_EMFILE && ! (gs->tp_errno & TP_EUNDEF) && gs->nthread > 1) {
        fprintf(stderr, "================================================================================\n");
        fprintf(stderr, "[W::%s] Last try (nthreads = %d) failed due to the issue of too many open files.\n",
                         __func__, gs->nthread);
        // reset gs, especially items related to thpool
        if (gs->tp) {
            thpool_destroy(gs->tp); gs->tp = NULL;
        }
        gs->tp_ntry++;
        gs->nthread = infer_nthread(gs);
        gs->tp_errno = 0;
        fprintf(stderr, "[W::%s] No.%d re-try: set nthread = %d to fix the issue.\n",
                         __func__, gs->tp_ntry, gs->nthread);
        fprintf(stderr, "================================================================================\n");
        return csp_pileup(gs);
    } else {
        return -1;
    }
  #else
    return -1;
  #endif
}

