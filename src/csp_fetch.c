/* cellsnp fetch method
 * Author: Xianejie Huang <hxj5@hku.hk>
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "thpool.h"
#include "htslib/sam.h"
#include "config.h"
#include "csp.h"
#include "jfile.h"
#include "jsam.h"
#include "jnumeric.h"
#include "jstring.h"
#include "mplp.h"
#include "snp.h"

#if CSP_FIT_MULTI_SMP
    #include <errno.h>
#endif

/*@abstract  Pileup one read obtained by sam_itr_next().
@param pos   Pos of the reference sequence. 0-based.
@param p     Pointer of csp_pileup_t structure coming from csp_pileup_init() or csp_pileup_reset().
@param gs    Pointer of global settings.
@return      0 if success, -1 if error, 1 if the reads extracted are not in proper format, 2 if not passing filters.

@note        1. This function is modified from cigar_resolve2() function in sam.c of htslib.
             2. Reads filtering is also applied inside this function, including:
                   UMI and cell tags, read mapping quality, mapping flag and length of bases within alignment.
             3. To speed up, parameters will not be checked, so the caller should guarantee the parameters are valid, i.e.
                && p != NULL && gs != NULL.

@TODO        1. Filter unmapped reads (the read itself unmapped or the mate read unmapped) ?
             2. All elements of csp_pileup_t other than 'base' and 'qual' would be overwritten in this function, while 'base' and
                'qual' would not be overwritten when the next op after finding the pos is BAM_CDEL or BAM_CREF_SKIP, in which case
                the read would be filtered as being DEL in current version. So 'base' and 'qual' would not be misused for the moment.
                But it's better to call csp_pileup_reset*() in case that we donot filter DEL.
 */
static int fetch_read(hts_pos_t pos, csp_pileup_t *p, global_settings *gs) {
    /* Filter reads in order. For example, filtering according to umi tag and cell tag would speed up in the case
       that do not use UMI or Cell-barcode at all. */
    if (use_umi(gs) && NULL == (p->umi = get_bam_aux_str(p->b, gs->umi_tag))) { return 1; }
    if (use_barcodes(gs) && NULL == (p->cb = get_bam_aux_str(p->b, gs->cell_tag))) { return 1; }
    bam1_core_t *c = &(p->b->core);
    if (c->tid < 0 || c->flag & BAM_FUNMAP) { return 2; }
    if (c->qual < gs->min_mapq) { return 2; }
    //if (c->flag > gs->max_flag) { return 2; }
    if (gs->rflag_filter && gs->rflag_filter & c->flag) { return 2; }
    if (gs->rflag_require && ! (gs->rflag_require & c->flag)) { return 2; }
    if (gs->no_orphan && c->flag & BAM_FPAIRED && ! (c->flag & BAM_FPROPER_PAIR)) { return 2; }
    uint32_t *cigar = bam_get_cigar(p->b);
    hts_pos_t x, px;       /* x is the coordinate of the reference. */
    int k, y, py, op, l;   /* y is the query coordinate. */
    uint32_t laln;
    assert(c->pos <= pos);   // otherwise a bug.
    /* find the pos. */
    p->qpos = 0; p->is_refskip = p->is_del = 0;
    for (k = 0, px = x = c->pos, py = y = 0, laln = 0; k < c->n_cigar; k++, px = x, py = y) {
        op = get_cigar_op(cigar[k]);
        l = get_cigar_len(cigar[k]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) { x += l; y += l; laln += l; }
        else if (op == BAM_CDEL || op == BAM_CREF_SKIP) { x += l; }
        else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) { y += l; }
        // else, do nothing.
        if (x > pos) { break; }
    }
    /* pileup */
    assert(k < c->n_cigar);   // otherwise a bug.
    if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
        p->qpos = py + (pos - px); 
        p->base = bam_seqi(bam_get_seq(p->b), p->qpos);
        p->qual = bam_get_qual(p->b)[p->qpos];
    } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
        p->is_del = 1; p->qpos = py; // FIXME: distinguish D and N!!!!!
        p->is_refskip = (op == BAM_CREF_SKIP);
    } // cannot be other operations; otherwise a bug
    if (p->is_del) { return 2; }
    if (p->is_refskip) { return 2; }
    /* continue processing cigar string. */
    for (k++; k < c->n_cigar; k++) {
        op = get_cigar_op(cigar[k]);
        l = get_cigar_len(cigar[k]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) { laln += l; }
    }
    if (laln < gs->min_len) { return 2; }
    else { p->laln = laln; }
    return 0;
}

/*@abstract    Pileup one SNP with method fetch.
@param snp     Pointer of snp_t structure.
@param fs      Pointer of array of pointers to the csp_bam_fs structures.
@param fp      Pointer of array of htsFile* of input files.
@param nfs     Size of @p fs.
@param pileup  Pointer of csp_pileup_t structure.
@param mplp    Pointer of csp_mplp_t structure.
@param gs      Pointer of global_settings structure.
@return        0 if success, -1 if error, 1 if pileup failure without error.

@note          1. This function is mainly called by csp_fetch_core(). Refer to csp_fetch_core() for notes.
               2. The statistics results of all pileuped reads for one SNP is stored in the csp_mplp_t after calling this function.
*/
static int fetch_snp(snp_t *snp, csp_bam_fs **fs, htsFile **fp, int nfs, csp_pileup_t *pileup, csp_mplp_t *mplp, global_settings *gs) 
{
    csp_bam_fs *bs = NULL;
    hts_itr_t *iter = NULL;
    int i, tid, r, ret, st, state = -1;
    size_t npushed = 0;
    kstring_t ks = KS_INITIALIZE, *s = &ks;
  #if DEBUG
    size_t npileup = 0;
  #endif
    mplp->ref_idx = snp->ref ? seq_nt16_char2int(snp->ref) : -1;
    mplp->alt_idx = snp->alt ? seq_nt16_char2int(snp->alt) : -1;
    for (i = 0; i < nfs; i++) {
        bs = fs[i];
        tid = csp_sam_hdr_name2id(bs->hdr, snp->chr, s);
        ks_clear(s);
        if (tid < 0) { state = 1; goto fail; }
        if (NULL == (iter = sam_itr_queryi(bs->idx, tid, snp->pos, snp->pos + 1))) { state = 1; goto fail; }
        while ((ret = sam_itr_next(fp[i], iter, pileup->b)) >= 0) {   // TODO: check if need to be reset in_fp?
          #if DEBUG
            npileup++;
          #endif
            if (0 == (st = fetch_read(snp->pos, pileup, gs))) { // no need to reset pileup as the values in it will be immediately overwritten.
                if (use_barcodes(gs)) { r = csp_mplp_push(pileup, mplp, -1, gs); }
                else if (use_sid(gs)) { r = csp_mplp_push(pileup, mplp, i, gs); }
                else { state = -1; goto fail; }
                if (r < 0) { state = -1; goto fail; }  // else if r == 1: pileuped barcode is not in the input barcode list.
                else if (r == 0) { npushed++; }
            } else if (st < 0) { state = -1; goto fail; }
        }
        if (ret < -1) { state = -1; goto fail; } 
        else { hts_itr_destroy(iter); iter = NULL; }  // TODO: check if could reset iter?
    }
  #if DEBUG
    fprintf(stderr, "[D::%s] before mplp statistics: npileup = %ld; npushed = %ld; the mplp is:\n", __func__, npileup, npushed);
    csp_mplp_print_(stderr, mplp, "\t");
  #endif
    if (npushed < gs->min_count) { state = 1; goto fail; }
    if ((ret = csp_mplp_stat(mplp, gs)) != 0) { state = (ret > 0) ? 1 : -1; goto fail; }
  #if DEBUG
    fprintf(stderr, "[D::%s] after mplp statistics: the mplp is:\n", __func__);
    csp_mplp_print_(stderr, mplp, "\t");
  #endif
    ks_free(s); s = NULL;
    return 0;
  fail:
    if (s) { ks_free(s); }
    if (iter) { hts_itr_destroy(iter); }
    return state;
}

/*@abstract  Pileup a region (a list of SNPs) with method of fetching.
@param args  Pointer to thread_data structure.
@return      Num of SNPs, including those filtered, that are processed.

@note        1. Differ from pileup method in samtools, this function fetches reads covering the SNPs and 
                pileups the reads by processing CIGAR strings with a self-defined resolver function.
             2. The internal variable "ret" in thread_data structure saves the running state of the function:
                  0 if success, 1 if error of other threads.
                  -1, undefined errno in this thread.
                  -2, open error in this thread.
             3. This function could be used by Mode1 and Mode3.		 
 */
static size_t csp_fetch_core(void *args) {
    thread_data *d = (thread_data*) args;
    global_settings *gs = d->gs;
    snp_t **a = gs->pl.a + d->n;  /* here we use directly the internal variables in snplist_t structure to speed up. */
    size_t n = 0;             /* n is the num of SNPs that are successfully processed. */
    csp_bam_fs **bam_fs = d->bfs;
    int nfs = d->nfs;
    htsFile **fp = NULL;
    int nfp = 0;
    csp_pileup_t *pileup = NULL;
    csp_mplp_t *mplp = NULL;
    int i, ret;
    kstring_t ks = KS_INITIALIZE, *s = &ks;
  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto fail; }
  #endif
  #if DEBUG
    fprintf(stderr, "[D::%s][Thread-%d] thread options:\n", __func__, d->i);
    thdata_print(stderr, d);
  #endif
    d->ret = -1;
    d->ns = d->nr_ad = d->nr_dp = d->nr_oth = 0;
    /* prepare data and structures. 
    */
  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto fail; }
  #endif
    if (jf_open(d->out_mtx_ad, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp mtx AD file '%s'.\n", __func__, d->out_mtx_ad->fn);
        d->ret = -2; goto fail;
    }
    if (jf_open(d->out_mtx_dp, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp mtx DP file '%s'.\n", __func__, d->out_mtx_dp->fn);
        d->ret = -2; goto fail;
    }
    if (jf_open(d->out_mtx_oth, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp mtx OTH file '%s'.\n", __func__, d->out_mtx_oth->fn);
        d->ret = -2; goto fail;
    }
    if (jf_open(d->out_vcf_base, NULL) <= 0) { 
        fprintf(stderr, "[E::%s] failed to open tmp vcf BASE file '%s'.\n", __func__, d->out_vcf_base->fn);
        d->ret = -2; goto fail;
    }
    if (gs->is_genotype) {
        if (jf_open(d->out_vcf_cells, NULL) <= 0) { 
            fprintf(stderr, "[E::%s] failed to open tmp vcf CELLS file '%s'.\n", __func__, d->out_vcf_cells->fn);
            d->ret = -2; goto fail;
        }
    }
    /* open input files */ 
  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto fail; }
  #endif
    fp = (htsFile**) calloc(gs->nin, sizeof(htsFile*));
    if (NULL == fp) { fprintf(stderr, "[E::%s] failed to open input files\n", __func__); goto fail; }                 
    for (; nfp < gs->nin; ) {
        if (d->i == 0) {               // the caller has opened input files for Thread-0 
            fp[nfp] = bam_fs[nfp]->fp; nfp++;
        } else if (NULL == (fp[nfp] = hts_open(gs->in_fns[nfp], "rb"))) {
            fprintf(stderr, "[E::%s] failed to open %s.\n", __func__, gs->in_fns[nfp]);
            d->ret = -2; goto fail;
        } else { nfp++; }
    }
    /* prepare mplp for pileup. */
  #if CSP_FIT_MULTI_SMP
    if (gs->tp_errno) { d->ret = 1; goto fail; }
  #endif
    if (NULL == (mplp = csp_mplp_init())) { fprintf(stderr, "[E::%s] could not init csp_mplp_t structure.\n", __func__); goto fail; }
    if (csp_mplp_prepare(mplp, gs) < 0) { fprintf(stderr, "[E::%s] could not prepare csp_mplp_t structure.\n", __func__); goto fail; }
    if (NULL == (pileup = csp_pileup_init())) { 
        fprintf(stderr, "[E::%s] Out of memory allocating csp_pileup_t struct.\n", __func__); 
        goto fail; 
    }
  #if VERBOSE
    double pos_m, pos_n, pos_r, nprints = 50;
    pos_n = pos_m = d->m / nprints;
    pos_r = 100.0 / d->m;
  #endif
    /* pileup each SNP. 
    */
    for (; n < d->m; n++) {
      #if CSP_FIT_MULTI_SMP
        if (gs->tp_errno) { d->ret = 1; goto fail; }
      #endif
      #if VERBOSE
        if (n >= pos_n) {
            fprintf(stderr, "[I::%s][Thread-%d] %.2f%% SNPs processed.\n", __func__, d->i, n * pos_r);
            pos_n += pos_m;
            pos_n = pos_n <= d->m ? pos_n : d->m;
        }
      #endif
      #if DEBUG
        fputc('\n', stderr);
        fprintf(stderr, "[D::%s] chr = %s; pos = %ld; ref = %c; alt = %c;\n", __func__, a[n]->chr, a[n]->pos + 1, a[n]->ref, a[n]->alt);
      #endif
        if ((ret = fetch_snp(a[n], bam_fs, fp, nfs, pileup, mplp, gs)) != 0) {
            if (ret < 0) {
                fprintf(stderr, "[E::%s] failed to pileup snp (%s:%ld)\n", __func__, a[n]->chr, a[n]->pos + 1);
                goto fail; 
            }
          #if DEBUG
            fprintf(stderr, "[W::%s] snp (%s:%ld) filtered, error code = %d\n", __func__, a[n]->chr, a[n]->pos + 1, ret);
          #endif
            csp_mplp_reset(mplp); ks_clear(s);
            continue;
        } else { d->ns++; }
        d->nr_ad += mplp->nr_ad; d->nr_dp += mplp->nr_dp; d->nr_oth += mplp->nr_oth;
        /* output mplp to mtx and vcf. */
        csp_mplp_to_mtx(mplp, d->out_mtx_ad, d->out_mtx_dp, d->out_mtx_oth, d->ns);
        ksprintf(s, "%s\t%ld\t.\t%c\t%c\t.\tPASS\tAD=%ld;DP=%ld;OTH=%ld", a[n]->chr, a[n]->pos + 1, \
                seq_nt16_int2char(mplp->ref_idx), seq_nt16_int2char(mplp->alt_idx), mplp->ad, mplp->dp, mplp->oth);
        jf_puts(ks_str(s), d->out_vcf_base); jf_putc('\n', d->out_vcf_base);
        if (gs->is_genotype) {
            jf_puts(ks_str(s), d->out_vcf_cells);
            jf_puts("\tGT:AD:DP:OTH:PL:ALL", d->out_vcf_cells);
            csp_mplp_to_vcf(mplp, d->out_vcf_cells);
            jf_putc('\n', d->out_vcf_cells);
        }
        csp_mplp_reset(mplp); ks_clear(s);
    }
    // clean
    ks_free(s); s = NULL;
    jf_close(d->out_mtx_ad); jf_close(d->out_mtx_dp); jf_close(d->out_mtx_oth);
    jf_close(d->out_vcf_base); if (gs->is_genotype) { jf_close(d->out_vcf_cells); }
    if (d->i > 0) {
        for (i = 0; i < nfp; i++) { hts_close(fp[i]); }
    } free(fp); fp = NULL;
    csp_pileup_destroy(pileup);
    csp_mplp_destroy(mplp);
    d->ret = 0;
    return n;
  fail:
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
    if (fp) {
        if (d->i > 0) {
            for (i = 0; i < nfp; i++) { hts_close(fp[i]); }
        } free(fp);
    }
    if (pileup) csp_pileup_destroy(pileup);
    if (mplp) { csp_mplp_destroy(mplp); }
    return n;
}

#if CSP_FIT_MULTI_SMP
/*@abstract  Infer appropriate number of threads for multi samples to
             avoid the issue that too many open files
@param gs    Pointer to global_settings
@return      Infered nthread value, no more than 1.
 */
static inline int infer_nthread(global_settings *gs) {
    if (gs->tp_ntry == 0) { return gs->mthread; }  // the first time to try, just use the value user specified
    if (gs->tp_ntry == 1) { 
        int n0 = 3;      // FIXME!!! the initial files opened by this program. eg. the program itself and lz, lhts etc.
        int n = gs->nin + 4 + gs->is_genotype;    // 4 is the 4 output files: base.vcf, ad.mtx, dp.mtx and oth.mtx
        int m = (gs->tp_max_open - n0) / n;
        if (m >= gs->mthread) { m = gs->mthread - 1; }
        if (m < 1) { m = 1; }
        return m;
    } else { return gs->nthread > 1 ? gs->nthread - 1 : 1; }
}
#endif

/*abstract  Run cellSNP Mode with method of fetching.
@param gs   Pointer to the global_settings structure.
@return     0 if success, -1 otherwise.
 */
int csp_fetch(global_settings *gs) {
    /* check options (input) */
    if (NULL == gs || gs->nin <= 0 || (gs->nbarcode <= 0 && gs->nsid <= 0) || snplist_size(gs->pl) <= 0 || \
            NULL == gs->out_dir) {
        fprintf(stderr, "[E::%s] error options for fetch modes.\n", __func__);
        return -1;
    }
    /* prepare running data & options for each thread based on the checked global parameters.*/
    if (gs->nthread > 1 && NULL == gs->tp && NULL == (gs->tp = thpool_init(gs->nthread))) {
        fprintf(stderr, "[E::%s] could not initialize the thread pool.\n", __func__);
        return -1;
    }
    /* core part. */
    int nsample = use_barcodes(gs) ? gs->nbarcode : gs->nsid;
    int nthread = gs->nthread;
    thread_data **td = NULL, *d = NULL;
    int ntd = 0, mtd; // ntd: num of thread-data structures that have been created. mtd: size of td array.
    csp_bam_fs **bam_fs = NULL;       /* use array instead of single element to compatible with multi-input-files. */
    int nfs = 0;
    csp_bam_fs *bs = NULL;
    int i, ret;
    size_t npos, mpos, rpos, tpos, ns, nr_ad, nr_dp, nr_oth, ns_merge, nr_merge;
    jfile_t **out_tmp_mtx_ad, **out_tmp_mtx_dp, **out_tmp_mtx_oth, **out_tmp_vcf_base, **out_tmp_vcf_cells;
    out_tmp_mtx_ad = out_tmp_mtx_dp = out_tmp_mtx_oth = out_tmp_vcf_base = out_tmp_vcf_cells = NULL;
    /* calc number of threads and number of SNPs for each thread. */
    mtd = min2(snplist_size(gs->pl), nthread);
    mpos = snplist_size(gs->pl) / mtd;
    rpos = snplist_size(gs->pl) - mpos * mtd;     // number of remaining positions
    /* create output tmp filenames. */
    if (NULL == (out_tmp_mtx_ad = create_tmp_files(gs->out_mtx_ad, mtd, CSP_TMP_ZIP))) {
        fprintf(stderr, "[E::%s] fail to create tmp files for mtx_AD.\n", __func__);
        goto fail;
    }
    if (NULL == (out_tmp_mtx_dp = create_tmp_files(gs->out_mtx_dp, mtd, CSP_TMP_ZIP))) {
        fprintf(stderr, "[E::%s] fail to create tmp files for mtx_DP.\n", __func__);
        goto fail;
    }
    if (NULL == (out_tmp_mtx_oth = create_tmp_files(gs->out_mtx_oth, mtd, CSP_TMP_ZIP))) {
        fprintf(stderr, "[E::%s] fail to create tmp files for mtx_OTH.\n", __func__);
        goto fail;
    }
    if (mtd > 1) {
        if (NULL == (out_tmp_vcf_base = create_tmp_files(gs->out_vcf_base, mtd, CSP_TMP_ZIP))) {
            fprintf(stderr, "[E::%s] fail to create tmp files for vcf_BASE.\n", __func__);
            goto fail;
        }
        if (gs->is_genotype && NULL == (out_tmp_vcf_cells = create_tmp_files(gs->out_vcf_cells, mtd, CSP_TMP_ZIP))) {
            fprintf(stderr, "[E::%s] fail to create tmp files for vcf_CELLS.\n", __func__);
            goto fail;
        }
    }
    /* construct bam_fs */
    bam_fs = (csp_bam_fs**) calloc(gs->nin, sizeof(csp_bam_fs*));
    if (NULL == bam_fs) { fprintf(stderr, "[E::%s] could not initialize csp_bam_fs array.\n", __func__); goto fail; }
    for (nfs = 0; nfs < gs->nin; nfs++) {
        if (NULL == (bs = csp_bam_fs_init())) { fprintf(stderr, "[E::%s] failed to create csp_bam_fs.\n", __func__); goto fail; }
        if (NULL == (bs->fp = hts_open(gs->in_fns[nfs], "rb"))) {
            fprintf(stderr, "[E::%s] failed to open %s.\n", __func__, gs->in_fns[nfs]); 
            goto fail;
        }
        if (NULL == (bs->hdr = sam_hdr_read(bs->fp))) {
            fprintf(stderr, "[E::%s] failed to read header for %s.\n", __func__, gs->in_fns[nfs]);
            goto fail; 
        }
        if (NULL == (bs->idx = sam_index_load(bs->fp, gs->in_fns[nfs]))) {
            fprintf(stderr, "[E::%s] failed to load index for %s.\n", __func__, gs->in_fns[nfs]);
            goto fail; 
        }
        bam_fs[nfs] = bs;
    } bs = NULL;
    /* prepare data for thread pool. */
    td = (thread_data**) calloc(mtd, sizeof(thread_data*));
    if (NULL == td) { fprintf(stderr, "[E::%s] could not initialize the array of thread_data structure.\n", __func__); goto fail; }
    for (npos = 0; ntd < mtd; ntd++, npos += tpos) {
        if (NULL == (d = thdata_init())) {
            fprintf(stderr, "[E::%s] could not initialize the thread_data structure.\n", __func__); 
            goto fail; 
        }
        tpos = ntd < rpos ? mpos + 1 : mpos;
        d->i = ntd; d->gs = gs; d->bfs = bam_fs; d->nfs = nfs; d->n = npos; d->m = tpos;
        d->out_mtx_ad = out_tmp_mtx_ad[ntd]; d->out_mtx_dp = out_tmp_mtx_dp[ntd]; d->out_mtx_oth = out_tmp_mtx_oth[ntd];
        if (mtd > 1) {
            d->out_vcf_base = out_tmp_vcf_base[ntd]; d->out_vcf_cells = gs->is_genotype ? out_tmp_vcf_cells[ntd] : NULL;
        } else {
            d->out_vcf_base = gs->out_vcf_base; d->out_vcf_cells = gs->is_genotype ? gs->out_vcf_cells : NULL;
        }
        td[ntd] = d;
    } d = NULL;
    // run the threads
    if (mtd > 1) {
        for (i = 0; i < ntd; i++) {
            if (thpool_add_work(gs->tp, (void*) csp_fetch_core, td[i]) < 0) {
                fprintf(stderr, "[E::%s] could not add thread work (No. %d)\n", __func__, i);
                goto fail;
            } 
        }
        thpool_wait(gs->tp);
    } else { csp_fetch_core(td[0]); }
    /* check running status of threads. */
  #if DEBUG
    for (i = 0; i < mtd; i++) { fprintf(stderr, "[D::%s] ret of thread-%d is %d\n", __func__, i, td[i]->ret); }
  #endif
    for (i = 0; i < mtd; i++) { if (td[i]->ret < 0) goto fail; }
    /* merge tmp files. */
    ns = nr_ad = nr_dp = nr_oth = 0;
    for (i = 0; i < mtd; i++) {
        nr_ad += td[i]->nr_ad; nr_dp += td[i]->nr_dp; nr_oth += td[i]->nr_oth;
        ns += td[i]->ns;
    }
    if (jf_open(gs->out_mtx_ad, NULL) < 0) { fprintf(stderr, "[E::%s] failed to open mtx AD.\n", __func__); goto fail; }
    jf_printf(gs->out_mtx_ad, "%ld\t%d\t%ld\n", ns, nsample, nr_ad);
    merge_mtx(gs->out_mtx_ad, out_tmp_mtx_ad, mtd, &ns_merge, &nr_merge, &ret);
    if (ret < 0 || ns_merge != ns || nr_merge != nr_ad) { fprintf(stderr, "[E::%s] failed to merge mtx AD.\n", __func__); goto fail; }
    jf_close(gs->out_mtx_ad);

    if (jf_open(gs->out_mtx_dp, NULL) < 0) { fprintf(stderr, "[E::%s] failed to open mtx DP.\n", __func__); goto fail; }
    jf_printf(gs->out_mtx_dp, "%ld\t%d\t%ld\n", ns, nsample, nr_dp);
    merge_mtx(gs->out_mtx_dp, out_tmp_mtx_dp, mtd, &ns_merge, &nr_merge, &ret);
    if (ret < 0 || ns_merge != ns || nr_merge != nr_dp) { fprintf(stderr, "[E::%s] failed to merge mtx DP.\n", __func__); goto fail; }
    jf_close(gs->out_mtx_dp);

    if (jf_open(gs->out_mtx_oth, NULL) < 0) { fprintf(stderr, "[E::%s] failed to open mtx OTH.\n", __func__); goto fail; }
    jf_printf(gs->out_mtx_oth, "%ld\t%d\t%ld\n", ns, nsample, nr_oth);
    merge_mtx(gs->out_mtx_oth, out_tmp_mtx_oth, mtd, &ns_merge, &nr_merge, &ret);
    if (ret < 0 || ns_merge != ns || nr_merge != nr_oth) { fprintf(stderr, "[E::%s] failed to merge mtx OTH.\n", __func__); goto fail; }
    jf_close(gs->out_mtx_oth);

    if (mtd > 1) {
        if (jf_open(gs->out_vcf_base, NULL) < 0) { fprintf(stderr, "[E::%s] failed to open vcf BASE.\n", __func__); goto fail; }
        merge_vcf(gs->out_vcf_base, out_tmp_vcf_base, mtd, &ret);
        if (ret < 0) { fprintf(stderr, "[E::%s] failed to merge vcf BASE.\n", __func__); goto fail; }
        jf_close(gs->out_vcf_base);

        if (gs->is_genotype) {
            if (jf_open(gs->out_vcf_cells, NULL) < 0) { fprintf(stderr, "[E::%s] failed to open vcf CELLS.\n", __func__); goto fail; }
            merge_vcf(gs->out_vcf_cells, out_tmp_vcf_cells, mtd, &ret);
            if (ret < 0) { fprintf(stderr, "[E::%s] failed to merge vcf CELLS.\n", __func__); goto fail; }    
            jf_close(gs->out_vcf_cells);     
        }
    }
    /* clean */
    for (i = 0; i < mtd; i++) { thdata_destroy(td[i]); }
    free(td); td = NULL;
    for (i = 0; i < nfs; i++) { csp_bam_fs_destroy(bam_fs[i]); }
    free(bam_fs); bam_fs = NULL;
    if (destroy_tmp_files(out_tmp_mtx_ad, mtd) < 0) {
        fprintf(stderr, "[W::%s] failed to remove tmp mtx AD files.\n", __func__);
    } out_tmp_mtx_ad = NULL;
    if (destroy_tmp_files(out_tmp_mtx_dp, mtd) < 0) {
        fprintf(stderr, "[W::%s] failed to remove tmp mtx DP files.\n", __func__);
    } out_tmp_mtx_dp = NULL;
    if (destroy_tmp_files(out_tmp_mtx_oth, mtd) < 0) {
        fprintf(stderr, "[W::%s] failed to remove tmp mtx OTH files.\n", __func__);
    } out_tmp_mtx_oth = NULL;
    if (mtd > 1) {
        if (destroy_tmp_files(out_tmp_vcf_base, mtd) < 0) {
            fprintf(stderr, "[W::%s] failed to remove tmp vcf BASE files.\n", __func__);
        } out_tmp_vcf_base = NULL;
        if (gs->is_genotype) {         
            if (destroy_tmp_files(out_tmp_vcf_cells, mtd) < 0) {
                fprintf(stderr, "[W::%s] failed to remove tmp vcf CELLS files.\n", __func__);
            } out_tmp_vcf_cells = NULL;
        }
    }
    return 0;
  fail:
    if (td) {
        for (i = 0; i < mtd; i++) { thdata_destroy(td[i]); }
        free(td);
    }
    if (d) { thdata_destroy(d); }
    if (bam_fs) {
        for (i = 0; i < nfs; i++) { csp_bam_fs_destroy(bam_fs[i]); }
        free(bam_fs);
    }
    if (bs) { csp_bam_fs_destroy(bs); }
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
        gs->tp_ntry++;
        gs->nthread = infer_nthread(gs);
        gs->tp_errno = 0;
        fprintf(stderr, "[W::%s] No.%d re-try: set nthread = %d to fix the issue.\n",
                         __func__, gs->tp_ntry, gs->nthread);
        fprintf(stderr, "================================================================================\n");
        return csp_fetch(gs);
    } else { return -1; }
  #else
    return -1;
  #endif
}

