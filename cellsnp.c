
/* TODO: try using multi_iter fetching method of bam/sam/cram for multi regions (SNPs) if it can in theory speed cellsnp up.
 */
#define DEBUG 0
#define VERBOSE 1
/* DEVELOP defined to 1 means some codes for future version of cellsnp will be included. */
#define DEVELOP 0

#define CSP_VERSION "0.1.0"
#define CSP_AUTHOR "hxj5"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <time.h>
#include <zlib.h>
#include "thpool.h"
#include "general_util.h"
#include "cellsnp_util.h"

/* Define default values of global parameters. */
#define CSP_CHROM_ALL  {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"}
#define CSP_NCHROM     22
#define CSP_CELL_TAG   "CB"
#define CSP_UMI_TAG    "UR"
#define CSP_NTHREAD    1
#define CSP_MIN_COUNT  20
#define CSP_MIN_MAF    0.0
#define CSP_MIN_LEN    30
#define CSP_MIN_MAPQ   20
#define CSP_MAX_FLAG   255
#define CSP_OUT_VCF_CELLS   "cellSNP.cells.vcf.gz"
#define CSP_OUT_VCF_BASE    "cellSNP.base.vcf.gz"
#define CSP_OUT_SAMPLES     "cellSNP.samples.tsv"
#define CSP_OUT_MTX_AD      "cellSNP.tag.AD.mtx"
#define CSP_OUT_MTX_DP      "cellSNP.tag.DP.mtx"
#define CSP_OUT_MTX_OTH     "cellSNP.tag.OTH.mtx"

/*Structure that stores global settings/options/parameters. 
Note:
1. In current version, one and only one of out_dir and out_fn would exist and work, the other 
   would be freed. Refer to check_global_args() for details.
2. In current version, one and only one of barcode(s) and sample-ID(s) would exist and work, the other
   would be freed. Refer to check_global_args() for details.
*/
struct _gll_settings {
    char *in_fn_file;      // Name of the file containing a list of input bam/sam/cram files, one input file per line.
    int nin;               // Num of input bam/sam/cram files.
    char **in_fns;         // Pointer to the array of names of input bam/sam/cram files.
    char *out_dir;         // Pointer to the path of dir containing the output files.
    char *out_fn;          // Pointer to the path of output cells_vcf.
    char *out_vcf_cells, *out_vcf_base;
    char *out_mtx_ad, *out_mtx_dp, *out_mtx_oth;
    int is_out_zip;
    char *pos_list_file;   // Name of file containing a list of SNPs, usually a vcf file.
    csp_snplist_t pl;      // List of the input SNPs.
    char *barcode_file;    // Name of the file containing a list of barcodes, one barcode per line.
    char **barcodes;       // Pointer to the array of barcodes.
    int nbarcode;          // Num of the barcodes.
    char *sid_list_file;   // Name of the file containing a list of sample IDs, one sample-ID per line.
    char **sample_ids;     // Pointer to the array of sample IDs.
    int nsid;              // Num of sample IDs.
    char **chrom_all;      // Pointer to the array of the chromosomes to use.
    int nchrom;            // Num of chromosomes.
    char *cell_tag;        // Tag for cell barcodes, NULL means no cell tags.
    char *umi_tag;         // Tag for UMI: UR, NULL. NULL means no UMI but read counts.
    int nthread;           // Num of threads.
    threadpool tp;         // Pointer to thread pool.
    int min_count;     // Minimum aggragated count.
    double min_maf;    // Minimum minor allele frequency.
    int double_gl;     // 0 or 1. 1: keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5. 0: not keep.
    int save_hdf5;     // 0 or 1. 1: save an output file in HDF5 format, 0: not save.
    int min_len;       // Minimum mapped length for read filtering.
    int min_mapq;      // Minimum MAPQ for read filtering.
    int max_flag;      // Maximum FLAG for read filtering.
};

/*@abstract  Whether to use barcodes for sample grouping during pileup.
@param gs    Pointer of global settings structure [global_settings*].
@return      1, yes; 0, no.
*/
#define use_barcodes(gs) ((gs)->cell_tag)

/*@abstract  Whether to use sample IDs for sample grouping during pileup.
@param gs    Pointer of global settings structure [global_settings*].
@return      1, yes; 0, no.
*/
#define use_sid(gs) ((gs)->sample_ids)

/*@abstract  Whether to use UMI for reads grouping during pileup.
@param gs    Pointer of global settings structure [global_settings*].
@return      1, yes; 0, no.
*/
#define use_umi(gs) ((gs)->umi_tag)

/*@note  Do not free gs pointer itself! the system would do that! */
static void gll_setting_free(global_settings *gs) { 
    if (gs) {
        if (gs->in_fn_file) { free(gs->in_fn_file); gs->in_fn_file = NULL; }
        if (gs->in_fns) { str_arr_destroy(gs->in_fns, gs->nin); gs->in_fns = NULL; }
        if (gs->out_dir) { free(gs->out_dir); gs->out_dir = NULL; }
        if (gs->out_fn) { free(gs->out_fn); gs->out_fn = NULL; }
        if (gs->out_vcf_base) { free(gs->out_vcf_base); gs->out_vcf_base = NULL; }
        if (gs->out_vcf_cells) { free(gs->out_vcf_cells); gs->out_vcf_cells = NULL; }
        if (gs->out_mtx_ad) { free(gs->out_mtx_ad); gs->out_mtx_ad = NULL; }
        if (gs->out_mtx_dp) { free(gs->out_mtx_dp); gs->out_mtx_dp = NULL; }
        if (gs->out_mtx_oth) { free(gs->out_mtx_oth); gs->out_mtx_oth = NULL; }
        if (gs->pos_list_file) { free(gs->pos_list_file); gs->pos_list_file = NULL; }
        csp_snplist_destroy(gs->pl);
        if (gs->barcode_file) { free(gs->barcode_file); gs->barcode_file = NULL; }
        if (gs->barcodes) { str_arr_destroy(gs->barcodes, gs->nbarcode); gs->barcodes = NULL; }
        if (gs->sid_list_file) { free(gs->sid_list_file); gs->sid_list_file = NULL; }
        if (gs->sample_ids) { str_arr_destroy(gs->sample_ids, gs->nsid); gs->sample_ids = NULL; }
        if (gs->chrom_all) { str_arr_destroy(gs->chrom_all, gs->nchrom); gs->chrom_all = NULL; }
        if (gs->cell_tag) { free(gs->cell_tag); gs->cell_tag = NULL; }
        if (gs->umi_tag) { free(gs->umi_tag); gs->umi_tag = NULL; }
        if (gs->tp) { thpool_destroy(gs->tp); gs->tp = NULL; }
    }
}

/*@abstract    Set default values for global_settings structure.
@param gs      Pointer to global_settings structure returned by gll_setting_init().
@return        Void.

@note          Internal use only!
 */
static void gll_set_default(global_settings *gs) {
    if (gs) {
        gs->in_fn_file = NULL; gs->in_fns = NULL; gs->nin = 0;
        gs->out_dir = NULL; gs->out_fn = NULL;
        gs->out_vcf_base = NULL; gs->out_vcf_cells = NULL;
        gs->out_mtx_ad = NULL; gs->out_mtx_dp = NULL; gs->out_mtx_oth = NULL;
        gs->is_out_zip = 1;
        gs->pos_list_file = NULL; csp_snplist_init(gs->pl);
        gs->barcode_file = NULL; gs->nbarcode = 0; gs->barcodes = NULL; 
        gs->sid_list_file = NULL; gs->sample_ids = NULL; gs->nsid = 0;
        char *chrom_tmp[] = CSP_CHROM_ALL;
        gs->chrom_all = (char**) calloc(CSP_NCHROM, sizeof(char*));
        for (gs->nchrom = 0; gs->nchrom < CSP_NCHROM; gs->nchrom++) { gs->chrom_all[gs->nchrom] = safe_strdup(chrom_tmp[gs->nchrom]); }
        gs->cell_tag = safe_strdup(CSP_CELL_TAG); gs->umi_tag = safe_strdup(CSP_UMI_TAG);
        gs->nthread = CSP_NTHREAD; gs->tp = NULL;
        gs->min_count = CSP_MIN_COUNT; gs->min_maf = CSP_MIN_MAF; 
        gs->double_gl = 0; gs->save_hdf5 = 0; 
        gs->min_len = CSP_MIN_LEN; gs->min_mapq = CSP_MIN_MAPQ;
        gs->max_flag = CSP_MAX_FLAG;
    }
}

/* print global settings. */
static void gll_setting_print(FILE *fp, global_settings *gs, char *prefix) {
    if (gs) {
        int i;
        fprintf(fp, "%snum of input files = %d\n", prefix, gs->nin);
        fprintf(fp, "%sout_dir = %s\n", prefix, gs->out_dir);
        fprintf(fp, "%sout_fn = %s\n", prefix, gs->out_fn);
        fprintf(fp, "%sout_vcf_cells = %s\n", prefix, gs->out_vcf_base);
        fprintf(fp, "%sis_out_zip = %d\n", prefix, gs->is_out_zip);
        fprintf(fp, "%snum_of_pos = %lu\n", prefix, csp_snplist_size(gs->pl));
        fprintf(fp, "%snum_of_barcodes = %d, num_of_samples = %d\n", prefix, gs->nbarcode, gs->nsid);
        fprintf(fp, "%s%d chroms: ", prefix, gs->nchrom);
        for (i = 0; i < gs->nchrom; i++) fprintf(fp, "%s ", gs->chrom_all[i]);
        fputc('\n', fp);
        fprintf(fp, "%scell-tag = %s, umi-tag = %s\n", prefix, gs->cell_tag, gs->umi_tag);
        fprintf(fp, "%snum_of_threads = %d\n", prefix, gs->nthread);
        fprintf(fp, "%smin_count = %d, min_maf = %.2f, double_gl = %d\n", prefix, gs->min_count, gs->min_maf, gs->double_gl);
        fprintf(fp, "%ssave_hdf5 = %d, min_len = %d, min_mapq = %d\n", prefix, gs->save_hdf5, gs->min_len, gs->min_mapq);
        fprintf(fp, "%smax_flag = %d\n", prefix, gs->max_flag);
    }
}

// SZ_NUMERIC_OP_INIT(csp_s, size_t);

#define CSP_VCF_HEADER "##fileformat=VCFv4.2\n" 																							\
    "##source=cellSNP_v" CSP_VERSION "\n"																								\
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"																			\
    "##FILTER=<ID=.,Description=\"Filter info not available\">\n"																		\
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"total counts for ALT and REF\">\n"             									\
    "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"total counts for ALT\">\n"														\
    "##INFO=<ID=OTH,Number=1,Type=Integer,Description=\"total counts for other bases from REF and ALT\">\n"							\
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"																\
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"							\
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"total counts for ALT and REF\">\n"											\
    "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"total counts for ALT\">\n"														\
    "##FORMAT=<ID=OTH,Number=1,Type=Integer,Description=\"total counts for other bases from REF and ALT\">\n"						\
    "##FORMAT=<ID=ALL,Number=5,Type=Integer,Description=\"total counts for all bases in order of A,C,G,T,N\">\n"

#define CSP_VCF_CONTIG "##contig=<ID=1>\n##contig=<ID=2>\n##contig=<ID=3>\n##contig=<ID=4>\n##contig=<ID=5>\n"                             \
    "##contig=<ID=6>\n##contig=<ID=7>\n##contig=<ID=8>\n##contig=<ID=9>\n##contig=<ID=10>\n"									\
    "##contig=<ID=11>\n##contig=<ID=12>\n##contig=<ID=13>\n##contig=<ID=14>\n##contig=<ID=15>\n"									\
    "##contig=<ID=16>\n##contig=<ID=17>\n##contig=<ID=18>\n##contig=<ID=19>\n##contig=<ID=20>\n"									\
    "##contig=<ID=21>\n##contig=<ID=22>\n##contig=<ID=X>\n##contig=<ID=Y>\n"

/*@abstract  Set values for internal variables of csp_mplp_t to prepare for pileup. 
@param mplp  Pointer of csp_mplp_t structure.
@param gs    Pointer of global_settings structure.
@return      0 if success, -1 otherwise.
 */
static int csp_mplp_prepare(csp_mplp_t *mplp, global_settings *gs) {
    char **sgnames;
    int nsg;
    csp_map_sg_iter k;
    csp_plp_t *plp;
    /* init HashMap, pool of ul, pool of uu for mplp. */
    mplp->h = csp_map_sg_init();
    if (NULL == mplp->h) { fprintf(stderr, "[E::%s] could not init csp_map_sg_t structure.\n", __func__); return -1; }
    if (use_umi(gs)) {
        #if DEVELOP
            mplp->pl = csp_pool_ul_init();
            if (NULL == mplp->pl) { fprintf(stderr, "[E::%s] could not init csp_pool_ul_t structure.\n", __func__); return -1; }
            mplp->pu = csp_pool_uu_init();
            if (NULL == mplp->pu) { fprintf(stderr, "[E::%s] could not init csp_pool_uu_t structure.\n", __func__); return -1; }
        #endif
        mplp->su = csp_pool_ps_init();
        if (NULL == mplp->su) { fprintf(stderr, "[E::%s] could not init csp_pool_su_t structure.\n", __func__); return -1; }
    }
    /* set sample names for sample groups. */
    if (use_barcodes(gs)) { sgnames = gs->barcodes; nsg = gs->nbarcode; }
    else if (use_sid(gs)) { sgnames = gs->sample_ids; nsg = gs->nsid; }
    else { fprintf(stderr, "[E::%s] failed to set sample names.\n", __func__); return -1; }  // should not come here!
    if (csp_mplp_set_sg(mplp, sgnames, nsg) < 0) { fprintf(stderr, "[E::%s] failed to set sample names.\n", __func__); return -1; }
    /* init plp for each sample group in mplp->h and init HashMap plp->h for UMI grouping. */
    for (k = csp_map_sg_begin(mplp->h); k != csp_map_sg_end(mplp->h); k++) {
        if (! csp_map_sg_exist(mplp->h, k)) { continue; }
        if (NULL == (plp = csp_map_sg_val(mplp->h, k))) { 
            if (NULL == (csp_map_sg_val(mplp->h, k) = plp = csp_plp_init())) {
                fprintf(stderr, "[E::%s] failed to init csp_plp_t for sg HashMap of csp_mplp_t.\n", __func__);
                return -1;
            }
        }
        if (use_umi(gs)) {
            plp->h = csp_map_ug_init();
            if (NULL == plp->h) { fprintf(stderr, "[E::%s] could not init csp_map_ug_t structure.\n", __func__); return -1; }
        }
    }
    return 0;
}

#if DEVELOP
/*@abstract    Push content of one csp_pileup_t structure into the csp_mplp_t structure for statistics.
@param pileup  Pointer of csp_pileup_t structure to be pushed.
@param mplp    Pointer of csp_mplp_t structure pushing into.
@param sid     Index of Sample ID in the input Sample IDs.
@param gs      Pointer of global_settings structure.
@return        0 if success; Negative numbers for error:
                 -1, cell-barcode is not in input barcode-list;  
                 -2, neither barcodes or Sample IDs are used.              
                 -3, failed to calculate qual vector.
                 -5, failed to convert qual matrix to genotype.
                 -7, khash_put error.

@note   1. When all csp_pileup_t structures of one query pos have been pushed into the csp_mplp_t, call
           csp_mplp_push(NULL, mplp, sample, gs) to inform the function. And then the csp_mplp_push() 
           function will do statistics of the query pos.
        2. To speed up, the caller should guarantee that:
           a) the parameters are valid, i.e. mplp and gs must not be NULL. In fact, this function is supposed to be 
              called after csp_mplp_t is created and set names of sample-groups, so mplp, mplp->h could not be NULL.
           b) the csp_pileup_t must have passed the read filtering, refer to pileup_read_with_fetch() for details.
           c) each key (sample group name) in csp_map_sg_t already has a valid, not NULL, value (csp_plp_t*);
              This usually can be done by calling csp_mplp_prepare().
        3. This function is expected to be used by Mode1 & Mode2 & Mode3.
 */
static int csp_mplp_push(csp_pileup_t *pileup, csp_mplp_t *mplp, int sid, global_settings *gs) {
    csp_map_sg_iter k;
    csp_map_ug_iter u;
    csp_plp_t *plp = NULL;
    csp_list_uu_t *ul = NULL;
    csp_umi_unit_t *uu = NULL;
    char **s;
    int i, j, r, n, idx;
    if (NULL == pileup) {  // all csp_pileup_t have been pushed, begin to do statistics.
        for (i = 0; i < mplp->nsg; i++) {
            plp = csp_map_sg_val(mplp->h, mplp->hiter[i]);
            if (use_umi(gs)) {  // statistics based on UMI. 
                /* In current version, only the result (i.e. base and qual) of the first read of one UMI group will be used.
                   TODO: consistency correction in each UMI group. */
                for (u = csp_map_ug_begin(plp->h); u != csp_map_ug_end(plp->h); u++) {
                    if (csp_map_ug_exist(plp->h, u)) { ul = csp_map_ug_val(plp->h, u); }
                    else { continue; }
                    if (csp_list_uu_size(ul) > 0) { uu = csp_list_uu_A(ul, 0); }
                    else { continue; }
                    idx = seq_nt16_idx2int(uu->base);
                    plp->bcount[idx]++;
                    if (get_qual_vector(uu->qual, 45, 0.25, mplp->qvec) < 0) { return -3; }
                    for (j = 0; j < 4; j++) plp->qmat[idx][j] += mplp->qvec[j];				
                }
            } // else do nothing. statistics without UMI has been done before.
            for (j = 0; j < 5; j++) { 
                plp->tcount += plp->bcount[j]; 
                mplp->bc[j] += plp->bcount[j];
            }	
        }
        for (i = 0; i < 5; i++) { mplp->tc += mplp->bc[i]; }
        csp_infer_allele(mplp->bc, &mplp->inf_rid, &mplp->inf_aid);   // must be called after mplp->bc are completely calculated.
        if (mplp->ref_idx < 0 || mplp->alt_idx < 0) {  // ref or alt is not valid. Refer to csp_mplp_t.
            mplp->ref_idx = mplp->inf_rid;
            mplp->alt_idx = mplp->inf_aid;
        }
        for (i = 0; i < mplp->nsg; i++) {
            plp = csp_map_sg_val(mplp->h, mplp->hiter[i]);
            if (qual_matrix_to_geno(plp->qmat, plp->bcount, mplp->ref_idx, mplp->alt_idx, gs->double_gl, plp->gl, &plp->ngl) < 0) { return -5; }
        }
        return 0;
    }
    /* Push one csp_pileup_t into csp_mplp_t.
    *  The pileup->cb, pileup->umi could not be NULL as the pileuped read has passed filtering.
    */
    if (use_barcodes(gs)) { 
        if ((k = csp_map_sg_get(mplp->h, pileup->cb)) == csp_map_sg_end(mplp->h)) { return -1; }
        plp = csp_map_sg_val(mplp->h, k);
    } else if (use_sid(gs)) { 
        plp = csp_map_sg_val(mplp->h, mplp->hiter[sid]); 
    } else { return -2; }  // should not come here!
    if (use_umi(gs)) {
        u = csp_map_ug_get(plp->h, pileup->umi);
        if (u == csp_map_ug_end(plp->h)) {
            s = csp_pool_ps_get(mplp->su);
            *s = strdup(pileup->umi);
            u = csp_map_ug_put(plp->h, *s, &r);
            if (r < 0) { return -7; }
            csp_map_ug_val(plp->h, u) = ul = csp_pool_ul_get(mplp->pl);
        } else { ul = csp_map_ug_val(plp->h, u); }
        uu = csp_pool_uu_get(mplp->pu);
        uu->base = pileup->base; uu->qual = pileup->qual;
        csp_list_uu_push(ul, uu);
    } else {
        idx = seq_nt16_idx2int(pileup->base);
        plp->bcount[idx]++;
        if (get_qual_vector(pileup->qual, 45, 0.25, mplp->qvec) < 0) { return -3; }
        for (i = 0; i < 4; i++) plp->qmat[idx][i] += mplp->qvec[i];
    }
    return 0;
}
#else
/*@abstract    Push content of one csp_pileup_t structure into the csp_mplp_t structure for statistics.
@param pileup  Pointer of csp_pileup_t structure to be pushed.
@param mplp    Pointer of csp_mplp_t structure pushing into.
@param sid     Index of Sample ID in the input Sample IDs.
@param gs      Pointer of global_settings structure.
@return        0 if success; Negative numbers for error:
                 -1, cell-barcode is not in input barcode-list;  
                 -2, neither barcodes or Sample IDs are used.              
                 -3, failed to calculate qual vector.
                 -5, failed to convert qual matrix to genotype.
                 -7, khash_put error.

@note   1. When all csp_pileup_t structures of one query pos have been pushed into the csp_mplp_t, call
           csp_mplp_push(NULL, mplp, sample, gs) to inform the function. And then the csp_mplp_push() 
           function will do statistics of the query pos.
        2. To speed up, the caller should guarantee that:
           a) the parameters are valid, i.e. mplp and gs must not be NULL. In fact, this function is supposed to be 
              called after csp_mplp_t is created and set names of sample-groups, so mplp, mplp->h could not be NULL.
           b) the csp_pileup_t must have passed the read filtering, refer to pileup_read_with_fetch() for details.
           c) each key (sample group name) in csp_map_sg_t already has a valid, not NULL, value (csp_plp_t*);
              This usually can be done by calling csp_mplp_prepare().
        3. This function is expected to be used by Mode1 & Mode2 & Mode3.

@discuss  In current version, only the result (base and qual) of the first read in one UMI group will be used for mplp statistics.
          TODO: store results of all reads in one UMI group (maybe could do consistency correction in each UMI group) and then 
          do mplp statistics.
 */
static int csp_mplp_push(csp_pileup_t *pileup, csp_mplp_t *mplp, int sid, global_settings *gs) {
    csp_map_sg_iter k;
    csp_map_ug_iter u;
    csp_plp_t *plp = NULL;
    char **s;
    int i, j, r, idx;
    if (NULL == pileup) {  // all csp_pileup_t have been pushed, begin to do statistics.
        for (i = 0; i < mplp->nsg; i++) {
            plp = csp_map_sg_val(mplp->h, mplp->hiter[i]);
            for (j = 0; j < 5; j++) { 
                plp->tcount += plp->bcount[j]; 
                mplp->bc[j] += plp->bcount[j];
            }	
        }
        for (i = 0; i < 5; i++) { mplp->tc += mplp->bc[i]; }
        csp_infer_allele(mplp->bc, &mplp->inf_rid, &mplp->inf_aid);   // must be called after mplp->bc are completely calculated.
        if (mplp->ref_idx < 0 || mplp->alt_idx < 0) {  // ref or alt is not valid. Refer to csp_mplp_t.
            mplp->ref_idx = mplp->inf_rid;
            mplp->alt_idx = mplp->inf_aid;
        }
        for (i = 0; i < mplp->nsg; i++) {
            plp = csp_map_sg_val(mplp->h, mplp->hiter[i]);
            if (qual_matrix_to_geno(plp->qmat, plp->bcount, mplp->ref_idx, mplp->alt_idx, gs->double_gl, plp->gl, &plp->ngl) < 0) { return -5; }
        }
        return 0;
    }
    /* Push one csp_pileup_t into csp_mplp_t.
    *  The pileup->cb, pileup->umi could not be NULL as the pileuped read has passed filtering.
    */
    if (use_barcodes(gs)) { 
        if ((k = csp_map_sg_get(mplp->h, pileup->cb)) == csp_map_sg_end(mplp->h)) { return -1; }
        plp = csp_map_sg_val(mplp->h, k);
    } else if (use_sid(gs)) { 
        plp = csp_map_sg_val(mplp->h, mplp->hiter[sid]);
    } else { return -2; }  // should not come here!
    if (use_umi(gs)) {
        u = csp_map_ug_get(plp->h, pileup->umi);
        if (u == csp_map_ug_end(plp->h)) {
            s = csp_pool_ps_get(mplp->su);
            *s = strdup(pileup->umi);
            u = csp_map_ug_put(plp->h, *s, &r);
            if (r < 0) { return -7; }
            idx = seq_nt16_idx2int(pileup->base);
            plp->bcount[idx]++;
            if (get_qual_vector(pileup->qual, 45, 0.25, mplp->qvec) < 0) { return -3; }
            for (j = 0; j < 4; j++) plp->qmat[idx][j] += mplp->qvec[j];				
        } // else: do nothing.
    } else {
        idx = seq_nt16_idx2int(pileup->base);
        plp->bcount[idx]++;
        if (get_qual_vector(pileup->qual, 45, 0.25, mplp->qvec) < 0) { return -3; }
        for (i = 0; i < 4; i++) plp->qmat[idx][i] += mplp->qvec[i];
    }
    return 0;
}
#endif

/*@abstract  Pileup one read obtained by sam_itr_next().
@param pos   Pos of the reference sequence. 0-based.
@param p     Pointer of csp_pileup_t structure coming from csp_pileup_init() or csp_pileup_reset().
@param gs    Pointer of global settings.
@return      0 if success, negative numbers otherwise:
               -1, invalid umi.
               -3, invalid cell-barcode.
               -5, low mapq.
               -7, high sam mapping flag.
               -9, is in del.
               -11, is in refskip.
               -13, not enough length of bases within alignment.

@note        1. This function is modified from cigar_resolve2() function in sam.c of htslib.
             2. Reads filtering is also applied inside this function, including:
                   UMI and cell tags, read mapping quality, mapping flag and length of bases within alignment.
             3. To speed up, parameters will not be checked, so the caller should guarantee the parameters are valid, i.e.
                && p != NULL && gs != NULL.

@TODO        Filter unmapped reads (the read itself unmapped or the mate read unmapped) ?
 */
static int pileup_read_with_fetch(hts_pos_t pos, csp_pileup_t *p, global_settings *gs) {
    /* Filter reads in order. For example, filtering according to umi tag and cell tag would speed up in the case
       that do not use UMI or Cell-barcode at all. */
    if (use_umi(gs) && NULL == (p->umi = get_bam_aux_str(p->b, gs->umi_tag))) { return -1; }
    if (use_barcodes(gs) && NULL == (p->cb = get_bam_aux_str(p->b, gs->cell_tag))) { return -3; }
    bam1_core_t *c = &(p->b->core);
    if (c->qual < gs->min_mapq) { return -5; }
    if (c->flag > gs->max_flag) { return -7; }
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
    if (p->is_del) { return -9; }
    if (p->is_refskip) { return -11; }
    /* continue processing cigar string. */
    for (k++; k < c->n_cigar; k++) {
        op = get_cigar_op(cigar[k]);
        l = get_cigar_len(cigar[k]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) { laln += l; }
    }
    if (laln < gs->min_len) { return -13; }
    else { p->laln = laln; }
    return 0;
}

/*@abstract    Pileup one SNP with method fetch.
@param snp     Pointer of csp_snp_t structure.
@param fs      Pointer of array of pointers to the csp_bam_fs structures.
@param nfs     Size of @p fs.
@param pileup  Pointer of csp_pileup_t structure.
@param mplp    Pointer of csp_mplp_t structure.
@param gs      Pointer of global_settings structure.
@param s       Pointer of kstring_t struture.
@return        0 if success, negative number for error:
                 -1, unknown error.
                 -2, cannot translate chr name to tid.
                 -3, sam_itr_queryi() failed.
                 -5, barcodes and sid are not provided.
                 -7, failed to push pileup into mplp.
                 -9, error reading input bam/sam/cram.
                 -11, no pileup pushed.
                 -13, failed to do statistics for mplp.
                 -15, not passing min_count.
                 -17, not passing min_maf.

@note          1. This function is mainly called by pileup_positions_with_fetch(). Refer to pileup_positions_with_fetch() for notes.
               2. The statistics result of all pileuped reads for one SNP is stored in the csp_mplp_t after calling this function.
*/
static int pileup_snp_with_fetch(csp_snp_t *snp, csp_bam_fs **fs, int nfs, csp_pileup_t *pileup, csp_mplp_t *mplp, global_settings *gs, kstring_t *s) 
{
    csp_bam_fs *bs = NULL;
    hts_itr_t *iter = NULL;
    int i, tid, r, ret, state = -1;
    size_t npushed = 0;
    #if DEBUG
        size_t npileup = 0;
    #endif
    mplp->ref_idx = snp->ref ? seq_nt16_char2int(snp->ref) : -1;
    mplp->alt_idx = snp->alt ? seq_nt16_char2int(snp->alt) : -1;
    for (i = 0; i < nfs; i++) {
        bs = fs[i];
        tid = csp_sam_hdr_name2id(bs->hdr, snp->chr, s);   // the kstring_t will be cleared inside the csp_sam_hdr_name2id().
        if (tid < 0) { state = -2; goto fail; }
        if (NULL == (iter = sam_itr_queryi(bs->idx, tid, snp->pos, snp->pos + 1))) { state = -3; goto fail; }
        while ((ret = sam_itr_next(bs->fp, iter, pileup->b)) >= 0) {   // TODO: check if need to be reset in_fp?
            #if DEBUG
                npileup++;
            #endif
            if (pileup_read_with_fetch(snp->pos, pileup, gs) >= 0) {
                if (use_barcodes(gs)) { r = csp_mplp_push(pileup, mplp, -1, gs); }
                else if (use_sid(gs)) { r = csp_mplp_push(pileup, mplp, i, gs); }
                else { state = -5; goto fail; }
                if (r < -1) { state = -7; goto fail; }  // else if r == -1: pileuped barcode is not in the input barcode list.
                else if (r >= 0) { npushed++; }
            } // no need to reset pileup as the values in it will be immediately overwritten.
        }
        if (ret < -1) { state = -9; goto fail; } 
        else { hts_itr_destroy(iter); iter = NULL; }  // TODO: check if could reset iter?
    }
    #if DEBUG
        fprintf(stderr, "[D::%s] before mplp statistics: npileup = %ld; npushed = %ld; the mplp is:\n", __func__, npileup, npushed);
        csp_mplp_print_(stderr, mplp, "\t");
    #endif
    if (npushed <= 0) { state = -11; goto fail; }
    if (csp_mplp_push(NULL, mplp, -1, gs) < 0) { state = -13; goto fail; }
    #if DEBUG
        fprintf(stderr, "[D::%s] after mplp statistics: the mplp is:\n", __func__);
        csp_mplp_print_(stderr, mplp, "\t");
    #endif
    if (mplp->tc < gs->min_count) { state = -15; goto fail; }
    if (mplp->bc[mplp->inf_aid] < mplp->tc * gs->min_maf) { state = -17; goto fail; }
    return 0;
  fail:
    if (iter) { hts_itr_destroy(iter); }
    return state;
}

/*@abstract  Pileup a region (a list of SNPs) with method of fetching.
@param args  Pointer to thread_data structure.
@return      Num of SNPs, including those filtered, that are processed.

@note        1. Differ from pileup method in samtools, this function fetches reads covering the SNPs and 
                pileups the reads by processing CIGAR strings with a self-defined resolver function.
             2. The internal variable "ret" in thread_data structure saves the running state of the function:
                  0 if success, 
                  -1 otherwise.
             3. This function could be used by Mode1 and Mode3.		 
 */
static size_t pileup_positions_with_fetch(void *args) {
#define TMP_BUFSIZE 8388608
    thread_data *d = (thread_data*) args;
    global_settings *gs = d->gs;
    csp_snp_t **a = gs->pl.a + d->n;  /* here we use directly the internal variables in csp_snplist_t structure to speed up. */
    size_t n = 0;                     /* The num of SNPs that are successfully processed. */
    csp_bam_fs **bam_fs = NULL;       /* use array instead of single element to compatible with multi-input-files. */
    int nfs = 0;
    csp_bam_fs *bs = NULL;
    gzFile out_zfp = NULL;
    FILE *out_fp = NULL;
    csp_pileup_t *pileup = NULL;
    csp_mplp_t *mplp = NULL;
    int i, ret;
    kstring_t ks = KS_INITIALIZE;
    kstring_t *s = &ks;
    kstring_t kbuf = KS_INITIALIZE;
    kstring_t *buf = &kbuf;
#if DEBUG
    fprintf(stderr, "[D::%s][Thread-%d] thread options:\n", __func__, d->i);
    thdata_print(stderr, d);
#endif
    d->ret = -1;
    /* prepare data and structures. 
    */
    if (d->is_out_zip) {
        out_zfp = gzopen(d->out_fn, d->out_fm);
        if (NULL == out_zfp) { fprintf(stderr, "[E::%s] could not open output tmp file '%s'\n", __func__, d->out_fn); return 0; }
    } else {
        out_fp = fopen(d->out_fn, d->out_fm);
        if (NULL == out_fp) { fprintf(stderr, "[E::%s] could not open output tmp file '%s'\n", __func__, d->out_fn); return 0; }
    }
    /* prepare mplp for pileup. */
    if (NULL == (mplp = csp_mplp_init())) { fprintf(stderr, "[E::%s] could not init csp_mplp_t structure.\n", __func__); goto fail; }
    if (csp_mplp_prepare(mplp, gs) < 0) { fprintf(stderr, "[E::%s] could not prepare csp_mplp_t structure.\n", __func__); goto fail; }
    /* create file structures for input bam/sam/cram. */
    bam_fs = (csp_bam_fs**) calloc(gs->nin, sizeof(csp_bam_fs*));  	
    if (NULL == bam_fs) { fprintf(stderr, "[E::%s] could not initialize csp_bam_fs array.\n", __func__); goto fail; }
    for (nfs = 0; nfs < gs->nin; nfs++) {
        if (NULL == (bs = csp_bam_fs_build(gs->in_fns[nfs], &ret))) {
            fprintf(stderr, "[E::%s] could not build csp_bam_fs structure.\n", __func__);
            goto fail;
        } else { bam_fs[nfs] = bs; }
    }
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
    for (n = 0; n < d->m; n++) {		
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
        if ((ret = pileup_snp_with_fetch(a[n], bam_fs, nfs, pileup, mplp, gs, s)) < 0) {
            #if DEBUG
                fprintf(stderr, "[W::%s] snp (%s:%ld) filtered, error code = %d\n", __func__, a[n]->chr, a[n]->pos + 1, ret);
            #endif
            csp_mplp_reset(mplp);
            continue;
        }
        size_t dp = mplp->bc[mplp->ref_idx] + mplp->bc[mplp->alt_idx];
        if (csp_mplp_to_vcf(mplp) < 0) { 
            fprintf(stderr, "[E::%s] failed to convert mplp to vcf format for snp (%s:%ld).\n", __func__, a[n]->chr, a[n]->pos + 1);
            csp_mplp_reset(mplp);
            continue;
        }
        ksprintf(buf, "%s\t%ld\t.\t%c\t%c\t.\tPASS\tAD=%ld;DP=%ld;OTH=%ld\tGT:AD:DP:OTH:PL:ALL%s\n", \
                    a[n]->chr, a[n]->pos + 1, seq_nt16_int2char(mplp->ref_idx), seq_nt16_int2char(mplp->alt_idx), \
                    mplp->bc[mplp->alt_idx], dp, mplp->tc - dp, mplp->s.s);
        if (ks_len(buf) >= TMP_BUFSIZE) { 
            if (d->is_out_zip) { gzwrite(out_zfp, ks_str(buf), ks_len(buf)); ks_clear(buf); }
            else { fwrite(ks_str(buf), 1, ks_len(buf), out_fp); ks_clear(buf); }
        }
        csp_mplp_reset(mplp);
    }
    if (ks_len(buf)) {
        if (d->is_out_zip) { gzwrite(out_zfp, ks_str(buf), ks_len(buf)); ks_clear(buf); }
        else { fwrite(ks_str(buf), 1, ks_len(buf), out_fp); }
    }
    ks_free(s);
    ks_free(buf);
    csp_pileup_destroy(pileup);
    for (i = 0; i < nfs; i++) csp_bam_fs_destroy(bam_fs[i]);
    free(bam_fs);
    csp_mplp_destroy(mplp);
    if (out_fp) fclose(out_fp);
    if (out_zfp) gzclose(out_zfp);
    d->ret = 0;
    return n;
  fail:
    ks_free(s);
    ks_free(buf);
    if (pileup) csp_pileup_destroy(pileup);
    if (bam_fs) {
        for (i = 0; i < nfs; i++) csp_bam_fs_destroy(bam_fs[i]);
        free(bam_fs);		
    }
    if (mplp) { csp_mplp_destroy(mplp); }
    if (out_fp) fclose(out_fp);
    if (out_zfp) gzclose(out_zfp);
    return n;
#undef TMP_BUFSIZE
}

/*abstract  Run cellSNP Mode with method of fetching.
@param gs   Pointer to the global_settings structure.
@return     0 if success, -1 otherwise.
 */
static int run_mode_with_fetch(global_settings *gs) {
    /* check options (input) */
    if (NULL == gs || gs->nin <= 0 || gs->nbarcode <= 0 || csp_snplist_size(gs->pl) <= 0 || \
        ! ((NULL == gs->out_dir) ^ (NULL == gs->out_fn)) || NULL == gs->out_vcf_cells) {
        fprintf(stderr, "[E::%s] error options for mode1.\n", __func__);
        return -1;
    }
    /* output vcf headers. */
    kstring_t out_buf = KS_INITIALIZE;
    kstring_t *buf = &out_buf;
    int k;
    char *out_fn = gs->out_vcf_cells;
    gzFile out_zfp = NULL;
    FILE *out_fp = NULL;
    if (gs->is_out_zip) { 
        out_zfp = gzopen(out_fn, "wb"); 
        if (NULL == out_zfp) { fprintf(stderr, "[E::%s] cannot open '%s'\n", __func__, out_fn); goto fail_main; }
    } else { 
        out_fp = fopen(out_fn, "w");
        if (NULL == out_fp) { fprintf(stderr, "[E::%s] cannot open '%s'\n", __func__, out_fn); goto fail_main; }
    }
    kputs(CSP_VCF_HEADER, buf);
    kputs(CSP_VCF_CONTIG, buf);
    kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", buf);
    if (use_barcodes(gs) && gs->barcodes) {
        for (k = 0; k < gs->nbarcode; k++) { kputc_('\t', buf); kputs(gs->barcodes[k], buf); }
    } else if (use_sid(gs) && gs->sample_ids) {
        for (k = 0; k < gs->nsid; k++) { kputc_('\t', buf); kputs(gs->sample_ids[k], buf); }
    } else { fprintf(stderr, "[E::%s] neither barcodes or sample IDs exist.\n", __func__); goto fail_main; }
    kputc('\n', buf);
    if (gs->is_out_zip) {
        if (0 == gzwrite(out_zfp, ks_str(buf), ks_len(buf))) { fprintf(stderr, "[E::%s] gzwrite error.\n", __func__); goto fail_main; }
    } else if (fwrite(ks_str(buf), 1, ks_len(buf), out_fp) != ks_len(buf)) { 
        fprintf(stderr, "[E::%s] fwrite error.\n", __func__); 
        goto fail_main; 
    }
    ks_free(buf); buf = NULL;
    /* core part of Mode 1. */
    if (gs->tp && gs->nthread > 1) {
        char **out_tmp_fn = NULL;
        int nout_tmp = 0, mout_tmp = 0;
        thread_data **td = NULL, *d = NULL;
        int ntd = 0, mtd = gs->nthread; // ntd: num of thread-data structures that have been created. mtd: size of td array.
        int i, n;
        size_t npos, mpos;
        kstring_t ks = KS_INITIALIZE;
        kstring_t *s = &ks;
        /* create output tmp filenames. */
        if (NULL == (out_tmp_fn = (char**) calloc(gs->nthread, sizeof(char*)))) {
            fprintf(stderr, "[E::%s] could not allocate space for out_tmp_fn.\n", __func__);
            ks_free(s);
            goto fail;
        } else { mout_tmp = gs->nthread; }
        for (nout_tmp = 0; nout_tmp < mout_tmp; nout_tmp++) {  /* file path delim '/' only works for Unix? */
            ksprintf(s, "%s.%d", out_fn, nout_tmp); out_tmp_fn[nout_tmp] = strdup(ks_str(s)); ks_clear(s);
        } ks_free(s);
        if (nout_tmp < mout_tmp) { 
            fprintf(stderr, "[E::%s] could not create all out_tmp files, NO. %d failed.\n", __func__, nout_tmp); 
            goto fail; 
        }	
        /* prepare data for thread pool and run. */
        td = (thread_data**) calloc(mtd, sizeof(thread_data*));
        if (NULL == td) { fprintf(stderr, "[E::%s] could not initialize the array of thread_data structure.\n", __func__); goto fail; }
        for (npos = 0, mpos = csp_snplist_size(gs->pl) / mtd; ntd < mtd - 1; ntd++, npos += mpos) { /* mtd is equal to mout_tmp */
            if (NULL == (d = thdata_init())) { 
                fprintf(stderr, "[E::%s] could not initialize the thread_data structure.\n", __func__); 
                goto fail; 
            }
            d->gs = gs; d->n = npos; d->m = mpos; d->i = ntd; 
            d->out_fn = out_tmp_fn[ntd]; d->out_fm = "w"; d->is_out_zip = 0;
            if (thpool_add_work(gs->tp, (void*) pileup_positions_with_fetch, d) < 0) {
                fprintf(stderr, "[E::%s] could not add thread work (No. %d)\n", __func__, ntd++);
                goto fail;
            } else { td[ntd] = d; }
        }
        if (csp_snplist_size(gs->pl) - npos > 0) { // still have some SNPs to be processed.
            if (NULL == (d = thdata_init())) { 
                fprintf(stderr, "[E::%s] could not initialize the thread_data structure.\n", __func__); 
                goto fail; 
            }
            d->gs = gs; d->n = npos; d->m = csp_snplist_size(gs->pl) - npos; d->i = ntd; 
            d->out_fn = out_tmp_fn[ntd]; d->out_fm = "w"; d->is_out_zip = 0;
            if (thpool_add_work(gs->tp, (void*) pileup_positions_with_fetch, d) < 0) {
                fprintf(stderr, "[E::%s] could not add thread work (No. %d)\n", __func__, ntd++);
                goto fail;
            } else { td[ntd++] = d; }
        }
        thpool_wait(gs->tp);
        /* check running status of threads. */
        #if DEBUG
            for (i = 0; i < ntd; i++) { fprintf(stderr, "[D::%s] ret of thread-%d is %d\n", __func__, i, td[i]->ret); }
        #endif
        for (i = 0; i < ntd; i++) { if (td[i]->ret < 0) goto fail; } 
        /* clean */ 
        for (i = 0; i < ntd; i++) { thdata_destroy(td[i]); }
        free(td); td = NULL;
        if (gs->is_out_zip) { n = merge_files_to_zip_fp(out_tmp_fn, nout_tmp, out_zfp); }   /* merge tmp files */
        else { n = merge_files_to_fp(out_tmp_fn, nout_tmp, out_fp); }
        if (n < nout_tmp) { fprintf(stderr, "[E::%s] could not merge files.\n", __func__); goto fail; }
        if (remove_files(out_tmp_fn, nout_tmp) < nout_tmp) { goto fail; }  /* clean tmp files. */
        str_arr_destroy(out_tmp_fn, nout_tmp);
        if (out_fp) fclose(out_fp);
        if (out_zfp) gzclose(out_zfp);
        return 0;
      fail:
        if (td) {
            for (i = 0; i < ntd; i++) { thdata_destroy(td[i]); }
            free(td);
        }
        if (out_tmp_fn) { remove_files(out_tmp_fn, nout_tmp); str_arr_destroy(out_tmp_fn, nout_tmp); } /* clean tmp files. */
        if (out_fp) fclose(out_fp);
        if (out_zfp) gzclose(out_zfp);
        return -1;
    } else if (1 == gs->nthread) {  // only one thread.
        thread_data *d = NULL;
        if (NULL == (d = thdata_init())) { 
            fprintf(stderr, "[E::%s] could not initialize the thread_data structure.\n", __func__); 
            goto fail_single; 
        }
        d->gs = gs; d->n = 0; d->m = csp_snplist_size(gs->pl); d->i = 0; 
        d->out_fn = out_fn; d->out_fm = "a"; d->is_out_zip = gs->is_out_zip;
        if (gs->is_out_zip) { gzclose(out_zfp); out_zfp = NULL; }
        else { fclose(out_fp); out_fp = NULL; }
        pileup_positions_with_fetch(d);
        if (d->ret < 0) goto fail_single;
        thdata_destroy(d); d = NULL;
        if (out_fp) fclose(out_fp);
        if (out_zfp) gzclose(out_zfp);
        return 0;
        fail_single:
        if (d) { thdata_destroy(d); }
        if (out_fp) fclose(out_fp);
        if (out_zfp) gzclose(out_zfp);
        return -1;
    } /* else: do nothing. should not come here! */
  fail_main:
    ks_free(buf);
    if (out_fp) fclose(out_fp);
    if (out_zfp) gzclose(out_zfp);
    return -1;
}

static inline int run_mode1(global_settings *gs) { return run_mode_with_fetch(gs); }

static int run_mode2(global_settings *gs) {
    return 0;
}

static inline int run_mode3(global_settings *gs) { return run_mode_with_fetch(gs); }

static void print_usage(FILE *fp) {
    fprintf(fp, 
"\n"
"Usage: cellsnp [options]\n"
"\n"
"Options:\n"
"  -h, --help           Show this help message and exit.\n"
"  -s, --samFile STR    Indexed sam/bam file(s), comma separated multiple samples.\n"
"                       Mode 1&2: one sam/bam file with single cell.\n"
"                       Mode 3: one or multiple bulk sam/bam files,\n"
"                       no barcodes needed, but sample ids and regionsVCF.\n"
"  -S, --samFileList FILE   A list file containing bam files, each per line, for Mode 3.\n"
"  -O, --outDir DIR         Output directory for VCF and sparse matrices: AD, DP, OTH.\n"
"  -o, --outVCF FILE        Output full path with file name for VCF file. Only use if not given outDir.\n"
"  -R, --regionsVCF FILE    A vcf file listing all candidate SNPs, for fetch each variants.\n" 
"                           If None, pileup the genome. Needed for bulk samples.\n"
"  -b, --barcodeFile FILE   A plain file listing all effective cell barcode.\n"
"  -i, --sampleList FILE    A list file containing sample IDs, each per line.\n"
"  -I, --sampleIDs STR      Comma separated sample ids.\n");
    fprintf(fp,
"\n"
"Optional arguments:\n"
"  -p, --nproc INT      Number of subprocesses [%d]\n", CSP_NTHREAD);
    fprintf(fp,
"  --chrom STR          The chromosomes to use, comma separated [1 to %d]\n", CSP_NCHROM);
    fprintf(fp,
"  --cellTAG STR        Tag for cell barcodes, turn off with None [%s]\n", CSP_CELL_TAG);
    fprintf(fp,
"  --UMItag STR         Tag for UMI: UR, Auto, None. For Auto mode, use UR if barcodes is inputted,\n"
"                       otherwise use None. None mode means no UMI but read counts [%s]\n", CSP_UMI_TAG);
    fprintf(fp,
"  --minCOUNT INT       Minimum aggragated count [%d]\n", CSP_MIN_COUNT);
    fprintf(fp,
"  --minMAF FLOAT       Minimum minor allele frequency [%.2f]\n", CSP_MIN_MAF);
    fprintf(fp,
"  --doubletGL          If use, keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5.\n"
"  --saveHDF5           If use, save an output file in HDF5 format.\n"
"\n"
"Read filtering:\n"
"  --minLEN INT         Minimum mapped length for read filtering [%d]\n", CSP_MIN_LEN);
    fprintf(fp,
"  --minMAPQ INT        Minimum MAPQ for read filtering [%d]\n", CSP_MIN_MAPQ);
    fprintf(fp,
"  --maxFLAG INT        Maximum FLAG for read filtering [%d]\n", CSP_MAX_FLAG);
    fputc('\n', fp);
}

static inline int csp_cmp_barcodes(const void *x, const void *y) {
    return strcmp(*((char**) x), *((char**) y));
}

/*@abstract    Perform basic check for global settings right after running getopt()/getopt_long() function.
@param gs      Pointer to the global settings.
@return        0 if no error, negative numbers otherwise:
                 -1, should print_usage after return.
                 -2, no action.

@note          This is just basic check for the shared parameters of different running modes.
               More careful and personalized check would be performed by each running mode.
 */
static int check_global_args(global_settings *gs) {
    int i;
    if (gs->in_fn_file) {
        if (gs->in_fns) { 
            fprintf(stderr, "[E::%s] should not specify -s/--samFile and -S/--samFileList options at the same time.\n", __func__); 
            return -1; 
        } else if (NULL == (gs->in_fns = hts_readlines(gs->in_fn_file, &gs->nin)) || gs->nin <= 0) {
            fprintf(stderr, "[E::%s] could not read '%s'\n", __func__, gs->in_fn_file); 
            return -2;
        }
    } else if (NULL == gs->in_fns) { 
        fprintf(stderr, "[E::%s] should specify -s/--samFile or -S/--samFileList option.\n", __func__); 
        return -1; 
    }
    for (i = 0; i < gs->nin; i++) {
        if (0 != access(gs->in_fns[i], F_OK)) { fprintf(stderr, "[E::%s] '%s' does not exist.\n", __func__, gs->in_fns[i]); return -2; }
    }
    /* In current version, one and only one of out_dir and out_fn would exist and work. */
    if (gs->out_dir) {
        if (gs->out_fn) { 
            fprintf(stderr, "[E::%s] should not specify -o/--outVCF and -O/--outDir options at the same time.\n", __func__); 
            return -1; 
        } else if (0 != access(gs->out_dir, F_OK) && 0 != mkdir(gs->out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) { 
            fprintf(stderr, "[E::%s] '%s' does not exist.\n", __func__, gs->out_dir); 
            return -2; 
        } else { 
            gs->out_vcf_base = join_path(gs->out_dir, CSP_OUT_VCF_BASE);  
            gs->out_vcf_cells = join_path(gs->out_dir, CSP_OUT_VCF_CELLS);
            gs->out_mtx_ad = join_path(gs->out_dir, CSP_OUT_MTX_AD); 
            gs->out_mtx_dp = join_path(gs->out_dir, CSP_OUT_MTX_DP);
            gs->out_mtx_oth = join_path(gs->out_dir, CSP_OUT_MTX_OTH);
        }
    } else if (gs->out_fn) {  // the pointer returned by dirname() should not be freed?
        gs->out_vcf_cells = strdup(gs->out_fn);
    } else { fprintf(stderr, "[E::%s] should specify -o/--outVCF or -O/--outDir option.\n", __func__); return -1; }
     /* 1. In current version, one and only one of barcodes and sample-ids would exist and work. Prefer barcodes. 
        2. For barcodes, the barcode file would not be read unless cell-tag is set, i.e. the barcodes and cell-tag are
           effective only when both of them are valid. */
    if (gs->cell_tag && (0 == strcmp(gs->cell_tag, "None") || 0 == strcmp(gs->cell_tag, "none"))) { 
        free(gs->cell_tag);  gs->cell_tag = NULL; 
    }
    if (gs->cell_tag && gs->barcode_file) {
        if (gs->sample_ids || gs->sid_list_file) { 
            fprintf(stderr, "[E::%s] should not specify barcodes and sample IDs at the same time.\n", __func__); 
            return -1; 
        } else if (NULL == (gs->barcodes = hts_readlines(gs->barcode_file, &gs->nbarcode))) {
            fprintf(stderr, "[E::%s] could not read barcode file '%s'\n", __func__, gs->barcode_file); 
            return -2;
        } else { qsort(gs->barcodes, gs->nbarcode, sizeof(char*), csp_cmp_barcodes); }
    } else if ((NULL == gs->cell_tag) ^ (NULL == gs->barcode_file)) {
        fprintf(stderr, "[E::%s] should not specify barcodes or cell-tag alone.\n", __func__); 
        return -1;
    } else {
        if (NULL == gs->sample_ids) {
            if (NULL == gs->sid_list_file) { 
                kstring_t ks = KS_INITIALIZE;
                kstring_t *s = &ks;
                for (i = 0; i < gs->nin; i++) { ksprintf(s, "Sample_%d", i); gs->sample_ids[i] = strdup(ks_str(s)); ks_clear(s); }
                ks_free(s);
            } else if (NULL == (gs->sample_ids = hts_readlines(gs->sid_list_file, &gs->nsid))) {
                fprintf(stderr, "[E::%s] could not read '%s'\n", __func__, gs->sid_list_file);  
                return -2;
            } // else: sort sample ids and corresponded input-bam-files?
        } else if (gs->sid_list_file) { 
            fprintf(stderr, "[E::%s] should not specify -i/--samileList and -I/--sampleIDs options at the same time.\n", __func__);
            return -1; 
        } // else do nothing.
        if (gs->nin != gs->nsid) {
            fprintf(stderr, "[E::%s] num of sample IDs (%d) is not equal with num of input bam/sam/cram files (%d).\n", __func__, gs->nsid, gs->nin);
            return -2;
        }
    }
    /* 1. In current version, one and only one of pos_list and chrom(s) would exist and work. Prefer pos_list. 
       2. Sometimes, pos_list_file and chrom_all are both not NULL as the chrom_all has been set to default value when
          global_settings structure was just created. In this case, free chrom_all and save pos_list_file. */
    if (NULL == gs->pos_list_file || 0 == strcmp(gs->pos_list_file, "None") || 0 == strcmp(gs->pos_list_file, "none")) { 
        if (NULL == gs->chrom_all) { fprintf(stderr, "[E::%s] should specify -R/--regionsVCF or --chrom option.\n", __func__); return -1; }
        if (gs->pos_list_file) { free(gs->pos_list_file); gs->pos_list_file = NULL; }
    } else if (gs->chrom_all) { str_arr_destroy(gs->chrom_all, gs->nchrom); gs->chrom_all = NULL; gs->nchrom = 0; }
    if (gs->umi_tag) {
        if (0 == strcmp(gs->umi_tag, "Auto")) {
            if (gs->barcodes) { free(gs->umi_tag); gs->umi_tag = strdup("UR"); }
            else { free(gs->umi_tag); gs->umi_tag = NULL; }
        } else if (0 == strcmp(gs->umi_tag, "None") || 0 == strcmp(gs->umi_tag, "none")) { free(gs->umi_tag); gs->umi_tag = NULL; }
    }
    return 0;
}

int main(int argc, char **argv) {
    /* timing */
    time_t start_time, end_time;
    struct tm *time_info;
    char time_str[30];
    time(&start_time);
    time_info = localtime(&start_time);
    strftime(time_str, 30, "%Y-%m-%d %H:%M:%S", time_info);
    /* Formal part */
    global_settings gs;
    gll_set_default(&gs);
    int c, ret, print_time = 1;
    struct option lopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"samFile", required_argument, NULL, 's'},
        {"samfile", required_argument, NULL, 's'},
        {"samFileList", required_argument, NULL, 'S'},			
        {"outDir", required_argument, NULL, 'O'},
        {"outdir", required_argument, NULL, 'O'},
        {"outVCF", required_argument, NULL, 'o'},
        {"outvcf", required_argument, NULL, 'o'},
        {"regionsVCF", required_argument, NULL, 'R'},
        {"regionsvcf", required_argument, NULL, 'R'},
        {"barcodeFile", required_argument, NULL, 'b'},
        {"barcodefile", required_argument, NULL, 'b'},
        {"sampleList", required_argument, NULL, 'i'},
        {"sampleIDs", required_argument, NULL, 'I'},
        {"sampleids", required_argument, NULL, 'I'},
        {"nproc", required_argument, NULL, 'p'},
        {"chrom", required_argument, NULL, 1},
        {"cellTAG", required_argument, NULL, 2},
        {"celltag", required_argument, NULL, 2},
        {"UMItag", required_argument, NULL, 3},
        {"umitag", required_argument, NULL, 3},
        {"minCOUNT", required_argument, NULL, 4},
        {"minCount", required_argument, NULL, 4},
        {"mincount", required_argument, NULL, 4},
        {"minMAF", required_argument, NULL, 5},
        {"doubleGL", no_argument, NULL, 6},
        {"saveHDF5", no_argument, NULL, 7},
        {"minLEN", required_argument, NULL, 8},
        {"minLen", required_argument, NULL, 8},
        {"minlen", required_argument, NULL, 8},
        {"minMAPQ", required_argument, NULL, 9},
        {"maxFLAG", required_argument, NULL, 10},
        {"maxFlag", required_argument, NULL, 10},
        {"maxflag", required_argument, NULL, 10},		
    };
    if (1 == argc) { print_usage(stderr); print_time = 0; goto fail; }
    while ((c = getopt_long(argc, argv, "hs:S:O:o:R:b:i:I:p:", lopts, NULL)) != -1) {
        switch (c) {
            case 'h': print_usage(stderr); print_time = 0; goto fail;
            case 's': 
                    if (gs.in_fns) { str_arr_destroy(gs.in_fns, gs.nin); }
                    if (NULL == (gs.in_fns = hts_readlist(optarg, 0, &gs.nin)) || gs.nin <= 0) {
                        fprintf(stderr, "[E::%s] could not read input-list '%s' or list empty.\n", __func__, optarg);
                        goto fail;
                    } else { break;	}
            case 'S': 
                    if (gs.in_fn_file) free(gs.in_fn_file);
                    gs.in_fn_file = strdup(optarg); break;
            case 'O': 
                    if (gs.out_dir) { free(gs.out_dir); }
                    gs.out_dir = strdup(optarg); break;
            case 'o':
                    if (gs.out_fn) free(gs.out_fn); 
                    gs.out_fn = strdup(optarg); break;
            case 'R': 
                    if (gs.pos_list_file) free(gs.pos_list_file);
                    gs.pos_list_file = strdup(optarg); break;
            case 'b': 
                    if (gs.barcode_file) free(gs.barcode_file);
                    gs.barcode_file = strdup(optarg); break;
            case 'i':
                    if (gs.sid_list_file) free(gs.sid_list_file);
                    gs.sid_list_file = strdup(optarg); break;
            case 'I': 
                    if (gs.sample_ids) { str_arr_destroy(gs.sample_ids, gs.nsid); }
                    if (NULL == (gs.sample_ids = hts_readlist(optarg, 0, &gs.nsid))) {
                        fprintf(stderr, "[E::%s] could not read sample-id file '%s'\n", __func__, optarg);
                        goto fail;
                    } else { break; }
            case 'p': gs.nthread = atoi(optarg); break;
            case 1:  
                    if (gs.chrom_all) { str_arr_destroy(gs.chrom_all, gs.nchrom); }
                    if (NULL == (gs.chrom_all = hts_readlist(optarg, 0, &gs.nchrom))) {
                        fprintf(stderr, "[E::%s] could not read chrom-list '%s'\n", __func__, optarg);
                        goto fail;
                    }  else { break; }
            case 2:  
                    if (gs.cell_tag) free(gs.cell_tag);
                    gs.cell_tag = strdup(optarg); break;
            case 3:  
                    if (gs.umi_tag) free(gs.umi_tag);
                    gs.umi_tag = strdup(optarg); break;
            case 4:  gs.min_count = atoi(optarg); break;
            case 5:  gs.min_maf = atof(optarg); break;
            case 6:  gs.double_gl = 1; break;
            case 7:  gs.save_hdf5 = 1; break;
            case 8:  gs.min_len = atoi(optarg); break;
            case 9:  gs.min_mapq = atoi(optarg); break;
            case 10: gs.max_flag = atoi(optarg); break;
            default:  fprintf(stderr,"Invalid option: '%c'\n", c); goto fail;													
        }
    }
    fprintf(stderr, "[I::%s] start time: %s\n", __func__, time_str);
#if DEBUG
    fprintf(stderr, "[D::%s] global settings before checking:\n", __func__);
    gll_setting_print(stderr, &gs, "\t");
#endif
    /* check global settings */
    if ((ret = check_global_args(&gs)) < 0) { 
        fprintf(stderr, "[E::%s] error global settings\n", __func__);
        if (ret == -1) print_usage(stderr);
        goto fail;
    }
#if DEBUG
    fprintf(stderr, "[D::%s] global settings after checking:\n", __func__);
    gll_setting_print(stderr, &gs, "\t");
#endif
    /* prepare running data & options for each thread based on the checked global parameters.*/
    if (gs.nthread > 1 && NULL == (gs.tp = thpool_init(gs.nthread))) {
        fprintf(stderr, "[E::%s] could not initialize the thread pool.\n", __func__);
        goto fail;
    }
    /* run based on the mode of input. 
        Mode1: pileup a list of SNPs for a single BAM/SAM file with barcodes.
        Mode2: pileup whole chromosome(s) for one or multiple BAM/SAM files
        Mode3: pileup a list of SNPs for one or multiple BAM/SAM files with sample IDs.
    */
    if (gs.pos_list_file) {
        fprintf(stderr, "[I::%s] loading the VCF file for given SNPs ...\n", __func__);
        if (get_snplist(gs.pos_list_file, &gs.pl, &ret) <= 0 || ret < 0) {
            fprintf(stderr, "[E::%s] get SNP list from '%s' failed.\n", __func__, gs.pos_list_file);
            goto fail;
        }
        if (gs.barcodes) { 
            fprintf(stderr, "[I::%s] mode 1: fetch given SNPs in %d single cells.\n", __func__, gs.nbarcode); 
            if (run_mode1(&gs) < 0) { fprintf(stderr, "[E::%s] running mode 1 failed.\n", __func__); goto fail; } 
        } else { 
            fprintf(stderr, "[I::%s] mode 3: fetch given SNPs in %d bulk samples.\n", __func__, gs.nsid);
            if (run_mode3(&gs) < 0) { fprintf(stderr, "[E::%s] running mode 3 failed.\n", __func__); goto fail; } 
        }
    } else if (gs.chrom_all) { 
        if (gs.barcodes) { fprintf(stderr, "[I::%s] mode2: pileup %d whole chromosomes in %d single cells.\n", __func__, gs.nchrom, gs.nbarcode); }
        else { fprintf(stderr, "[I::%s] mode2: pileup %d whole chromosomes in one bulk sample.\n", __func__, gs.nchrom); }
        if (run_mode2(&gs) < 0) { fprintf(stderr, "[E::%s] running mode 2 failed.\n", __func__); goto fail; }
    } else {
        fprintf(stderr, "[E::%s] no proper mode to run, check input options.\n", __func__);
        print_usage(stderr);
        goto fail;
    }
    /* clean */
    gll_setting_free(&gs);
    /* calc time spent */
    if (print_time) {
        time(&end_time);
        time_info = localtime(&end_time);
        strftime(time_str, 30, "%Y-%m-%d %H:%M:%S", time_info);
        fprintf(stderr, "[I::%s] end time: %s\n", __func__, time_str);
        fprintf(stderr, "[I::%s] time spent: %ld seconds.\n", __func__, end_time - start_time);
    }
    return 0;
  fail:
    gll_setting_free(&gs);
    if (print_time) {
        time(&end_time);
        time_info = localtime(&end_time);
        strftime(time_str, 30, "%Y-%m-%d %H:%M:%S", time_info);
        fprintf(stderr, "[I::%s] end time: %s\n", __func__, time_str);
        fprintf(stderr, "[I::%s] time spent: %ld seconds.\n", __func__, end_time - start_time);
    }
    return 1;
}