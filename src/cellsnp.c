/* cellsnp main function
 * Author: Xianejie Huang <hxj5@hku.hk>
 */

/* TODO: 
- Try multi-process (process pool) for multi input samples
- Output vcf header according to input bam header
- merge csp_fetch() and csp_pileup() for the two functions share most codes?
- Try using multi_iter fetching method of bam/sam/cram for multi regions (SNPs) if it can in theory speed cellsnp up.
- Write test scripts for some key functions.
- Consistency correction could be done in UMI groups with the help of @p pu & @p pl inside mplp structure.
- More filters could be applied when extracting/fetching reads.
- Improve the jfile_t structure, for example, adding @p is_error.
- Improve the SZ_POOL structure, for example, adding @p base_init_f.
- Use optional sparse matrices tags with the help of function pointers.
- Output optionally qual values/letters to mtx file.
- Deal with the problem that some UMIs have the letter 'N'.
 */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <time.h>
#include "thpool.h"
#include "htslib/sam.h"
#include "config.h"
#include "csp.h"
#include "jfile.h"
#include "jsam.h"
#include "jstring.h"
#include "mplp.h"
#include "snp.h"

/*@abstract    Set default values for global_settings structure.
@param gs      Pointer to global_settings structure returned by gll_setting_init().
@return        Void.

@note          Internal use only!
 */
static void gll_set_default(global_settings *gs) {
    if (gs) {
        gs->in_fn_file = NULL; gs->in_fns = NULL; gs->nin = 0;
        gs->out_dir = NULL; 
        gs->out_vcf_base = NULL; gs->out_vcf_cells = NULL; gs->out_samples = NULL;
        gs->out_mtx_ad = NULL; gs->out_mtx_dp = NULL; gs->out_mtx_oth = NULL;
        gs->is_genotype = 0; gs->is_out_zip = 0;
        gs->snp_list_file = NULL; csp_snplist_init(gs->pl);
        gs->barcode_file = NULL; gs->nbarcode = 0; gs->barcodes = NULL; 
        gs->sid_list_file = NULL; gs->sample_ids = NULL; gs->nsid = 0;
        char *chrom_tmp[] = CSP_CHROM_ALL;
        gs->chroms = (char**) calloc(CSP_NCHROM, sizeof(char*));
        for (gs->nchrom = 0; gs->nchrom < CSP_NCHROM; gs->nchrom++) { gs->chroms[gs->nchrom] = safe_strdup(chrom_tmp[gs->nchrom]); }
        gs->cell_tag = safe_strdup(CSP_CELL_TAG); gs->umi_tag = safe_strdup(CSP_UMI_TAG);
        gs->nthread = CSP_NTHREAD; gs->tp = NULL;
        gs->min_count = CSP_MIN_COUNT; gs->min_maf = CSP_MIN_MAF; 
        gs->double_gl = 0;
        gs->min_len = CSP_MIN_LEN; gs->min_mapq = CSP_MIN_MAPQ;
        //gs->max_flag = -1;
        gs->rflag_filter = -1; gs->rflag_require = CSP_INCL_FMASK;
        gs->plp_max_depth = CSP_PLP_MAX_DEPTH; gs->no_orphan = CSP_NO_ORPHAN;
    }
}

static inline int run_mode1(global_settings *gs) { return csp_fetch(gs); }
static inline int run_mode2(global_settings *gs) { return csp_pileup(gs); }
static inline int run_mode3(global_settings *gs) { return csp_fetch(gs); }

static void print_usage(FILE *fp) {
    char *tmp_require = bam_flag2str(CSP_INCL_FMASK);
    char *tmp_filter_umi  = bam_flag2str(CSP_EXCL_FMASK_UMI);
    char *tmp_filter_noumi = bam_flag2str(CSP_EXCL_FMASK_NOUMI);

    fprintf(fp, 
"\n"
"Usage: %s [options]\n", CSP_NAME);
    fprintf(fp,
"\n"
"Options:\n"
"  -s, --samFile STR    Indexed sam/bam file(s), comma separated multiple samples.\n"
"                       Mode 1&2: one sam/bam file with single cell.\n"
"                       Mode 3: one or multiple bulk sam/bam files,\n"
"                       no barcodes needed, but sample ids and regionsVCF.\n"
"  -S, --samFileList FILE   A list file containing bam files, each per line, for Mode 3.\n"
"  -O, --outDir DIR         Output directory for VCF and sparse matrices.\n"
"  -R, --regionsVCF FILE    A vcf file listing all candidate SNPs, for fetch each variants.\n" 
"                           If None, pileup the genome. Needed for bulk samples.\n"
"  -b, --barcodeFile FILE   A plain file listing all effective cell barcode.\n"
"  -i, --sampleList FILE    A list file containing sample IDs, each per line.\n"
"  -I, --sampleIDs STR      Comma separated sample ids.\n"
"  -V, --version            Print software version and exit.\n"
"  -h, --help               Show this help message and exit.\n");
    fprintf(fp,
"\n"
"Optional arguments:\n"
"  --genotype           If use, do genotyping in addition to counting.\n"
"  --gzip               If use, the output files will be zipped into BGZF format.\n"
"  --printSkipSNPs      If use, the SNPs skipped when loading VCF will be printed.\n"
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
"\n"
"Read filtering:\n");
    fprintf(fp,
"  --inclFLAG STR|INT   Required flags: skip reads with all mask bits unset [%s]\n", tmp_require);
    fprintf(fp,
"  --exclFLAG STR|INT   Filter flags: skip reads with any mask bits set [%s\n"
"                       (when use UMI) or %s (otherwise)]\n", tmp_filter_umi, tmp_filter_noumi);
    fprintf(fp,
"  --minLEN INT         Minimum mapped length for read filtering [%d]\n", CSP_MIN_LEN);
    fprintf(fp,
"  --minMAPQ INT        Minimum MAPQ for read filtering [%d]\n", CSP_MIN_MAPQ);
    fprintf(fp,
"  --countORPHAN        If use, do not skip anomalous read pairs.\n");
/*
    fprintf(fp,
"  --maxFLAG INT        Maximum FLAG for read filtering [%d (when use UMI) or %d (otherwise)]\n", \
                        CSP_MAX_FLAG_WITH_UMI, CSP_MAX_FLAG_WITHOUT_UMI);
*/
    fputs("\n"
"Note that the \"--maxFLAG\" option is now deprecated, please use \"--inclFLAG\" or \"--exclFLAG\" instead.\n"
"You can easily aggregate and convert the flag mask bits to an integer by refering to:\n"
"https://broadinstitute.github.io/picard/explain-flags.html\n", fp);
    fputc('\n', fp);

    free(tmp_require); 
    free(tmp_filter_umi); 
    free(tmp_filter_noumi);
}

static inline int cmp_barcodes(const void *x, const void *y) {
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
    if (gs->out_dir) {
        if (0 != access(gs->out_dir, F_OK) && 0 != mkdir(gs->out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) { 
            fprintf(stderr, "[E::%s] '%s' does not exist.\n", __func__, gs->out_dir); 
            return -2; 
        }
    } else { fprintf(stderr, "[E::%s] should specify -O/--outDir option.\n", __func__); return -1; }
     /* 1. In current version, one and only one of barcodes and sample-ids would exist and work. Prefer barcodes. 
        2. For barcodes, the barcode file would not be read unless cell-tag is set, i.e. the barcodes and cell-tag are
           effective only when both of them are valid. 
        3. Codes below are a little repetitive and redundant, but it works well, maybe improve them in future.
    */
    if (gs->cell_tag && (0 == strcmp(gs->cell_tag, "None") || 0 == strcmp(gs->cell_tag, "none"))) { 
        free(gs->cell_tag);  gs->cell_tag = NULL; 
    }
    if (gs->sample_ids || gs->sid_list_file) {
        if (gs->barcode_file) { fprintf(stderr, "[E::%s] should not specify barcodes and sample IDs at the same time.\n", __func__); return -1; }
        else if (gs->cell_tag) { free(gs->cell_tag); gs->cell_tag = NULL; } 
    }
    if (gs->cell_tag && gs->barcode_file) {
        if (gs->sample_ids || gs->sid_list_file) { 
            fprintf(stderr, "[E::%s] should not specify barcodes and sample IDs at the same time.\n", __func__); 
            return -1; 
        } else if (NULL == (gs->barcodes = hts_readlines(gs->barcode_file, &gs->nbarcode))) {
            fprintf(stderr, "[E::%s] could not read barcode file '%s'\n", __func__, gs->barcode_file); 
            return -2;
        } else { qsort(gs->barcodes, gs->nbarcode, sizeof(char*), cmp_barcodes); }
    } else if ((NULL == gs->cell_tag) ^ (NULL == gs->barcode_file)) {
        fprintf(stderr, "[E::%s] should not specify barcodes or cell-tag alone.\n", __func__); 
        return -1;
    } else {
        if (NULL == gs->sample_ids) {
            if (NULL == gs->sid_list_file) { 
                if (NULL == (gs->sample_ids = (char**) calloc(gs->nin, sizeof(char*)))) { 
                    fprintf(stderr, "[E::%s] failed to allocate space for sample_ids\n", __func__); 
                    return -2; 
                }
                kstring_t ks = KS_INITIALIZE, *s = &ks;
                for (i = 0; i < gs->nin; i++) { ksprintf(s, "Sample_%d", i); gs->sample_ids[i] = strdup(ks_str(s)); ks_clear(s); }
                gs->nsid = i; ks_free(s);
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
       2. Sometimes, snp_list_file and chroms are both not NULL as the chroms has been set to default value when
          global_settings structure was just created. In this case, free chroms and save snp_list_file. */
    if (NULL == gs->snp_list_file || 0 == strcmp(gs->snp_list_file, "None") || 0 == strcmp(gs->snp_list_file, "none")) { 
        if (NULL == gs->chroms) { fprintf(stderr, "[E::%s] should specify -R/--regionsVCF or --chrom option.\n", __func__); return -1; }
        if (gs->snp_list_file) { free(gs->snp_list_file); gs->snp_list_file = NULL; }
    } else if (gs->chroms) { str_arr_destroy(gs->chroms, gs->nchrom); gs->chroms = NULL; gs->nchrom = 0; }
    if (gs->umi_tag) {
        if (0 == strcmp(gs->umi_tag, "Auto")) {
            if (gs->barcodes) { free(gs->umi_tag); gs->umi_tag = strdup("UR"); }
            else { free(gs->umi_tag); gs->umi_tag = NULL; }
        } else if (0 == strcmp(gs->umi_tag, "None") || 0 == strcmp(gs->umi_tag, "none")) { free(gs->umi_tag); gs->umi_tag = NULL; }
    }
    //if (gs->max_flag < 0) { gs->max_flag = gs->umi_tag ? CSP_MAX_FLAG_WITH_UMI : CSP_MAX_FLAG_WITHOUT_UMI; }
    if (gs->rflag_filter < 0) gs->rflag_filter = use_umi(gs) ? CSP_EXCL_FMASK_UMI : CSP_EXCL_FMASK_NOUMI;
    return 0;
}

/*@abstract    Output headers to files (vcf, mtx etc.)
@param fs      Pointer of jfile_t that the header will be writen into.
@param fm      File mode; if NULL, use default file mode inside jfile_t.
@param header  Header to be outputed, ends with '\0'.
@param len     Size of header.
@return        0 if success, negative numbers otherwise:
                 -1, open error; -2, write error; -3, close error.
 */
static inline int output_headers(jfile_t *fs, char *fm, char *header, size_t len) {
    int ret;
    if (jf_open(fs, fm) <= 0) { return -1; }
    if (jf_puts(header, fs) != len) { ret = -2; goto fail; }
    if (jf_close(fs) < 0) { ret = -3; goto fail; }
    return 0;
  fail:
    if (jf_isopen(fs)) { jf_close(fs); }
    return ret; 
}

static inline char* format_fn(char *fn, int is_zip, kstring_t *s) {
    if (is_zip) {
        kputs(fn, s); kputs(".gz", s);
        return strdup(ks_str(s));
    } else { return fn; }
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
    kstring_t ks = KS_INITIALIZE, *s = &ks;
    int c, k, ret, print_time = 0, print_skip_snp = 0;
    struct option lopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {"samFile", required_argument, NULL, 's'},
        {"samfile", required_argument, NULL, 's'},
        {"samFileList", required_argument, NULL, 'S'},			
        {"outDir", required_argument, NULL, 'O'},
        {"outdir", required_argument, NULL, 'O'},
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
        {"minLEN", required_argument, NULL, 8},
        {"minLen", required_argument, NULL, 8},
        {"minlen", required_argument, NULL, 8},
        {"minMAPQ", required_argument, NULL, 9},
        //{"maxFLAG", required_argument, NULL, 10},
        //{"maxFlag", required_argument, NULL, 10},
        //{"maxflag", required_argument, NULL, 10},
        {"genotype", no_argument, NULL, 11},
        {"gzip", no_argument, NULL, 12},
        {"printSkipSNPs", no_argument, NULL, 13},
        {"inclFLAG", required_argument, NULL, 14},
        {"exclFLAG", required_argument, NULL, 15},
        {"countORPHAN", no_argument, NULL, 16}
    };
    if (1 == argc) { print_usage(stderr); goto fail; }
    while ((c = getopt_long(argc, argv, "hVs:S:O:R:b:i:I:p:", lopts, NULL)) != -1) {
        switch (c) {
            case 'h': print_usage(stderr); goto fail;
            case 'V': printf("%s\n", CSP_VERSION); goto fail;
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
            case 'R': 
                    if (gs.snp_list_file) free(gs.snp_list_file);
                    gs.snp_list_file = strdup(optarg); break;
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
                    if (gs.chroms) { str_arr_destroy(gs.chroms, gs.nchrom); gs.nchrom = 0; }
                    if (NULL == (gs.chroms = hts_readlist(optarg, 0, &gs.nchrom))) {
                        fprintf(stderr, "[E::%s] could not read chrom-list '%s'\n", __func__, optarg);
                        goto fail;
                    }  else { break; }
            case 2:  
                    if (gs.cell_tag) free(gs.cell_tag);
                    gs.cell_tag = strdup(optarg); break;
            case 3:  
                    if (gs.umi_tag) free(gs.umi_tag);
                    gs.umi_tag = strdup(optarg); break;
            case 4: gs.min_count = atoi(optarg); break;
            case 5: gs.min_maf = atof(optarg); break;
            case 6: gs.double_gl = 1; break;
            case 8: gs.min_len = atoi(optarg); break;
            case 9: gs.min_mapq = atoi(optarg); break;
            //case 10: gs.max_flag = atoi(optarg); break;
            case 11: gs.is_genotype = 1; break;
            case 12: gs.is_out_zip = 1; break;
            case 13: print_skip_snp = 1; break;
            case 14: 
                    if ((gs.rflag_require = bam_str2flag(optarg)) < 0) {
                        fprintf(stderr, "[E::%s] could not parse --inclFLAG '%s'\n", __func__, optarg);
                        goto fail;
                    } else { break; }
            case 15:
                    if ((gs.rflag_filter = bam_str2flag(optarg)) < 0) {
                        fprintf(stderr, "[E::%s] could not parse --exclFLAG '%s'\n", __func__, optarg);
                        goto fail;
                    } else { break; }
            case 16: gs.no_orphan = 0; break;
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
        if (ret == -1) { print_usage(stderr); }
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
    /* prepare output files. */
    if (NULL == (gs.out_mtx_ad = jf_init()) || NULL == (gs.out_mtx_dp = jf_init()) || \
        NULL == (gs.out_mtx_oth = jf_init()) || NULL == (gs.out_samples = jf_init()) || \
        NULL == (gs.out_vcf_base = jf_init()) || (gs.is_genotype && NULL == (gs.out_vcf_cells = jf_init()))) {
        fprintf(stderr, "[E::%s] fail to create jfile_t.\n", __func__);
        goto fail;
    }
    gs.out_mtx_ad->is_zip = 0; gs.out_mtx_ad->is_tmp = 0;
    gs.out_mtx_ad->fn = format_fn(join_path(gs.out_dir, CSP_OUT_MTX_AD), gs.out_mtx_ad->is_zip, s); ks_clear(s);
    gs.out_mtx_dp->is_zip = 0; gs.out_mtx_dp->is_tmp = 0;
    gs.out_mtx_dp->fn = format_fn(join_path(gs.out_dir, CSP_OUT_MTX_DP), gs.out_mtx_dp->is_zip, s); ks_clear(s); 
    gs.out_mtx_oth->is_zip = 0; gs.out_mtx_oth->is_tmp = 0;
    gs.out_mtx_oth->fn = format_fn(join_path(gs.out_dir, CSP_OUT_MTX_OTH), gs.out_mtx_oth->is_zip, s); ks_clear(s);
    gs.out_vcf_base->is_zip = gs.is_out_zip; gs.out_vcf_base->is_tmp = 0;
    gs.out_vcf_base->fn = format_fn(join_path(gs.out_dir, CSP_OUT_VCF_BASE), gs.out_vcf_base->is_zip, s); ks_clear(s);
    gs.out_samples->is_zip = 0; gs.out_samples->is_tmp = 0;
    gs.out_samples->fn = format_fn(join_path(gs.out_dir, CSP_OUT_SAMPLES), gs.out_samples->is_zip, s); ks_clear(s);
    if (gs.is_genotype) { 
        gs.out_vcf_cells->is_zip = gs.is_out_zip; gs.out_vcf_cells->is_tmp = 0;
        gs.out_vcf_cells->fn = format_fn(join_path(gs.out_dir, CSP_OUT_VCF_CELLS), gs.out_vcf_cells->is_zip, s); ks_clear(s);
    } // no need to set is_tmp for these out files.
    /* output headers to files. */
    kputs(CSP_MTX_HEADER, s);
    if (output_headers(gs.out_mtx_ad, "wb", ks_str(s), ks_len(s)) < 0) {   // output header to mtx_AD
        fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_mtx_ad->fn);
        goto fail;
    }
    if (output_headers(gs.out_mtx_dp, "wb", ks_str(s), ks_len(s)) < 0) {   // output header to mtx_DP
        fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_mtx_dp->fn);
        goto fail;
    }
    if (output_headers(gs.out_mtx_oth, "wb", ks_str(s), ks_len(s)) < 0) {  // output header to mtx_OTH
        fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_mtx_oth->fn);
        goto fail;
    } ks_clear(s);
    if (use_barcodes(&gs)) {                     // output samples.
        for (k = 0; k < gs.nbarcode; k++) { kputs(gs.barcodes[k], s); kputc('\n', s); }
    } else if (use_sid(&gs)) {
        for (k = 0; k < gs.nsid; k++) { kputs(gs.sample_ids[k], s); kputc('\n', s); }
    } // else: should not come here!
    if (output_headers(gs.out_samples, "wb", ks_str(s), ks_len(s)) < 0) {
        fprintf(stderr, "[E::%s] fail to write samples to '%s'\n", __func__, gs.out_samples->fn);
        goto fail;
    } ks_clear(s);
    kputs(CSP_VCF_BASE_HEADER, s);             // output header to vcf base.
    kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", s);
    if (output_headers(gs.out_vcf_base, "wb", ks_str(s), ks_len(s)) < 0) {
        fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_vcf_base->fn);
        goto fail;
    } ks_clear(s);
    if (gs.is_genotype) {
        kputs(CSP_VCF_CELLS_HEADER, s);           // output header to vcf cells.
        kputs(CSP_VCF_CELLS_CONTIG, s);
        kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", s);
        if (use_barcodes(&gs) && gs.barcodes) {
            for (k = 0; k < gs.nbarcode; k++) { kputc_('\t', s); kputs(gs.barcodes[k], s); }
        } else if (use_sid(&gs) && gs.sample_ids) {
            for (k = 0; k < gs.nsid; k++) { kputc_('\t', s); kputs(gs.sample_ids[k], s); }
        } else { fprintf(stderr, "[E::%s] neither barcodes or sample IDs exist.\n", __func__); goto fail; }
        kputc('\n', s);
        if (output_headers(gs.out_vcf_cells, "wb", ks_str(s), ks_len(s)) < 0) {
            fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_vcf_cells->fn);
            goto fail;
        }
    }
    /* set file modes. */
    gs.out_mtx_ad->fm = gs.out_mtx_dp->fm = gs.out_mtx_oth->fm = "ab";
    gs.out_vcf_base->fm = "ab";
    if (gs.is_genotype) { gs.out_vcf_cells->fm = "ab"; }
    /* run based on the mode of input. 
        Mode1: pileup a list of SNPs for a single BAM/SAM file with barcodes.
        Mode2: pileup whole chromosome(s) for one or multiple BAM/SAM files
        Mode3: pileup a list of SNPs for one or multiple BAM/SAM files with sample IDs.
    */
    if (gs.snp_list_file) {
        fprintf(stderr, "[I::%s] loading the VCF file for given SNPs ...\n", __func__);
        if (get_snplist(gs.snp_list_file, &gs.pl, &ret, print_skip_snp) <= 0 || ret < 0) {
            fprintf(stderr, "[E::%s] get SNP list from '%s' failed.\n", __func__, gs.snp_list_file);
            print_time = 1; goto fail;
        } else { fprintf(stderr, "[I::%s] fetching %ld candidate variants ...\n", __func__, csp_snplist_size(gs.pl)); }
        if (gs.barcodes) { 
            fprintf(stderr, "[I::%s] mode 1: fetch given SNPs in %d single cells.\n", __func__, gs.nbarcode); 
            if (run_mode1(&gs) < 0) { fprintf(stderr, "[E::%s] running mode 1 failed.\n", __func__); print_time = 1; goto fail; } 
        } else { 
            fprintf(stderr, "[I::%s] mode 3: fetch given SNPs in %d bulk samples.\n", __func__, gs.nsid);
            if (run_mode3(&gs) < 0) { fprintf(stderr, "[E::%s] running mode 3 failed.\n", __func__); print_time = 1; goto fail; } 
        }
    } else if (gs.chroms) { 
        if (gs.barcodes) { fprintf(stderr, "[I::%s] mode2: pileup %d whole chromosomes in %d single cells.\n", __func__, gs.nchrom, gs.nbarcode); }
        else { fprintf(stderr, "[I::%s] mode2: pileup %d whole chromosomes in one bulk sample.\n", __func__, gs.nchrom); }
        if (run_mode2(&gs) < 0) { fprintf(stderr, "[E::%s] running mode 2 failed.\n", __func__); print_time = 1; goto fail; }
    } else {
        fprintf(stderr, "[E::%s] no proper mode to run, check input options.\n", __func__);
        print_usage(stderr);
        goto fail;
    }
    /* clean */
    ks_free(s); s = NULL;
    gll_setting_free(&gs);
    fprintf(stderr, "[I::%s] All Done!\n", __func__);
    /* calc time spent */
    time(&end_time);
    time_info = localtime(&end_time);
    strftime(time_str, 30, "%Y-%m-%d %H:%M:%S", time_info);
    fprintf(stderr, "[I::%s] end time: %s\n", __func__, time_str);
    fprintf(stderr, "[I::%s] time spent: %ld seconds.\n", __func__, end_time - start_time);
    return 0;
  fail:
    if (s) { ks_free(s); }
    gll_setting_free(&gs);
    if (print_time) {
        fprintf(stderr, "[E::%s] Quiting...\n", __func__);
        time(&end_time);
        time_info = localtime(&end_time);
        strftime(time_str, 30, "%Y-%m-%d %H:%M:%S", time_info);
        fprintf(stderr, "[I::%s] end time: %s\n", __func__, time_str);
        fprintf(stderr, "[I::%s] time spent: %ld seconds.\n", __func__, end_time - start_time);
    }
    return 1;
}
