/* cellsnp.c - cellsnp main function
 * Author: Xianejie Huang <hxj5@hku.hk>
 */

/* TODO: 
- Add --max-base-qual option to set max per base quality for genotyping
  * in htslib, per base qualities are stored in the Phred scale with no +33 offset.
  * get_qual_vector() in mplp.c.
- Write test scripts for some key functions.
- Consistency correction could be done in UMI groups with the help of @p pu & @p pl inside mplp structure.
  * update map_ug_t first!
- Output vcf header according to input bam header
- More filters could be applied when extracting/fetching reads.
- Add fetch and pileup sub-commands?
  * merge csp_fetch() and csp_pileup() for the two functions share most codes?
- Support calling germline SNPs for multiple bam files?
  * in bulk mode
- ?Change -T method, use qsort & linear search instead of regidx_t of htslib:
  * note the bug of qsort in lower version of glibc, refer to https://sourceware.org/bugzilla/show_bug.cgi?id=11655)
  * as the linear searching assume that the bam has been sortted by start pos, then how to deal with the partly aligned reads
    and one read in paired reads aligned to query chrom.
- Save memory for snp list
  * use fixed 32-bit index for chrom string.
  * use shared chrom string for all snps in that chrom, like the way in the develop branch.
- Use optional sparse matrices tags with the help of function pointers.
- Wrap up the steps in the whole workflow into proper data structures.
  * e.g., the snp iteration, the read extraction from input files (especially the redundant hts_itr_xxx), 
    mplp statistics, temporary file creation & output, etc. 
- Update mplp_t::su and plp_t::hug
- Try multi-process (process pool) for multi input samples
- Separate htsFile from csp_bam_fs as it cannot be shared among threads
- Try using multi_iter fetching method of bam/sam/cram for multi regions (SNPs) if it can in theory speed cellsnp up.
- Improve the jfile_t structure, for example, adding @p is_error.
- Improve the JMEMPOOL structure, for example, adding @p base_init_f.
- Output optionally qual values/letters to mtx file.
- Deal with the problem that some UMIs have the letter 'N'.
 */

/* Discuss
- Cellsnp-lite now only considers biallele (related to -f/--refseq).
  * !! WARNING !! in mode 1, the output could contain Homozygous SNV even with --minMAF 0.3. e.g., assuming one input SNV has
    REF/ALT - A/C, while the two real alleles are C/G with AF 0.6/0.4, then this SNV would pass --minMAF 0.3
    and the genotype is 1/1 (as REF is A, ALT is C), while the real genotype should be 1/2 (as two alt alleles C,G).
    (SNVs of this kind are not so many in practice? - in a recent case, only 158 out of 133k SNVs)
 */

/* Reference
- htslib header files: https://github.com/samtools/htslib/tree/develop/htslib
  mainly the sam.{h,c}, regidx.{h,c}, vcf.{h,c}, kstring.h, hts.{h,c} files
- mpileup.c in bcftools: https://github.com/samtools/bcftools/blob/develop/mpileup.c
- bam_plcmd.c in samtools: https://github.com/samtools/samtools/blob/develop/bam_plcmd.c
  refer to the cmdline options in this file too.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#include <libgen.h>
#include <time.h>
#include "thpool.h"
#include "htslib/sam.h"
#include "htslib/regidx.h"
#include "htslib/hts.h"
#include "config.h"
#include "csp.h"
#include "jfile.h"
#include "jsam.h"
#include "jstring.h"
#include "mplp.h"
#include "snp.h"

// Set default values for global_settings structure.
// Internal use only!
static void gll_set_default(global_settings *gs) {
    if (gs) {
        // Input Output
        gs->in_fn_file = NULL; gs->in_fns = NULL; gs->nin = 0;
        gs->out_dir = NULL; 
        gs->out_vcf_base = NULL; gs->out_vcf_cells = NULL; gs->out_samples = NULL;
        gs->out_mtx_ad = NULL; gs->out_mtx_dp = NULL; gs->out_mtx_oth = NULL;
        gs->snp_list_file = NULL; kv_init(gs->pl); gs->is_target = 0; gs->targets = NULL;
        gs->barcode_file = NULL; gs->nbarcode = 0; gs->barcodes = NULL; 
        gs->sid_list_file = NULL; gs->sample_ids = NULL; gs->nsid = 0;

        // Options
        gs->is_genotype = 0; gs->is_out_zip = 0;
        gs->refseq_file = NULL;
        char *chrom_tmp[] = CSP_CHROM_ALL;
        gs->chroms = (char**) calloc(CSP_NCHROM, sizeof(char*));
        for (gs->nchrom = 0; gs->nchrom < CSP_NCHROM; gs->nchrom++) {
            gs->chroms[gs->nchrom] = safe_strdup(chrom_tmp[gs->nchrom]);
        }
        gs->cell_tag = safe_strdup(CSP_CELL_TAG); gs->umi_tag = safe_strdup(CSP_UMI_TAG);
        gs->nthread = CSP_NTHREAD; gs->tp = NULL; gs->tp_max_open = TP_MAX_OPEN;
        gs->mthread = CSP_NTHREAD; gs->tp_errno = 0; gs->tp_ntry = 0;
        gs->min_count = CSP_MIN_COUNT; gs->min_maf = CSP_MIN_MAF; 
        gs->doublet_gl = 0;

        // Read Filtering
        gs->min_len = CSP_MIN_LEN; gs->min_mapq = CSP_MIN_MAPQ;
        gs->rflag_filter = -1; gs->rflag_require = CSP_INCL_FMASK;
        gs->max_depth = CSP_MAX_DEPTH; gs->max_pileup = CSP_MAX_PILEUP;
        gs->no_orphan = CSP_NO_ORPHAN;
    }
}

static inline int run_mode1a(global_settings *gs)   { return csp_fetch(gs); }
static inline int run_mode1a_T(global_settings *gs) { return csp_pileup(gs); }
static inline int run_mode1b(global_settings *gs)   { return csp_fetch(gs); }
static inline int run_mode1b_T(global_settings *gs) { return csp_pileup(gs); }
static inline int run_mode2a(global_settings *gs)   { return csp_pileup(gs); }
static inline int run_mode2b(global_settings *gs)   { return csp_pileup(gs); }

static void print_usage(FILE *fp) {
    char *tmp_require = bam_flag2str(CSP_INCL_FMASK);
    char *tmp_filter_umi  = bam_flag2str(CSP_EXCL_FMASK_UMI);
    char *tmp_filter_noumi = bam_flag2str(CSP_EXCL_FMASK_NOUMI);

    fprintf(fp, "\n");
    fprintf(fp, "Version: %s (htslib %s)\n", CSP_VERSION, hts_version());
    fprintf(fp, "Usage:   %s [options]\n", CSP_NAME);
    fprintf(fp, "\n");
    fprintf(fp, "Options:\n");
    fprintf(fp, "  -s, --samFile STR    Indexed sam/bam file(s), comma separated multiple samples.\n");
    fprintf(fp, "                       Mode 1a & 2a: one sam/bam file with single cell.\n");
    fprintf(fp, "                       Mode 1b & 2b: one or multiple bulk sam/bam files,\n");
    fprintf(fp, "                       no barcodes needed, but sample ids and regionsVCF.\n");
    fprintf(fp, "  -S, --samFileList FILE   A list file containing bam files, each per line, for Mode 1b & 2b.\n");
    fprintf(fp, "  -O, --outDir DIR         Output directory for VCF and sparse matrices.\n");
    fprintf(fp, "  -R, --regionsVCF FILE    A vcf file listing all candidate SNPs, for fetch each variants.\n");
    fprintf(fp, "                           If None, pileup the genome. Needed for bulk samples.\n");
    fprintf(fp, "  -T, --targetsVCF FILE    Similar as -R, but the next position is accessed by streaming rather\n");
    fprintf(fp, "                           than indexing/jumping (like -T in samtools/bcftools mpileup).\n");
    fprintf(fp, "  -b, --barcodeFile FILE   A plain file listing all effective cell barcode.\n");
    fprintf(fp, "  -i, --sampleList FILE    A list file containing sample IDs, each per line.\n");
    fprintf(fp, "  -I, --sampleIDs STR      Comma separated sample ids.\n");
    fprintf(fp, "  -V, --version            Print software version and exit.\n");
    fprintf(fp, "  -h, --help               Show this help message and exit.\n");
    fprintf(fp, "\n");
    fprintf(fp, "Optional arguments:\n");
    fprintf(fp, "  --genotype           If use, do genotyping in addition to counting.\n");
    fprintf(fp, "  --gzip               If use, the output files will be zipped into BGZF format.\n");
    fprintf(fp, "  --printSkipSNPs      If use, the SNPs skipped when loading VCF will be printed.\n");
    fprintf(fp, "  -p, --nproc INT      Number of subprocesses [%d]\n", CSP_NTHREAD);
    fprintf(fp, "  -f, --refseq FILE    Faidx indexed reference sequence file. If set, the real (genomic)\n");
    fprintf(fp, "                       ref extracted from this file would be used for Mode 2 or for the\n");
    fprintf(fp, "                       missing REFs in the input VCF for Mode 1.\n");
    fprintf(fp, "  --chrom STR          The chromosomes to use, comma separated [1 to %d]\n", CSP_NCHROM);
    fprintf(fp, "  --cellTAG STR        Tag for cell barcodes, turn off with None [%s]\n", CSP_CELL_TAG);
    fprintf(fp, "  --UMItag STR         Tag for UMI: UB, Auto, None. For Auto mode, use UB if barcodes are inputted,\n");
    fprintf(fp, "                       otherwise use None. None mode means no UMI but read counts [%s]\n", CSP_UMI_TAG);
    fprintf(fp, "  --minCOUNT INT       Minimum aggragated count [%d]\n", CSP_MIN_COUNT);
    fprintf(fp, "  --minMAF FLOAT       Minimum minor allele frequency [%.2f]\n", CSP_MIN_MAF);
    fprintf(fp, "  --doubletGL          If use, keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5.\n");
    fprintf(fp, "\n");
    fprintf(fp, "Read filtering:\n");
    fprintf(fp, "  --inclFLAG STR|INT   Required flags: skip reads with all mask bits unset [%s]\n", tmp_require);
    fprintf(fp, "  --exclFLAG STR|INT   Filter flags: skip reads with any mask bits set [%s\n", tmp_filter_umi);
    fprintf(fp, "                       (when use UMI) or %s (otherwise)]\n", tmp_filter_noumi);
    fprintf(fp, "  --minLEN INT         Minimum mapped length for read filtering [%d]\n", CSP_MIN_LEN);
    fprintf(fp, "  --minMAPQ INT        Minimum MAPQ for read filtering [%d]\n", CSP_MIN_MAPQ);
    fprintf(fp, "  --maxPILEUP INT      Maximum pileup for one site of one file (including those filtered reads),\n");
    fprintf(fp, "                       avoids excessive memory usage; 0 means highest possible value [%d]\n", CSP_MAX_PILEUP);
    fprintf(fp, "  --maxDEPTH INT       Maximum depth for one site of one file (excluding those filtered reads),\n");
    fprintf(fp, "                       avoids excessive memory usage; 0 means highest possible value [%d]\n", CSP_MAX_DEPTH);
    fprintf(fp, "  --countORPHAN        If use, do not skip anomalous read pairs.\n");
    fprintf(fp, "\n");
    fprintf(fp, "Note that the \"--maxFLAG\" option is now deprecated, please use \"--inclFLAG\" or \"--exclFLAG\"\n");
    fprintf(fp, "instead. You can easily aggregate and convert the flag mask bits to an integer by refering to:\n");
    fprintf(fp, "https://broadinstitute.github.io/picard/explain-flags.html\n");
    fprintf(fp, "\n");

    free(tmp_require); 
    free(tmp_filter_umi); 
    free(tmp_filter_noumi);
}

static inline int cmp_barcodes(const void *x, const void *y) {
    return strcmp(*((char**) x), *((char**) y));
}

/*!@func
@abstract  Perform basic check for global settings right after running getopt()/getopt_long() function.
@param gs  Pointer to the global settings.
@return    0 if no error, negative numbers otherwise:
             -1, should print_usage after return.
             -2, no action.
*/
static int check_args(global_settings *gs) {
    int i;
    struct rlimit r;
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
        if (0 != access(gs->in_fns[i], F_OK)) {
            fprintf(stderr, "[E::%s] '%s' does not exist.\n", __func__, gs->in_fns[i]);
            return -2;
        }
    }

    if (gs->out_dir) {
        if (0 != access(gs->out_dir, F_OK) && 0 != mkdir(gs->out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) { 
            fprintf(stderr, "[E::%s] '%s' does not exist.\n", __func__, gs->out_dir); 
            return -2; 
        }
    } else {
        fprintf(stderr, "[E::%s] should specify -O/--outDir option.\n", __func__);
        return -1;
    }

     /* 1. In current version, one and only one of barcodes and sample-ids would exist and work. Prefer barcodes. 
        2. For barcodes, the barcode file would not be read unless cell-tag is set, i.e. the barcodes and cell-tag are
           effective only when both of them are valid. 
        3. Codes below are a little repetitive and redundant, but it works well, maybe improve them in future.
    */
    if (gs->cell_tag && (0 == strcmp(gs->cell_tag, "None") || 0 == strcmp(gs->cell_tag, "none"))) { 
        free(gs->cell_tag);  gs->cell_tag = NULL; 
    }
    if (gs->sample_ids || gs->sid_list_file) {
        if (gs->barcode_file) {
            fprintf(stderr, "[E::%s] should not specify barcodes and sample IDs at the same time.\n", __func__);
            return -1;
        } else if (gs->cell_tag) {
            free(gs->cell_tag); gs->cell_tag = NULL;
        } 
    }
    if (gs->cell_tag && gs->barcode_file) {
        if (gs->sample_ids || gs->sid_list_file) { 
            fprintf(stderr, "[E::%s] should not specify barcodes and sample IDs at the same time.\n", __func__); 
            return -1; 
        } else if (NULL == (gs->barcodes = hts_readlines(gs->barcode_file, &gs->nbarcode))) {
            fprintf(stderr, "[E::%s] could not read barcode file '%s'\n", __func__, gs->barcode_file); 
            return -2;
        } else {
            qsort(gs->barcodes, gs->nbarcode, sizeof(char*), cmp_barcodes);
        }
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
                for (i = 0; i < gs->nin; i++) {
                    ksprintf(s, "Sample_%d", i); 
                    gs->sample_ids[i] = strdup(ks_str(s));
                    ks_clear(s); 
                }
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
          global_settings structure was just created. In this case, free chroms and save snp_list_file.
     */
    if (NULL == gs->snp_list_file || 0 == strcmp(gs->snp_list_file, "None") || 0 == strcmp(gs->snp_list_file, "none")) { 
        if (NULL == gs->chroms) {
            fprintf(stderr, "[E::%s] should specify -R/--regionsVCF or --chrom option.\n", __func__);
            return -1;
        }
        if (gs->snp_list_file) {
            free(gs->snp_list_file); gs->snp_list_file = NULL;
        }
    } else if (gs->chroms) {
        str_arr_destroy(gs->chroms, gs->nchrom);
        gs->chroms = NULL; gs->nchrom = 0;
    }
    if (gs->umi_tag) {
        if (0 == strcmp(gs->umi_tag, "Auto")) {
            if (gs->barcodes) {
                free(gs->umi_tag); gs->umi_tag = strdup("UB"); 
            } else {
                free(gs->umi_tag); gs->umi_tag = NULL;
            }
        } else if (0 == strcmp(gs->umi_tag, "None") || 0 == strcmp(gs->umi_tag, "none")) {
            free(gs->umi_tag); gs->umi_tag = NULL;
        }
    }

    if (gs->rflag_filter < 0) 
        gs->rflag_filter = use_umi(gs) ? CSP_EXCL_FMASK_UMI : CSP_EXCL_FMASK_NOUMI;

    // increase number of max open files
    if (gs->nin > 1) {
        if (getrlimit(RLIMIT_NOFILE, &r) < 0) {
            fprintf(stderr, "[E::%s] getrlimit error.\n", __func__);
            return -2;
        }
        fprintf(stderr, "[I::%s] original limits of max open, soft = %d, hard = %d\n", __func__, (int) r.rlim_cur, (int) r.rlim_max);
        r.rlim_cur = r.rlim_max;
        if (setrlimit(RLIMIT_NOFILE, &r) < 0) {
            fprintf(stderr, "[E::%s] setrlimit error.\n", __func__);
            return -2;
        }
        if (getrlimit(RLIMIT_NOFILE, &r) < 0) {
            fprintf(stderr, "[E::%s] getrlimit error.\n", __func__);
            return -2;
        }
        fprintf(stderr, "[I::%s] new limits of max open, soft = %d, hard = %d\n", __func__, (int) r.rlim_cur, (int) r.rlim_max);
        gs->tp_max_open = (int) r.rlim_cur;
    } else {
        if (getrlimit(RLIMIT_NOFILE, &r) < 0) {
            fprintf(stderr, "[E::%s] getrlimit error.\n", __func__);
            return -2;
        }
        gs->tp_max_open = (int) r.rlim_cur;
    }

    // check max-depth
    if (gs->max_depth <= 0) {
        gs->max_depth = INT_MAX;
        fprintf(stderr, "[W::%s] Max depth set to maximum value (%d)\n", __func__, INT_MAX);
    }

    // check max-pileup
    if (gs->max_pileup <= 0) {
        gs->max_pileup = INT_MAX;
        fprintf(stderr, "[W::%s] Max pileup set to maximum value (%d)\n", __func__, INT_MAX);
    }

    return 0;
}

/*!@func
@abstract      Output headers to files (vcf, mtx etc.)
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

static inline void regidx_payload_free(void *payload) {
    biallele_t **data = (biallele_t**) payload;
    biallele_destroy(*data);
}

int main(int argc, char **argv) {
    int is_ok = 0;

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
    int i, c, k, ret, print_info = 0, print_skip_snp = 0;

    struct option lopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"version", no_argument, NULL, 'V'},
        {"samFile", required_argument, NULL, 's'},
        {"samfile", required_argument, NULL, 's'},
        {"samFileList", required_argument, NULL, 'S'},			
        {"samfilelist", required_argument, NULL, 'S'},			
        {"outDir", required_argument, NULL, 'O'},
        {"outdir", required_argument, NULL, 'O'},
        {"regionsVCF", required_argument, NULL, 'R'},
        {"regionsvcf", required_argument, NULL, 'R'},
        {"targetsVCF", required_argument, NULL, 'T'},
        {"targetsvcf", required_argument, NULL, 'T'},
        {"barcodeFile", required_argument, NULL, 'b'},
        {"barcodefile", required_argument, NULL, 'b'},
        {"sampleList", required_argument, NULL, 'i'},
        {"samplelist", required_argument, NULL, 'i'},
        {"sampleIDs", required_argument, NULL, 'I'},
        {"sampleids", required_argument, NULL, 'I'},
        {"nproc", required_argument, NULL, 'p'},
        {"refseq", required_argument, NULL, 'f'},
        {"chrom", required_argument, NULL, 1},
        {"cellTAG", required_argument, NULL, 2},
        {"celltag", required_argument, NULL, 2},
        {"UMItag", required_argument, NULL, 3},
        {"umitag", required_argument, NULL, 3},
        {"minCOUNT", required_argument, NULL, 4},
        {"minCount", required_argument, NULL, 4},
        {"mincount", required_argument, NULL, 4},
        {"minMAF", required_argument, NULL, 5},
        {"minMaf", required_argument, NULL, 5},
        {"minmaf", required_argument, NULL, 5},
        {"doubletGL", no_argument, NULL, 6},
        {"doubletGl", no_argument, NULL, 6},
        {"doubletgl", no_argument, NULL, 6},
        {"minLEN", required_argument, NULL, 8},
        {"minLen", required_argument, NULL, 8},
        {"minlen", required_argument, NULL, 8},
        {"minMAPQ", required_argument, NULL, 9},
        {"minMapq", required_argument, NULL, 9},
        {"minmapq", required_argument, NULL, 9},
        //{"maxFLAG", required_argument, NULL, 10},
        //{"maxFlag", required_argument, NULL, 10},
        //{"maxflag", required_argument, NULL, 10},
        {"genotype", no_argument, NULL, 11},
        {"gzip", no_argument, NULL, 12},
        {"printSkipSNPs", no_argument, NULL, 13},
        {"printskipsnps", no_argument, NULL, 13},
        {"inclFLAG", required_argument, NULL, 14},
        {"inclflag", required_argument, NULL, 14},
        {"exclFLAG", required_argument, NULL, 15},
        {"exclflag", required_argument, NULL, 15},
        {"countORPHAN", no_argument, NULL, 16},
        {"countorphan", no_argument, NULL, 16},
        {"maxDEPTH", required_argument, NULL, 17},
        {"maxDepth", required_argument, NULL, 17},
        {"maxdepth", required_argument, NULL, 17},
        {"maxPILEUP", required_argument, NULL, 18},
        {"maxPileup", required_argument, NULL, 18},
        {"maxpileup", required_argument, NULL, 18},
        {NULL, 0, NULL, 0}
    };
    if (1 == argc) { print_usage(stderr); goto clean; }
    while ((c = getopt_long(argc, argv, "hVs:S:O:R:T:b:i:I:p:f:", lopts, NULL)) != -1) {
        switch (c) {
            case 'h': print_usage(stderr); goto clean;
            case 'V': printf("%s %s (htslib %s)\n", CSP_NAME, CSP_VERSION, hts_version()); goto clean;
            case 's': 
                    if (gs.in_fns) { str_arr_destroy(gs.in_fns, gs.nin); }
                    if (NULL == (gs.in_fns = hts_readlist(optarg, 0, &gs.nin)) || gs.nin <= 0) {
                        fprintf(stderr, "[E::%s] could not read input-list '%s' or list empty.\n", __func__, optarg);
                        goto clean;
                    } else { break; }
            case 'S': 
                    if (gs.in_fn_file) free(gs.in_fn_file);
                    gs.in_fn_file = strdup(optarg); break;
            case 'O': 
                    if (gs.out_dir) { free(gs.out_dir); }
                    gs.out_dir = strdup(optarg); break;
            case 'R': 
                    if (gs.snp_list_file) free(gs.snp_list_file);
                    gs.snp_list_file = strdup(optarg); gs.is_target = 0; break;
            case 'T':
                    if (gs.snp_list_file) free(gs.snp_list_file);
                    gs.snp_list_file = strdup(optarg); gs.is_target = 1; break;
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
                        goto clean;
                    } else { break; }
            case 'p': gs.mthread = gs.nthread = atoi(optarg); break;
            case 'f':
                    if (gs.refseq_file) free(gs.refseq_file);
                    gs.refseq_file = strdup(optarg); break;
            case 1:  
                    if (gs.chroms) { str_arr_destroy(gs.chroms, gs.nchrom); gs.nchrom = 0; }
                    if (NULL == (gs.chroms = hts_readlist(optarg, 0, &gs.nchrom))) {
                        fprintf(stderr, "[E::%s] could not read chrom-list '%s'\n", __func__, optarg);
                        goto clean;
                    }  else { break; }
            case 2:  
                    if (gs.cell_tag) free(gs.cell_tag);
                    gs.cell_tag = strdup(optarg); break;
            case 3:  
                    if (gs.umi_tag) free(gs.umi_tag);
                    gs.umi_tag = strdup(optarg); break;
            case 4: gs.min_count = atoi(optarg); break;
            case 5: gs.min_maf = atof(optarg); break;
            case 6: gs.doublet_gl = 1; break;
            case 8: gs.min_len = atoi(optarg); break;
            case 9: gs.min_mapq = atoi(optarg); break;
            //case 10: gs.max_flag = atoi(optarg); break;
            case 11: gs.is_genotype = 1; break;
            case 12: gs.is_out_zip = 1; break;
            case 13: print_skip_snp = 1; break;
            case 14: 
                    if ((gs.rflag_require = bam_str2flag(optarg)) < 0) {
                        fprintf(stderr, "[E::%s] could not parse --inclFLAG '%s'\n", __func__, optarg);
                        goto clean;
                    } else { break; }
            case 15:
                    if ((gs.rflag_filter = bam_str2flag(optarg)) < 0) {
                        fprintf(stderr, "[E::%s] could not parse --exclFLAG '%s'\n", __func__, optarg);
                        goto clean;
                    } else { break; }
            case 16: gs.no_orphan = 0; break;
            case 17: gs.max_depth = atoi(optarg); break;
            case 18: gs.max_pileup = atoi(optarg); break;
            default: fprintf(stderr, "Invalid option: '%c'\n", c); goto clean;	
        }
    }

    fprintf(stderr, "[I::%s] start time: %s\n", __func__, time_str);
    print_info = 1;

    /* check global settings */
  #if DEBUG
    fprintf(stderr, "[D::%s] global settings before checking:\n", __func__);
    gll_setting_print(stderr, &gs, "\t");
  #endif
    if ((ret = check_args(&gs)) < 0) { 
        fprintf(stderr, "[E::%s] error global settings\n", __func__);
        if (ret == -1) { print_usage(stderr); }
        goto clean;
    }
    fprintf(stderr, "[I::%s] global settings after checking:\n", __func__);
    gll_setting_print(stderr, &gs, "\t");

    /* prepare output files. */
    if (NULL == (gs.out_mtx_ad = jf_init()) || NULL == (gs.out_mtx_dp = jf_init()) || \
        NULL == (gs.out_mtx_oth = jf_init()) || NULL == (gs.out_samples = jf_init()) || \
        NULL == (gs.out_vcf_base = jf_init()) || (gs.is_genotype && NULL == (gs.out_vcf_cells = jf_init()))) {
        fprintf(stderr, "[E::%s] fail to create jfile_t.\n", __func__);
        goto clean;
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
        goto clean;
    }
    if (output_headers(gs.out_mtx_dp, "wb", ks_str(s), ks_len(s)) < 0) {   // output header to mtx_DP
        fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_mtx_dp->fn);
        goto clean;
    }
    if (output_headers(gs.out_mtx_oth, "wb", ks_str(s), ks_len(s)) < 0) {  // output header to mtx_OTH
        fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_mtx_oth->fn);
        goto clean;
    } ks_clear(s);

    if (use_barcodes(&gs)) {                     // output samples.
        for (k = 0; k < gs.nbarcode; k++) {
            kputs(gs.barcodes[k], s);
            kputc('\n', s);
        }
    } else if (use_sid(&gs)) {
        for (k = 0; k < gs.nsid; k++) {
            kputs(gs.sample_ids[k], s);
            kputc('\n', s);
        }
    } // else: should not come here!

    if (output_headers(gs.out_samples, "wb", ks_str(s), ks_len(s)) < 0) {
        fprintf(stderr, "[E::%s] fail to write samples to '%s'\n", __func__, gs.out_samples->fn);
        goto clean;
    } ks_clear(s);
    kputs(CSP_VCF_BASE_HEADER, s);             // output header to vcf base.
    kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", s);
    if (output_headers(gs.out_vcf_base, "wb", ks_str(s), ks_len(s)) < 0) {
        fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_vcf_base->fn);
        goto clean;
    } ks_clear(s);
    if (gs.is_genotype) {
        kputs(CSP_VCF_CELLS_HEADER, s);           // output header to vcf cells.
        kputs(CSP_VCF_CELLS_CONTIG, s);
        kputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", s);
        if (use_barcodes(&gs) && gs.barcodes) {
            for (k = 0; k < gs.nbarcode; k++) {
                kputc_('\t', s);
                kputs(gs.barcodes[k], s);
            }
        } else if (use_sid(&gs) && gs.sample_ids) {
            for (k = 0; k < gs.nsid; k++) {
                kputc_('\t', s);
                kputs(gs.sample_ids[k], s);
            }
        } else {
            fprintf(stderr, "[E::%s] neither barcodes or sample IDs exist.\n", __func__);
            goto clean;
        }
        kputc('\n', s);
        if (output_headers(gs.out_vcf_cells, "wb", ks_str(s), ks_len(s)) < 0) {
            fprintf(stderr, "[E::%s] fail to write header to '%s'\n", __func__, gs.out_vcf_cells->fn);
            goto clean;
        }
    }

    /* set file modes. */
    gs.out_mtx_ad->fm = gs.out_mtx_dp->fm = gs.out_mtx_oth->fm = "ab";
    gs.out_vcf_base->fm = "ab";
    if (gs.is_genotype) { gs.out_vcf_cells->fm = "ab"; }

    /* run based on the mode of input. 
        Mode 1a: fetch a list of SNPs for a single BAM/SAM file with barcodes.
        Mode 1a(-T): pileup a list of SNPs for a single BAM/SAM file with barcodes.
        Mode 1b: fetch a list of SNPs for one or multiple BAM/SAM files with sample IDs.
        Mode 1b(-T): pileup a list of SNPs for one or multiple BAM/SAM files with sample IDs.
        Mode 2a: pileup whole chromosome(s) for a single BAM/SAM file with barcodes.
        Mode 2b: pileup whole chromosome(s) for one or multiple BAM/SAM files with sample IDs.
    */
    if (gs.snp_list_file) {
        fprintf(stderr, "[I::%s] loading the VCF file for given SNPs ...\n", __func__);
        if (gs.is_target) {
            char **tmp; 
            int ntmp;
            size_t i, n;
            snp_t *snp = NULL;
            biallele_t **ale = NULL;
            if (NULL == (ale = (biallele_t**) calloc(1, sizeof(biallele_t*)))) {
                fprintf(stderr, "[E::%s] out of space for biallele_t**!\n", __func__);
                goto clean;
            }
            if (get_snplist(gs.snp_list_file, &gs.pl, &ret, print_skip_snp) <= 0 || ret < 0) {
                fprintf(stderr, "[E::%s] get SNP list from '%s' failed.\n", __func__, gs.snp_list_file);
                goto clean;
            } else { 
                n = kv_size(gs.pl);
                fprintf(stderr, "[I::%s] pileuping %ld candidate variants ...\n", __func__, n); 
            }
            if (NULL == (gs.targets = regidx_init(NULL, NULL, regidx_payload_free, sizeof(biallele_t*), NULL))) {
                fprintf(stderr, "[E::%s] failed to init regidx for '%s'.\n", __func__, gs.snp_list_file);
                goto clean;
            } 
            for (i = 0; i < n; i++) {
                snp = kv_A(gs.pl, i);
                if (NULL == (*ale = biallele_init())) { 
                    fprintf(stderr, "[E::%s] failed to create biallele_t.\n", __func__);
                    goto clean;
                } else {
                    (*ale)->ref = snp->ref; (*ale)->alt = snp->alt;
                }
                if (regidx_push(gs.targets, snp->chr, snp->chr + strlen(snp->chr) - 1, snp->pos, snp->pos, ale) < 0) {
                    fprintf(stderr, "[E::%s] failed to push regidx.\n", __func__);
                    free(ale); ale = NULL;    // FIXME!! *ale should be safely freed.
                    goto clean;
                } 
            } free(ale); ale = NULL;
            snplist_destroy(gs.pl);
            tmp = regidx_seq_names(gs.targets, &ntmp);
            if (gs.chroms) {
                str_arr_destroy(gs.chroms, gs.nchrom);
                gs.nchrom = 0;
            }
            if (NULL == (gs.chroms = (char**) calloc(ntmp, sizeof(char*)))) {
                fprintf(stderr, "[E::%s] failed to allocate space for gs.chroms\n", __func__);
                goto clean;
            }
            for (gs.nchrom = 0; gs.nchrom < ntmp; gs.nchrom++) {
                gs.chroms[gs.nchrom] = strdup(tmp[gs.nchrom]);
            }
            if (gs.barcodes) { 
                fprintf(stderr, "[I::%s] mode 1a(-T): pileup %d whole chromosomes in %d single cells.\n", __func__, gs.nchrom, gs.nbarcode);
                if (run_mode1a_T(&gs) < 0) {
                    fprintf(stderr, "[E::%s] running mode 1a(-T) failed.\n", __func__);
                    goto clean;
                }
            } else { 
                fprintf(stderr, "[I::%s] mode 1b(-T): pileup %d whole chromosomes in %d sample(s).\n", __func__, gs.nchrom, gs.nin);
                if (run_mode1b_T(&gs) < 0) {
                    fprintf(stderr, "[E::%s] running mode 1b(-T) failed.\n", __func__);
                    goto clean;
                }
            }
        } else {
            if (get_snplist(gs.snp_list_file, &gs.pl, &ret, print_skip_snp) <= 0 || ret < 0) {
                fprintf(stderr, "[E::%s] get SNP list from '%s' failed.\n", __func__, gs.snp_list_file);
                goto clean;
            } else {
                fprintf(stderr, "[I::%s] fetching %ld candidate variants ...\n", __func__, kv_size(gs.pl));
            }
            if (gs.barcodes) { 
                fprintf(stderr, "[I::%s] mode 1a: fetch given SNPs in %d single cells.\n", __func__, gs.nbarcode); 
                if (run_mode1a(&gs) < 0) {
                    fprintf(stderr, "[E::%s] running mode 1a failed.\n", __func__);
                    goto clean;
                } 
            } else { 
                fprintf(stderr, "[I::%s] mode 1b: fetch given SNPs in %d bulk samples.\n", __func__, gs.nsid);
                if (run_mode1b(&gs) < 0) {
                    fprintf(stderr, "[E::%s] running mode 1b failed.\n", __func__);
                    goto clean;
                } 
            }
        }
    } else if (gs.chroms) { 
        if (gs.barcodes) { 
            fprintf(stderr, "[I::%s] mode 2a: pileup %d whole chromosomes in %d single cells.\n", __func__, gs.nchrom, gs.nbarcode);
            if (run_mode2a(&gs) < 0) {
                fprintf(stderr, "[E::%s] running mode 2a failed.\n", __func__);
                goto clean;
            }
        } else {
            fprintf(stderr, "[I::%s] mode 2b: pileup %d whole chromosomes in %d sample(s).\n", __func__, gs.nchrom, gs.nin);
            if (run_mode2b(&gs) < 0) {
                fprintf(stderr, "[E::%s] running mode 2b failed.\n", __func__);
                goto clean;
            }
        }
    } else {
        fprintf(stderr, "[E::%s] no proper mode to run, check input options.\n", __func__);
        print_usage(stderr);
        goto clean;
    }

    is_ok = 1;

  clean:
    /* clean */
    if (s) { ks_free(s); }
    gll_setting_free(&gs);
    if (print_info) {
        if (is_ok)
            fprintf(stderr, "[I::%s] All Done!\n", __func__);
        else
            fprintf(stderr, "[E::%s] Quiting...\n", __func__);

        /* command line */
        fprintf(stderr, "[I::%s] Version: %s (htslib %s)\n", __func__, CSP_VERSION, hts_version());
        fprintf(stderr, "[I::%s] CMD: %s", __func__, argv[0]);
        for (i = 1; i < argc; i++)
            fprintf(stderr, " %s", argv[i]);
        fputc('\n', stderr);

        /* calc time spent */
        time(&end_time);
        time_info = localtime(&end_time);
        strftime(time_str, 30, "%Y-%m-%d %H:%M:%S", time_info);
        fprintf(stderr, "[I::%s] end time: %s\n", __func__, time_str);
        fprintf(stderr, "[I::%s] time spent: %ld seconds.\n", __func__, end_time - start_time);
    }

    if (is_ok)
        return 0;
    else
        return 1;
}

