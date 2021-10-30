/* Utils
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_CSP_H
#define CSP_CSP_H

#include <stdio.h>
#include <limits.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/regidx.h"
#include "config.h"
#include "mplp.h"
#include "jfile.h"
#include "snp.h"
#include "thpool.h"

/* Define default values of global parameters. */
#define CSP_CHROM_ALL  {"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"}
#define CSP_NCHROM     22
#define CSP_CELL_TAG   "CB"
#define CSP_UMI_TAG    "Auto"
#define CSP_NTHREAD    1
#define CSP_MIN_COUNT  20
#define CSP_MIN_MAF    0.0
#define CSP_MIN_LEN    30
#define CSP_MIN_MAPQ   20
/* the following three macros are deprecated, use CSP_EXCL_FMASK* and CSP_INCL_FMASK* instead.
#define CSP_MAX_SAM_FLAG         4096                
#define CSP_MAX_FLAG_WITH_UMI    CSP_MAX_SAM_FLAG
#define CSP_MAX_FLAG_WITHOUT_UMI 255 
*/

#define CSP_OUT_VCF_CELLS   "cellSNP.cells.vcf"
#define CSP_OUT_VCF_BASE    "cellSNP.base.vcf"
#define CSP_OUT_SAMPLES     "cellSNP.samples.tsv"
#define CSP_OUT_MTX_AD      "cellSNP.tag.AD.mtx"
#define CSP_OUT_MTX_DP      "cellSNP.tag.DP.mtx"
#define CSP_OUT_MTX_OTH     "cellSNP.tag.OTH.mtx"

/* default values of pileup */
// default excluding flag mask, reads with any flag mask bit set would be filtered.
#define CSP_EXCL_FMASK_UMI  (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL)
#define CSP_EXCL_FMASK_NOUMI  (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)
// default including flag mask, reads with all flag mask bit unset would be filtered.
#define CSP_INCL_FMASK  0
// default max depth for one site of one bam file, 0 means highest possible value.
// It will be used by csp_fetch and the bam_mplp_set_maxcnt() in csp_pileup.
#define CSP_MAX_DEPTH   0
// if discard orphan reads
#define CSP_NO_ORPHAN   1

// if the tmp files to be zipped: 0: no, 1: yes.
#define CSP_TMP_ZIP 1

// output settings
#define CSP_VCF_CELLS_HEADER "##fileformat=VCFv4.2\n" 			\
    "##source=cellSNP_v" CSP_VERSION "\n"				\
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"		\
    "##FILTER=<ID=.,Description=\"Filter info not available\">\n"	\
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"total counts for ALT and REF\">\n" 	\
    "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"total counts for ALT\">\n"		\
    "##INFO=<ID=OTH,Number=1,Type=Integer,Description=\"total counts for other bases from REF and ALT\">\n"	\
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"						\
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"	\
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"total counts for ALT and REF\">\n"			\
    "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"total counts for ALT\">\n"				\
    "##FORMAT=<ID=OTH,Number=1,Type=Integer,Description=\"total counts for other bases from REF and ALT\">\n"		\
    "##FORMAT=<ID=ALL,Number=5,Type=Integer,Description=\"total counts for all bases in order of A,C,G,T,N\">\n"

#define CSP_VCF_CELLS_CONTIG "##contig=<ID=1>\n##contig=<ID=2>\n##contig=<ID=3>\n##contig=<ID=4>\n##contig=<ID=5>\n"        \
    "##contig=<ID=6>\n##contig=<ID=7>\n##contig=<ID=8>\n##contig=<ID=9>\n##contig=<ID=10>\n"					\
    "##contig=<ID=11>\n##contig=<ID=12>\n##contig=<ID=13>\n##contig=<ID=14>\n##contig=<ID=15>\n"				\
    "##contig=<ID=16>\n##contig=<ID=17>\n##contig=<ID=18>\n##contig=<ID=19>\n##contig=<ID=20>\n"				\
    "##contig=<ID=21>\n##contig=<ID=22>\n##contig=<ID=X>\n##contig=<ID=Y>\n"

#define CSP_MTX_HEADER "%%MatrixMarket matrix coordinate integer general\n"           \
    "%\n"

#define CSP_VCF_BASE_HEADER "##fileformat=VCFv4.2\n"

#if CSP_FIT_MULTI_SMP
    #define TP_EUNDEF       1        // error that undefined
    #define TP_EMFILE      (1 << 1)  // error that too many open files
#endif
#define TP_MAX_OPEN         1024     // default max number of open files

/* 
 * Global settings
 */
typedef struct _gll_settings global_settings;

/*Structure that stores global settings/options/parameters.
Note:
1. In current version, one and only one of barcode(s) and sample-ID(s) would exist and work, the other
   would be freed. Refer to check_global_args() for details.
*/
struct _gll_settings {
    char *in_fn_file;      // Name of the file containing a list of input bam/sam/cram files, one input file per line.
    int nin;               // Num of input bam/sam/cram files.
    char **in_fns;         // Pointer to the array of names of input bam/sam/cram files.
    char *out_dir;         // Pointer to the path of dir containing the output files.
    jfile_t *out_vcf_cells, *out_vcf_base, *out_samples;
    jfile_t *out_mtx_ad, *out_mtx_dp, *out_mtx_oth;
    int is_out_zip;        // If output files need to be zipped.
    int is_genotype;       // If need to do genotyping in addition to counting.
    char *snp_list_file;   // Name of file containing a list of SNPs, usually a vcf file.
    snplist_t pl;      // List of the input SNPs. TODO: local variable.
    int is_target;         // If the provided snp list should be used as target (like -T in samtools/bcftools mpileup). 1, yes; 0, no
    regidx_t *targets;     // Target regions.
    char *barcode_file;    // Name of the file containing a list of barcodes, one barcode per line.
    char **barcodes;       // Pointer to the array of barcodes.
    int nbarcode;          // Num of the barcodes.
    char *sid_list_file;   // Name of the file containing a list of sample IDs, one sample-ID per line.
    char **sample_ids;     // Pointer to the array of sample IDs.
    int nsid;              // Num of sample IDs.
    char **chroms;      // Pointer to the array of the chromosomes to use.
    int nchrom;            // Num of chromosomes.
    char *cell_tag;        // Tag for cell barcodes, NULL means no cell tags.
    char *umi_tag;         // Tag for UMI: UR, NULL. NULL means no UMI but read counts.
    int nthread;           // Num of threads to be used.
    threadpool tp;         // Pointer to thread pool.
    int mthread;           // Num of threads that user specified.
    int tp_errno;          // Error number, each bit could be used.
    int tp_ntry;           // Num of try
    int tp_max_open;       // Max num of open files for one process
    int min_count;     // Minimum aggragated count.
    double min_maf;    // Minimum minor allele frequency.
    int double_gl;     // 0 or 1. 1: keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5. 0: not keep.
    int min_len;       // Minimum mapped length for read filtering.
    int min_mapq;      // Minimum MAPQ for read filtering.
    /* max_flag is deprecated. use rflag_filter and rflag_require instead.
 *     int max_flag;      // Maximum FLAG for read filtering.
 *         */
    int rflag_filter;   // excluding flag mask, reads with any flag mask bit set would be filtered.
    int rflag_require;  // including flag mask, reads with all flag mask bit unset would be filtered.
    int max_depth;      // max depth for one site of one file, 0 means highest possible value.
    int no_orphan;     // 0 or 1. 1: donot use orphan reads; 0: use orphan reads.
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

/*@abstract  Whether to use target regions (-T).
@param gs    Pointer of global settings structure [global_settings*]
@return      1, yes; 0, no.
*/
#define use_target(gs) ((gs)->is_target)

void gll_setting_free(global_settings *gs); 
void gll_setting_print(FILE *fp, global_settings *gs, char *prefix);

/*
 * Mpileup processing
 */

/*@abstract  Set values for internal variables of csp_mplp_t to prepare for pileup. 
@param mplp  Pointer of csp_mplp_t structure.
@param gs    Pointer of global_settings structure.
@return      0 if success, -1 otherwise.
*/
int csp_mplp_prepare(csp_mplp_t *mplp, global_settings *gs);

/*@abstract    Push content of one csp_pileup_t structure into the csp_mplp_t structure.
@param pileup  Pointer of csp_pileup_t structure to be pushed.
@param mplp    Pointer of csp_mplp_t structure pushing into.
@param sid     Index of Sample ID in the input Sample IDs.
@param gs      Pointer of global_settings structure.
@return        0 if success;
               Negative numbers for error:
                 -1, neither barcodes or Sample IDs are used.
                 -2, khash_put error.
               Positive numbers for warning:
                 1, cell-barcode is not in input barcode-list;

@note   1. To speed up, the caller should guarantee that:
           a) the parameters are valid, i.e. mplp and gs must not be NULL. In fact, this function is supposed to be 
              called after csp_mplp_t is created and set names of sample-groups, so mplp, mplp->hsg could not be NULL.
           b) the csp_pileup_t must have passed the read filtering, refer to pileup_read_with_fetch() for details.
           c) each key (sample group name) in map_sg_t already has a valid, not NULL, value (csp_plp_t*);
              This usually can be done by calling csp_mplp_prepare().
        2. This function is expected to be used by Mode1 & Mode2.

@discuss  In current version, only the result (base and qual) of the first read in one UMI group will be used for mplp statistics.
          TODO: store results of all reads in one UMI group (maybe could do consistency correction in each UMI group) and then 
          do mplp statistics.
 */
int csp_mplp_push(csp_pileup_t *pileup, csp_mplp_t *mplp, int sid, global_settings *gs);

/*@abstract    Do statistics and filtering after all pileup results have been pushed.
@param mplp    Pointer of csp_mplp_t structure.
@param gs      Pointer of global_settings structure.
@return        0 if success; -1 if error; 1 if not passing filters.

@discuss  In current version, only the result (base and qual) of the first read in one UMI group will be used for mplp statistics.
          TODO: store results of all reads in one UMI group (maybe could do consistency correction in each UMI group) and then 
          do mplp statistics.
 */
int csp_mplp_stat(csp_mplp_t *mplp, global_settings *gs);

/*
* BAM/SAM/CRAM File
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
csp_bam_fs* csp_bam_fs_init(void);
void csp_bam_fs_destroy(csp_bam_fs* p);

/* 
 * Thread operatoins API/routine
 */

/* 
* Thread API
*/
/*@abstract    The data structure used as thread-func parameter.
@param gs      Pointer to the global_settings structure.
@param bfs     Array of csp_bam_fs.
@param nfs     Size (Number of elements) of @p bfs.
@param iter    Array of hts_itr_t**. 
@param niter   Size of @p iter.
@param nitr    Size of one element of @p iter.
@param n       Pos of next element in the snp-list/chrom-list to be used by certain thread.
@param m       Total size of elements to be used by certain thread, must not be changed.
@param i       Id of the thread data.
@param ret     Running state of the thread.
@param ns      Num of SNPs that passed all filters.
@param nr_*    Num of records for each output matrix file. 
@param out_*   Pointers of output files.
 */
typedef struct {
    global_settings *gs;
    csp_bam_fs **bfs;
    int nfs;
    hts_itr_t ***iter;
    int niter, nitr;
    size_t m, n;   // for snp-list or chrom-list.
    int i;
    int ret;
    size_t ns, nr_ad, nr_dp, nr_oth;
    jfile_t *out_mtx_ad, *out_mtx_dp, *out_mtx_oth, *out_vcf_base, *out_vcf_cells;
} thread_data;

/*@abstract  Create the thread_data structure.
@return      Pointer to the structure if success, NULL otherwise.
@note        The pointer returned successfully by thdata_init() should be freed
             by thdata_destroy() when no longer used.
 */
thread_data* thdata_init(void);
void thdata_destroy(thread_data *p);
void thdata_print(FILE *fp, thread_data *p);

/*
 * File Routine
 */

/*@abstract    Create jfile_t structure for tmp file.
@param fs      The file struct that the tmp file is based on.
@param idx     A number as suffix.
@param is_zip  If the tmp files should be zipped.
@param s       Pointer of kstring_t.
@return        Pointer to jfile_t for tmp file if success, NULL otherwise.
 */
jfile_t* create_tmp_fs(jfile_t *fs, int idx, int is_zip, kstring_t *s);

/*@abstract    Create array of tmp filen structures based on the given file structure.
@param fs      The file struct that the tmp file structs are based on.
@param n       Number of tmp file structs to be created.
@param is_zip  If the tmp files should be zipped.
@return        Pointer to the array of tmp file structs if success, NULL otherwise. 
 */
jfile_t** create_tmp_files(jfile_t *fs, int n, int is_zip);

/*@abstract  Remove tmp files and free memory.
@param fs    Pointer of array of jfile_t structures to be removed and freed.
@param n     Size of array.
@return      Num of tmp files that are removed if no error, -1 otherwise.
 */
int destroy_tmp_files(jfile_t **fs, const int n);

/*@abstract   Merge several tmp sparse matrices files.
@param out    Pointer of file structure merged into.
@param in     Pointer of array of tmp mtx files to be merged.
@param n      Num of tmp mtx files.
@param ns     Pointer to num of SNPs in all input mtx files.
@param nr     Pointer to num of records in all input mtx files.
@param ret    Pointer to the running state. 0 if success, negative numbers for error:
                -1, unknown error;
                -2, I/O error.
@return       Num of tmp mtx files that are successfully merged.
*/
int merge_mtx(jfile_t *out, jfile_t **in, const int n, size_t *ns, size_t *nr, int *ret);

/*@abstract   Merge several tmp vcf files.
@param out    Pointer of file structure merged into.
@param in     Pointer of array of tmp vcf files to be merged.
@param n      Num of tmp vcf files.
@param ret    Pointer to the running state. 0 if success, negative numbers for error:
                -1, unknown error;
                -2, I/O error.
@return       Num of tmp vcf files that are successfully merged.
*/
int merge_vcf(jfile_t *out, jfile_t **in, const int n, int *ret);

/*@abstract  Rewrite mtx file to fill in the stat info.
@param fs    Pointer of jfile_t that to be rewriten.
@param ns    Num of SNPs.
@param nsmp  Num of samples.
@param nr    Num of records.
@return      0 if success, -1 if error, 1 if the original file has no records while nr != 0.

@note        1. When proc = 1, the origial outputed mtx file was not filled with stat info:
                (totol SNPs, total samples, total records),
                so use this function to fill and rewrite.
             2. @p fs is not open when this function is just called and will keep not open when this function ends.
 */
int rewrite_mtx(jfile_t *fs, size_t ns, int nsmp, size_t nr);

/*
 * Subcommands
 */
int csp_fetch(global_settings *gs);
int csp_pileup(global_settings *gs);

#endif
