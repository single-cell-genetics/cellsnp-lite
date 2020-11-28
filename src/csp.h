/* Utils
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_CSP_H
#define CSP_CSP_H

#include <stdio.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "config.h"
#include "mplp.h"
#include "jfile.h"
#include "snp.h"
#include "thpool.h"


/* 
 * Global settings
 */

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
    csp_snplist_t pl;      // List of the input SNPs. TODO: local variable.
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
    int nthread;           // Num of threads.
    threadpool tp;         // Pointer to thread pool.
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
    int plp_max_depth;      // max depth for one site of one file, 0 means highest possible value.
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
           c) each key (sample group name) in csp_map_sg_t already has a valid, not NULL, value (csp_plp_t*);
              This usually can be done by calling csp_mplp_prepare().
        2. This function is expected to be used by Mode1 & Mode2 & Mode3.

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
inline csp_bam_fs* csp_bam_fs_init(void);
inline void csp_bam_fs_destroy(csp_bam_fs* p);

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
inline thread_data* thdata_init(void);
inline void thdata_destroy(thread_data *p);
inline void thdata_print(FILE *fp, thread_data *p);

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
inline jfile_t* create_tmp_fs(jfile_t *fs, int idx, int is_zip, kstring_t *s);

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
inline int destroy_tmp_files(jfile_t **fs, const int n);

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
