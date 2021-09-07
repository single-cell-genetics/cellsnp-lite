/* Pileup and MPileup operatoins API/routine
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_MPLP_H
#define CSP_MPLP_H

#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "kvec.h"
#include "jfile.h"
#include "jmempool.h"
#include "config.h"

/*
* Pileup and MPileup API
 */
/*@abstract          Store statistics result of one read matching the query pos.
@param b             Pointer of bam1_t structure.
@param qpos          The index of the query pos in the read seq.
@param base          The base corresponding to the query pos, equal to read_seq[qpos].
@param qual          The qual corresponding to the query pos, equal to qual_seq[qpos].
@param is_refskip    0/1. if the query pos is in res-skip region.
@param is_del        0/1. if the query pos is in deletion region (also set 1 when is_refskip to be compatible with sam.c).
@param umi           Pointer to UMI tag.
@param cb            Pointer to cell barcode.
@param laln          Length of the read part that aligned to reference.

@note  1. The umi and cb in the structure would be extracted from bam1_t, so do not need to free it as the bam1_t will do that!
          Refer to bam_get_aux(), bam_aux_get() and bam_aux2Z() in sam.h.
 */
typedef struct { 
    bam1_t *b;
    int32_t qpos;
    int8_t base, qual;
    uint8_t is_refskip:1, is_del:1;
    char *umi, *cb;
    uint32_t laln;
} csp_pileup_t;

/*@abstract   Initialize the csp_pileup_t structure.
@return       Pointer to csp_pileup_t structure if success, NULL otherwise.

@note         1. Memory for bam1_t is allocated in this function.
              2. The pointer returned successfully by csp_pileup_init() should be freed
                 by csp_pileup_destroy() when no longer used.
 */
inline csp_pileup_t* csp_pileup_init(void); 

inline void csp_pileup_destroy(csp_pileup_t *p); 

/* reset the csp_pileup_t structure without reallocating memory.
   return 0 if success, -1 otherwise. */
inline int csp_pileup_reset(csp_pileup_t *p);

/* only reset part of the csp_pileup_t as values of parameters in other parts will be immediately overwritten after
   calling this function. It's often called by pileup_read_with_fetch(). */
inline void csp_pileup_reset_(csp_pileup_t *p); 

inline void csp_pileup_print(FILE *fp, csp_pileup_t *p);

/*@abstract   Pool that stores char*. The real elements in the pool are char**.
@example      Refer to the example of JMEMPOOL in general_util.h.
@note         This structure is aimed to store dynamically allocated strings(char*).
 */
#define pool_str_free(s) free(*(s))
#define pool_str_reset(s) free(*(s))
JMEMPOOL_INIT(ps, char*, pool_str_free, pool_str_reset)
typedef jmempool_t(ps) pool_ps_t;
#define pool_ps_init() jmempool_init(ps)
#define pool_ps_destroy(p) jmempool_destroy(ps, p)
#define pool_ps_get(p) jmempool_get(ps, p)
#define pool_ps_reset(p) jmempool_reset(ps, p)

/*@abstract    This structure stores stat info of one read of one UMI group for certain query pos.
@param base    The base for the query pos in the read of the UMI gruop.
               A 4-bit integer returned by bam_seqi(), which is related to bam_nt16_table(now called seq_nt16_str).
@param qual    Qual value for the query pos in the read of the UMI gruop. The value is extracted by calling bam_get_qual() and 
               could be translated to qual char by plusing 33.
 */
typedef struct {
    int8_t base, qual;
} umi_unit_t;

inline umi_unit_t* umi_unit_init(void);
#define umi_unit_reset(p)
inline void umi_unit_destroy(umi_unit_t *p);

/*@abstract   Pool that stores umi_unit_t structures. The real elements in the pool are pointers to the umi_unit_t structures.
              The pool is aimed to save the overhead of reallocating memories for umi_unit_t structures.
@example      Refer to the example of JMEMPOOL in general_util.h.
 */
JMEMPOOL_INIT(uu, umi_unit_t, umi_unit_destroy, umi_unit_reset)
typedef jmempool_t(uu) pool_uu_t;
#define pool_uu_init() jmempool_init(uu)
#define pool_uu_destroy(p) jmempool_destroy(uu, p)
#define pool_uu_get(p) jmempool_get(uu, p)
#define pool_uu_reset(p) jmempool_reset(uu, p)

/* Struct list_uu_t APIs 
@abstract  The structure stores stat info of all reads of one UMI group for certain query pos.
@param v   Pointer to the list_uu_t structure [list_uu_t*].

@note      1. The elements in the list are pointers of umi_unit_t, resetting the list only sets
              the list's internal var n, which stands for the next available pos in the list, to be 0.
              After calling reset function, the values (i.e. pointers of umi_unit_t from pool_uu) 
              newly pushed into the list would overwrite the original ones, which will not cause memory error. 
           2. The parameter in list_uu_xxx functions is pointer of list_uu_t.

@example   Refer to the example in kvec.h, but the parameter in list_uu_xxx functions is pointer, 
           which is different from kv_xxx functions.
 */
typedef kvec_t(umi_unit_t*) list_uu_t;
inline list_uu_t* list_uu_init(void);
#define list_uu_resize(v, size) kv_resize(umi_unit_t*, *(v), size)
#define list_uu_push(v, x) kv_push(umi_unit_t*, *(v), x)
#define list_uu_A(v, i) kv_A(*(v), i)
#define list_uu_size(v) kv_size(*(v))
#define list_uu_max(v) kv_max(*(v))
#define list_uu_destroy(v) kv_destroy(*(v))
#define list_uu_reset(v) ((v)->n = 0)

/*@abstract   Pool that stores list_uu_t structures. The real elements in the pool are pointers to the list_uu_t structures.
              The pool is aimed to save the overhead of reallocating memories for list_uu_t structures.
@example      Refer to the example of JMEMPOOL in general_util.h.
 */
JMEMPOOL_INIT(ul, list_uu_t, list_uu_destroy, list_uu_reset)
typedef jmempool_t(ul) pool_ul_t;
#define pool_ul_init() jmempool_init(ul)
#define pool_ul_destroy(p) jmempool_destroy(ul, p)
#define pool_ul_get(p) jmempool_get(ul, p)
#define pool_ul_reset(p) jmempool_reset(ul, p)

/*@abstract    The HashMap maps UMI-group-name (char*) to list_uu_t (*).

@example (A simple example from khash.h)
KHASH_MAP_INIT_INT(32, char)
int main() {
    int ret, is_missing;
    khiter_t k;
    khash_t(32) *h = kh_init(32);
    k = kh_put(32, h, 5, &ret);
    kh_value(h, k) = 10;
    k = kh_get(32, h, 10);
    is_missing = (k == kh_end(h));
    k = kh_get(32, h, 5);
    kh_del(32, h, k);
    for (k = kh_begin(h); k != kh_end(h); ++k)
        if (kh_exist(h, k)) kh_value(h, k) = 1;
    kh_destroy(32, h);
    return 0;
}
 */

// TODO: use umi:array(list_uu_t*) pair to save memory
// specifically, the array is a[5] with each element being
// a list_uu_t* and the umi_unit_t could be changed to
// struct { int8_t qual; }.
KHASH_MAP_INIT_STR(ug, list_uu_t*)
typedef khash_t(ug) map_ug_t;
#define map_ug_iter khiter_t
#define map_ug_init() kh_init(ug)
#define map_ug_resize(h, s) kh_resize(ug, h, s)
#define map_ug_put(h, k, r)  kh_put(ug, h, k, r)
#define map_ug_get(h, k) kh_get(ug, h, k)
#define map_ug_del(h, k) kh_del(ug, h, k)
#define map_ug_exist(h, x) kh_exist(h, x)
#define map_ug_key(h, x) kh_key(h, x)
#define map_ug_val(h, x) kh_val(h, x)
#define map_ug_begin(h) kh_begin(h)
#define map_ug_end(h) kh_end(h)
#define map_ug_size(h) kh_size(h)
#define map_ug_reset(h) kh_clear(ug, h)
#define map_ug_destroy(h) {								\
    if (h) {										\
        map_ug_iter __k;								\
        for (__k = map_ug_begin(h); __k != map_ug_end(h); __k++) { 			\
            if (map_ug_exist(h, __k)) list_uu_destroy(map_ug_val(h, __k));		\
        }											\
        kh_destroy(ug, h);									\
    }												\
}

/* Struct list_qu_t APIs 
@abstract  The structure stores all qual value of one sample for certain query pos.
@param v   The list_qu_t structure [list_qu_t].

@example   Refer to the example in kvec.h.
 */
typedef kvec_t(int8_t) list_qu_t;
#define list_qu_init(v) kv_init(v)
#define list_qu_resize(v, size) kv_resize(int8_t, v, size)
#define list_qu_push(v, x) kv_push(int8_t, v, x)
#define list_qu_A(v, i) kv_A(v, i)
#define list_qu_size(v) kv_size(v)
#define list_qu_max(v) kv_max(v)
#define list_qu_destroy(v) kv_destroy(v)
#define list_qu_reset(v) ((v).n = 0)

/*@abstract    Internal function to convert the base call quality score to related values for different genotypes.
@param qual    Qual value for the query pos in the read of the UMI gruop. The value is extracted by calling bam_get_qual() and 
               could be translated to qual char by plusing 33.
@param cap_bq   
@param min_bq
@param rv      Pointer of vector of loglikelihood for AA, AA+AB (doublet), AB, B or E (see Demuxlet paper online methods)
               Size of RV must be no less than 4. Be careful that no size check inside the function
@return        0 if success, -1 otherwise.
@reference     1. http://emea.support.illumina.com/bulletins/2016/04/fastq-files-explained.html
               2. https://linkinghub.elsevier.com/retrieve/pii/S0002-9297(12)00478-8 
*/
inline int get_qual_vector(double qual, double cap_bq, double min_bq, double *rv);

/*@abstract     Internal function to translate qual matrix to vector of GL.
@param qm       Qual matrix: 5-by-4, columns are [1-Q, 3/4-2/3Q, 1/2-1/3Q, Q]. See csp_plp_t::qmat
@param bc       Base count: (5, ). See csp_plp_t::bc
@param ref_idx  Index of ref base in 'ACGTN'. 
@param alt_idx  Index of alt base in 'ACGTN'.
@param db       doublet_GL. 1/0.
@param gl       Array of GL: loglikelihood for 
                     GL1: L(rr|qual_matrix, base_count), 
                     GL2-GL5: L(ra|..), L(aa|..), L(rr+ra|..), L(ra+aa|..)
                Size must be no less than 5. Be careful that no size check inside the function
@param n        Pointer of num of elements in GL array.
@return         0 if success, -1 otherwise.
@note           TODO: In some special cases, ref=A and alt=AG for example, the ref_idx would be equal with alt_idx.
                Should be fixed in future.
 */
int qual_matrix_to_geno(double qm[][4], size_t *bc, int8_t ref_idx, int8_t alt_idx, int db, double *gl, int *n);

/*@abstract     Infer ref and alt based on the read count of each possible alleles (ACGTN).
@param bc       Pointer of array that stores the read count of each possible alleles, the read count is
                in the order of A/C/G/T/N.
@param ref_idx  Pointer of index of ref in "ACGTN".
@param alt_idx  Pointer of index of alt in "ACGTN".
@return         Void.

@note           It's usually called when the input pos has no ref or alt.
 */
inline void infer_allele(size_t *bc, int8_t *ref_idx, int8_t *alt_idx);

// infer alt index when the ref is known.
inline void infer_alt(size_t *bc, int8_t ref_idx, int8_t *alt_idx);

/*@abstract  The structure that store pileup info of one cell/sample for certain query pos.
@param bc    Total read count for each base in the order of 'ACGTN'.
@param tc    Total read count for all bases.
@param ad    Read count of alt.
@param dp    Read count of alt + ref.
@param oth   Read count of bases except alt and ref.
@param qu    All qual values of each base for one sample in the order of 'ACGTN'.
@param qmat  Matrix of qual with 'ACGTN' vs. [1-Q, 3/4-2/3Q, 1/2-1/3Q, Q].
@param gl    Array of GL: loglikelihood for 
               GL1: L(rr|qual_matrix, base_count), 
               GL2-GL5: L(ra|..), L(aa|..), L(rr+ra|..), L(ra+aa|..).
@param ngl   Num of valid elements in the array gl.
@param hug   Pointer of hash table that stores stat info of UMI groups.
 */
typedef struct {
    size_t bc[5];
    size_t tc, ad, dp, oth;
    list_qu_t qu[5];
    double qmat[5][4];
    double gl[5];
    int ngl;
    map_ug_t *hug;
} csp_plp_t;

/* note that the @p qu is also initialized after calling calloc(). */
inline csp_plp_t* csp_plp_init(void); 
inline void csp_plp_destroy(csp_plp_t *p); 
inline void csp_plp_reset(csp_plp_t *p);

/*@abstract    Print the content to csp_plp_t to stream.
@param fp      Pointer of stream.
@param p       Pointer of csp_plpt_t to be printed.
@param prefix  Pointer of prefix string. Set to "" if no prefix.
@return        Void.
 */
void csp_plp_print(FILE *fp, csp_plp_t *p, char *prefix);

/*@abstract     Format the content of csp_plp_t (the pileup stat info) of one cell/sample for certain query pos to 
                string in the output vcf file. 
@param p        Pointer of csp_plp_t structure corresponding to the pos.
@param s        Pointer of kstring_t which would store the formatted string.
@return         0 if success, -1 otherwise.
 */
int csp_plp_str_vcf(csp_plp_t *p, kstring_t *s);

int csp_plp_to_vcf(csp_plp_t *p, jfile_t *s);

/*@abstract    The HashMap maps sample-group-name (char*) to csp_plp_t (*).
@example       Refer to a simple example in khash.h.
 */
KHASH_MAP_INIT_STR(sg, csp_plp_t*)
typedef khash_t(sg) map_sg_t;
#define map_sg_iter khiter_t
#define map_sg_init() kh_init(sg)
#define map_sg_clear(h) kh_clear(sg, h)
#define map_sg_resize(h, s) kh_resize(sg, h, s)
#define map_sg_put(h, k, r)  kh_put(sg, h, k, r)
#define map_sg_get(h, k) kh_get(sg, h, k)
#define map_sg_del(h, k) kh_del(sg, h, k)
#define map_sg_exist(h, x) kh_exist(h, x)
#define map_sg_key(h, x) kh_key(h, x)
#define map_sg_val(h, x) kh_val(h, x)
#define map_sg_begin(h) kh_begin(h)
#define map_sg_end(h) kh_end(h)
#define map_sg_size(h) kh_size(h)
#define map_sg_destroy(h) {								\
    if (h) {											\
        map_sg_iter __k;									\
        for (__k = map_sg_begin(h); __k != map_sg_end(h); __k++) { 			\
            if (map_sg_exist(h, __k)) csp_plp_destroy(map_sg_val(h, __k)); 	\
        }										\
        kh_destroy(sg, h);									\
    }												\
}
#define map_sg_reset_val(h) {								\
    if (h) {											\
        map_sg_iter __k;									\
        for (__k = map_sg_begin(h); __k != map_sg_end(h); __k++) {			\
            if (map_sg_exist(h, __k)) csp_plp_reset(map_sg_val(h, __k)); 		\
        }											\
    }												\
}

//TODO: use kstring_t instead of char* in mplp_t::su as its base element type.
//Previously, we use a pool (actually as a array) to store all UMI strings (char*)
//pushed to plp_t::hug for easy management (mainly for easier free of the 
//UMI strings after each round of pileup). The overhead of frequent create (with @func 
//strdup and free (with @func free()) may be high. Using kstring_t is expected to
//largely reduce this overhead

/*@abstract  The structure stores the stat info of all sample groups for certain query pos.
@param ref_idx  Index of ref in "ACGTN". Negative number means not valid value.
@param alt_idx  Index of alt in "ACGTN". Negative number means not valid value.
@param inf_rid  Infered index of ref in "ACGTN". Negative number means not valid value.
@param inf_aid  Infered index of alt in "ACGTN". Negative number means not valid value.
@param bc    Read count of each base summarizing all sample groups for the pos, in the order of 'ACGTN'.
@param tc    Total read count of all bases for the pos.
@param ad    Read count of alt.
@param dp    Read count of alt + ref.
@param oth   Read count of bases except alt and ref.
@param nr_*  Num of records/lines outputed to mtx file for certain SNP/pos.
@param hsg   HashMap that stores the stat info of all sample groups for the pos.
@param hsg_iter Pointer of array of map_sg_iter. The iter in the array is in the same order of sg names.
@param nsg   Size of map_sg_iter array hsg_iter.
@param pu    Pool of umi_unit_t structures.
@param pl    Pool of list_uu_t structures.
@param su    Pool of UMI strings.
@param qvec  A container for the qual vector returned by get_qual_vector().
 */
typedef struct {
    int8_t ref_idx, alt_idx, inf_rid, inf_aid;
    size_t bc[5];
    size_t tc, ad, dp, oth;
    size_t nr_ad, nr_dp, nr_oth;
    map_sg_t *hsg;
    map_sg_iter *hsg_iter;
    int nsg;
    pool_uu_t *pu;
    pool_ul_t *pl;
    pool_ps_t *su;
    double qvec[4];
} csp_mplp_t;

/*@abstract  Initialize the csp_mplp_t structure.
@return      Pointer to the csp_mplp_t structure if success, NULL otherwise.

@note        1. The kstring_t s is also initialized inside this function.   
             2. The valid pointer returned by this function should be freed by csp_mplp_destroy() function
                   when no longer used.
 */
inline csp_mplp_t* csp_mplp_init(void); 
inline void csp_mplp_destroy(csp_mplp_t *p); 
inline void csp_mplp_reset(csp_mplp_t *p);

/*@abstract    Print the content to csp_mplp_t to stream.
@param fp      Pointer of stream.
@param p       Pointer of csp_mplpt_t to be printed.
@param prefix  Pointer of prefix string. Set to "" if no prefix.
@return        Void.
 */
void csp_mplp_print(FILE *fp, csp_mplp_t *p, char *prefix);

inline void csp_mplp_print_(FILE *fp, csp_mplp_t *p, char *prefix);

/*@abstract  Set sample group names for the csp_mplp_t structure.
@param p     Pointer to the csp_mplp_t structure.
@parma s     Pointer to the array of names of sample groups.
@param n     Num of sample groups.
@return      0, no error; -1 otherwise.

@note        1. This function should be called just one time right after csp_mplp_t structure was created
                becuase the sgname wouldn't change once set.
             2. The HashMap (for sgnames) in csp_mplp_t should be empty or NULL.
             3. The keys of HashMap are exactly pointers to sg names coming directly from @p s.
 */
int csp_mplp_set_sg(csp_mplp_t *p, char **s, const int n);

/*@abstract  Format the content of csp_mplp_t (the pileup stat info) of certain query pos to string in the output vcf file.
@param mplp  Pointer of the csp_mplp_t structure corresponding to the pos.
@param s     Pointer of kstring_t which stores the formatted string.
@return      0 if success, -1 otherwise.
 */
inline int csp_mplp_str_vcf(csp_mplp_t *mplp, kstring_t *s);

inline int csp_mplp_to_vcf(csp_mplp_t *mplp, jfile_t *s);

/*@abstract    Format the content of csp_mplp_t of certain query pos to string in the output sparse matrices file.
@param mplp    Pointer of the csp_mplp_t structure corresponding to the pos.
@param ks_ad   Pointer of kstring_t which is to store formatted AD string.
@param ks_dp   Pointer of kstring_t which is to store formatted DP string.
@param ks_oth  Pointer of kstring_t which is to store formatted OTH string.
@param idx     Index of the SNP/mplp (1-based).
@return        0 if success, -1 otherwise.
 */
inline int csp_mplp_str_mtx(csp_mplp_t *mplp, kstring_t *ks_ad, kstring_t *ks_dp, kstring_t *ks_oth, size_t idx);

/*@abstract    Format the content of csp_mplp_t of certain query pos to string in the tmp output sparse matrices file.
@param mplp    Pointer of the csp_mplp_t structure corresponding to the pos.
@param ks_ad   Pointer of kstring_t which is to store formatted AD string.
@param ks_dp   Pointer of kstring_t which is to store formatted DP string.
@param ks_oth  Pointer of kstring_t which is to store formatted OTH string.
@return        0 if success, -1 otherwise.

@note          This function is used for tmp files.
 */
inline int csp_mplp_str_mtx_tmp(csp_mplp_t *mplp, kstring_t *ks_ad, kstring_t *ks_dp, kstring_t *ks_oth);

int csp_mplp_to_mtx(csp_mplp_t *mplp, jfile_t *fs_ad, jfile_t *fs_dp, jfile_t *fs_oth, size_t idx); 

#if DEVELOP
/* 
* Tags for sparse matrices
* It's useful when the input tags are optional.
*/
typedef size_t mtx_value_t;
typedef int mtx_iter_t;
typedef mtx_value_t (*mtx_value_func_t)(csp_mplp_t*, csp_plp_t*);

inline mtx_value_t mtx_value_AD(csp_mplp_t *mplp, csp_plp_t *plp) {
    return plp->ad;
}

inline mtx_value_t mtx_value_DP(csp_mplp_t *mplp, csp_plp_t *plp) {
    return plp->dp;
}

inline mtx_value_t mtx_value_OTH(csp_mplp_t *mplp, csp_plp_t *plp) {
    return plp->oth;
}

static mtx_iter_t mtx_ntags = 3;
static char *mtx_tags[] = {"AD", "DP", "OTH"};
static mtx_value_func_t mtx_value_funcs[] = { &mtx_value_AD, &mtx_value_DP, &mtx_value_OTH };

inline char* csp_get_mtx_fn(char *tag) {
    kstring_t ks = KS_INITIALIZE;
    ksprintf(&ks, "cellSNP.tag.%s.mtx", tag);
    char *p = strdup(ks_str(&ks));
    ks_free(&ks);
    return p;
}

inline mtx_iter_t csp_get_mtx_idx(char *tag) {
    mtx_iter_t i;
    for (i = 0; i < mtx_ntags; i++) {
        if (0 == strcmp(tag, mtx_tags[i])) { return i; }
    } 
    return -1;
}

inline mtx_value_func_t csp_get_mtx_value_func(mtx_iter_t i) { return mtx_value_funcs[i]; }

typedef struct {
    char *out_fn;
    mtx_value_func_t vfunc;
} mtx_tag_fs;

inline mtx_tag_fs* mtx_tag_fs_init(void) {
    return (mtx_tag_fs*) calloc(1, sizeof(mtx_tag_fs));
}

inline void mtx_tag_fs_destroy(mtx_tag_fs *p) {
    if (p) { free(p->out_fn); free(p); }
}
#endif

#endif
