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
#include "jmemory.h"
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
static inline csp_pileup_t* csp_pileup_init(void); 

static inline void csp_pileup_destroy(csp_pileup_t *p); 

/* reset the csp_pileup_t structure without reallocating memory.
   return 0 if success, -1 otherwise. */
static inline int csp_pileup_reset(csp_pileup_t *p);

/* only reset part of the csp_pileup_t as values of parameters in other parts will be immediately overwritten after
   calling this function. It's often called by pileup_read_with_fetch(). */
static inline void csp_pileup_reset_(csp_pileup_t *p); 

static inline void csp_pileup_print(FILE *fp, csp_pileup_t *p);

/*@abstract   Pool that stores char*. The real elements in the pool are char**.
@example      Refer to the example of SZ_POOL in general_util.h.
@note         This structure is aimed to store dynamically allocated strings(char*).
 */
#define csp_pool_str_free(s) free(*(s))
#define csp_pool_str_reset(s) free(*(s))
SZ_POOL_INIT(ps, char*, csp_pool_str_free, csp_pool_str_reset)
typedef sz_pool_t(ps) csp_pool_ps_t;
#define csp_pool_ps_init() sz_pool_init(ps)
#define csp_pool_ps_destroy(p) sz_pool_destroy(ps, p)
#define csp_pool_ps_get(p) sz_pool_get(ps, p)
#define csp_pool_ps_reset(p) sz_pool_reset(ps, p)

/*@abstract    This structure stores stat info of one read of one UMI group for certain query pos.
@param base    The base for the query pos in the read of the UMI gruop.
               A 4-bit integer returned by bam_seqi(), which is related to bam_nt16_table(now called seq_nt16_str).
@param qual    Qual value for the query pos in the read of the UMI gruop. The value is extracted by calling bam_get_qual() and 
               could be translated to qual char by plusing 33.
 */
typedef struct {
    int8_t base, qual;
} csp_umi_unit_t;

static inline csp_umi_unit_t* csp_umi_unit_init(void);
#define csp_umi_unit_reset(p)
static inline void csp_umi_unit_destroy(csp_umi_unit_t *p);

/*@abstract   Pool that stores csp_umi_unit_t structures. The real elements in the pool are pointers to the csp_umi_unit_t structures.
              The pool is aimed to save the overhead of reallocating memories for csp_umi_unit_t structures.
@example      Refer to the example of SZ_POOL in general_util.h.
 */
SZ_POOL_INIT(uu, csp_umi_unit_t, csp_umi_unit_destroy, csp_umi_unit_reset)
typedef sz_pool_t(uu) csp_pool_uu_t;
#define csp_pool_uu_init() sz_pool_init(uu)
#define csp_pool_uu_destroy(p) sz_pool_destroy(uu, p)
#define csp_pool_uu_get(p) sz_pool_get(uu, p)
#define csp_pool_uu_reset(p) sz_pool_reset(uu, p)

/* Struct csp_list_uu_t APIs 
@abstract  The structure stores stat info of all reads of one UMI group for certain query pos.
@param v   Pointer to the csp_list_uu_t structure [csp_list_uu_t*].

@note      1. The elements in the list are pointers of csp_umi_unit_t, resetting the list only sets
              the list's internal var n, which stands for the next available pos in the list, to be 0.
              After calling reset function, the values (i.e. pointers of csp_umi_unit_t from csp_pool_uu) 
              newly pushed into the list would overwrite the original ones, which will not cause memory error. 
           2. The parameter in csp_list_uu_xxx functions is pointer of csp_list_uu_t.

@example   Refer to the example in kvec.h, but the parameter in csp_list_uu_xxx functions is pointer, 
           which is different from kv_xxx functions.
 */
typedef kvec_t(csp_umi_unit_t*) csp_list_uu_t;
static inline csp_list_uu_t* csp_list_uu_init(void);
#define csp_list_uu_resize(v, size) kv_resize(csp_umi_unit_t*, *(v), size)
#define csp_list_uu_push(v, x) kv_push(csp_umi_unit_t*, *(v), x)
#define csp_list_uu_A(v, i) kv_A(*(v), i)
#define csp_list_uu_size(v) kv_size(*(v))
#define csp_list_uu_max(v) kv_max(*(v))
#define csp_list_uu_destroy(v) kv_destroy(*(v))
#define csp_list_uu_reset(v) ((v)->n = 0)

/*@abstract   Pool that stores csp_list_uu_t structures. The real elements in the pool are pointers to the csp_list_uu_t structures.
              The pool is aimed to save the overhead of reallocating memories for csp_list_uu_t structures.
@example      Refer to the example of SZ_POOL in general_util.h.
 */
SZ_POOL_INIT(ul, csp_list_uu_t, csp_list_uu_destroy, csp_list_uu_reset)
typedef sz_pool_t(ul) csp_pool_ul_t;
#define csp_pool_ul_init() sz_pool_init(ul)
#define csp_pool_ul_destroy(p) sz_pool_destroy(ul, p)
#define csp_pool_ul_get(p) sz_pool_get(ul, p)
#define csp_pool_ul_reset(p) sz_pool_reset(ul, p)

/*@abstract    The HashMap maps UMI-group-name (char*) to csp_list_uu_t (*).

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
KHASH_MAP_INIT_STR(ug, csp_list_uu_t*)
typedef khash_t(ug) csp_map_ug_t;
#define csp_map_ug_iter khiter_t
#define csp_map_ug_init() kh_init(ug)
#define csp_map_ug_resize(h, s) kh_resize(ug, h, s)
#define csp_map_ug_put(h, k, r)  kh_put(ug, h, k, r)
#define csp_map_ug_get(h, k) kh_get(ug, h, k)
#define csp_map_ug_del(h, k) kh_del(ug, h, k)
#define csp_map_ug_exist(h, x) kh_exist(h, x)
#define csp_map_ug_key(h, x) kh_key(h, x)
#define csp_map_ug_val(h, x) kh_val(h, x)
#define csp_map_ug_begin(h) kh_begin(h)
#define csp_map_ug_end(h) kh_end(h)
#define csp_map_ug_size(h) kh_size(h)
#define csp_map_ug_reset(h) kh_clear(ug, h)
#define csp_map_ug_destroy(h) {								\
    if (h) {										\
        csp_map_ug_iter __k;								\
        for (__k = csp_map_ug_begin(h); __k != csp_map_ug_end(h); __k++) { 			\
            if (csp_map_ug_exist(h, __k)) csp_list_uu_destroy(csp_map_ug_val(h, __k));		\
        }											\
        kh_destroy(ug, h);									\
    }												\
}

/* Struct csp_list_qu_t APIs 
@abstract  The structure stores all qual value of one sample for certain query pos.
@param v   The csp_list_qu_t structure [csp_list_qu_t].

@example   Refer to the example in kvec.h.
 */
typedef kvec_t(int8_t) csp_list_qu_t;
#define csp_list_qu_init(v) kv_init(v)
#define csp_list_qu_resize(v, size) kv_resize(int8_t, v, size)
#define csp_list_qu_push(v, x) kv_push(int8_t, v, x)
#define csp_list_qu_A(v, i) kv_A(v, i)
#define csp_list_qu_size(v) kv_size(v)
#define csp_list_qu_max(v) kv_max(v)
#define csp_list_qu_destroy(v) kv_destroy(v)
#define csp_list_qu_reset(v) ((v).n = 0)

/*@abstract     Infer ref and alt based on the read count of each possible alleles (ACGTN).
@param bc       Pointer of array that stores the read count of each possible alleles, the read count is
                in the order of A/C/G/T/N.
@param ref_idx  Pointer of index of ref in "ACGTN".
@param alt_idx  Pointer of index of alt in "ACGTN".
@return         Void.

@note           It's usually called when the input pos has no ref or alt.
 */
static inline void csp_infer_allele(size_t *bc, int8_t *ref_idx, int8_t *alt_idx);

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
    csp_list_qu_t qu[5];
    double qmat[5][4];
    double gl[5];
    int ngl;
    csp_map_ug_t *hug;
} csp_plp_t;

/* note that the @p qu is also initialized after calling calloc(). */
static inline csp_plp_t* csp_plp_init(void); 
static inline void csp_plp_destroy(csp_plp_t *p); 
static inline void csp_plp_reset(csp_plp_t *p);

/*@abstract    Print the content to csp_plp_t to stream.
@param fp      Pointer of stream.
@param p       Pointer of csp_plpt_t to be printed.
@param prefix  Pointer of prefix string. Set to "" if no prefix.
@return        Void.
 */
static void csp_plp_print(FILE *fp, csp_plp_t *p, char *prefix);

/*@abstract     Format the content of csp_plp_t (the pileup stat info) of one cell/sample for certain query pos to 
                string in the output vcf file. 
@param p        Pointer of csp_plp_t structure corresponding to the pos.
@param s        Pointer of kstring_t which would store the formatted string.
@return         0 if success, -1 otherwise.
 */
static int csp_plp_str_vcf(csp_plp_t *p, kstring_t *s);

static int csp_plp_to_vcf(csp_plp_t *p, jfile_t *s);

/*@abstract    The HashMap maps sample-group-name (char*) to csp_plp_t (*).
@example       Refer to a simple example in khash.h.
 */
KHASH_MAP_INIT_STR(sg, csp_plp_t*)
typedef khash_t(sg) csp_map_sg_t;
#define csp_map_sg_iter khiter_t
#define csp_map_sg_init() kh_init(sg)
#define csp_map_sg_clear(h) kh_clear(sg, h)
#define csp_map_sg_resize(h, s) kh_resize(sg, h, s)
#define csp_map_sg_put(h, k, r)  kh_put(sg, h, k, r)
#define csp_map_sg_get(h, k) kh_get(sg, h, k)
#define csp_map_sg_del(h, k) kh_del(sg, h, k)
#define csp_map_sg_exist(h, x) kh_exist(h, x)
#define csp_map_sg_key(h, x) kh_key(h, x)
#define csp_map_sg_val(h, x) kh_val(h, x)
#define csp_map_sg_begin(h) kh_begin(h)
#define csp_map_sg_end(h) kh_end(h)
#define csp_map_sg_size(h) kh_size(h)
#define csp_map_sg_destroy(h) {								\
    if (h) {											\
        csp_map_sg_iter __k;									\
        for (__k = csp_map_sg_begin(h); __k != csp_map_sg_end(h); __k++) { 			\
            if (csp_map_sg_exist(h, __k)) csp_plp_destroy(csp_map_sg_val(h, __k)); 	\
        }										\
        kh_destroy(sg, h);									\
    }												\
}
#define csp_map_sg_reset_val(h) {								\
    if (h) {											\
        csp_map_sg_iter __k;									\
        for (__k = csp_map_sg_begin(h); __k != csp_map_sg_end(h); __k++) {			\
            if (csp_map_sg_exist(h, __k)) csp_plp_reset(csp_map_sg_val(h, __k)); 		\
        }											\
    }												\
}

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
@param hsg_iter Pointer of array of csp_map_sg_iter. The iter in the array is in the same order of sg names.
@param nsg   Size of csp_map_sg_iter array hsg_iter.
@param pu    Pool of csp_umi_unit_t structures.
@param pl    Pool of csp_list_uu_t structures.
@param su    Pool of UMI strings.
@param qvec  A container for the qual vector returned by get_qual_vector().
 */
typedef struct {
    int8_t ref_idx, alt_idx, inf_rid, inf_aid;
    size_t bc[5];
    size_t tc, ad, dp, oth;
    size_t nr_ad, nr_dp, nr_oth;
    csp_map_sg_t *hsg;
    csp_map_sg_iter *hsg_iter;
    int nsg;
    csp_pool_uu_t *pu;
    csp_pool_ul_t *pl;
    csp_pool_ps_t *su;
    double qvec[4];
} csp_mplp_t;

/*@abstract  Initialize the csp_mplp_t structure.
@return      Pointer to the csp_mplp_t structure if success, NULL otherwise.

@note        1. The kstring_t s is also initialized inside this function.   
             2. The valid pointer returned by this function should be freed by csp_mplp_destroy() function
                   when no longer used.
 */
static inline csp_mplp_t* csp_mplp_init(void); 
static inline void csp_mplp_destroy(csp_mplp_t *p); 
static inline void csp_mplp_reset(csp_mplp_t *p);

/*@abstract    Print the content to csp_mplp_t to stream.
@param fp      Pointer of stream.
@param p       Pointer of csp_mplpt_t to be printed.
@param prefix  Pointer of prefix string. Set to "" if no prefix.
@return        Void.
 */
static void csp_mplp_print(FILE *fp, csp_mplp_t *p, char *prefix);

static inline void csp_mplp_print_(FILE *fp, csp_mplp_t *p, char *prefix);

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
static int csp_mplp_set_sg(csp_mplp_t *p, char **s, const int n);

/*@abstract  Format the content of csp_mplp_t (the pileup stat info) of certain query pos to string in the output vcf file.
@param mplp  Pointer of the csp_mplp_t structure corresponding to the pos.
@param s     Pointer of kstring_t which stores the formatted string.
@return      0 if success, -1 otherwise.
 */
static inline int csp_mplp_str_vcf(csp_mplp_t *mplp, kstring_t *s);

static inline int csp_mplp_to_vcf(csp_mplp_t *mplp, jfile_t *s);

/*@abstract    Format the content of csp_mplp_t of certain query pos to string in the output sparse matrices file.
@param mplp    Pointer of the csp_mplp_t structure corresponding to the pos.
@param ks_ad   Pointer of kstring_t which is to store formatted AD string.
@param ks_dp   Pointer of kstring_t which is to store formatted DP string.
@param ks_oth  Pointer of kstring_t which is to store formatted OTH string.
@param idx     Index of the SNP/mplp (1-based).
@return        0 if success, -1 otherwise.
 */
static inline int csp_mplp_str_mtx(csp_mplp_t *mplp, kstring_t *ks_ad, kstring_t *ks_dp, kstring_t *ks_oth, size_t idx);

/*@abstract    Format the content of csp_mplp_t of certain query pos to string in the tmp output sparse matrices file.
@param mplp    Pointer of the csp_mplp_t structure corresponding to the pos.
@param ks_ad   Pointer of kstring_t which is to store formatted AD string.
@param ks_dp   Pointer of kstring_t which is to store formatted DP string.
@param ks_oth  Pointer of kstring_t which is to store formatted OTH string.
@return        0 if success, -1 otherwise.

@note          This function is used for tmp files.
 */
static inline int csp_mplp_str_mtx_tmp(csp_mplp_t *mplp, kstring_t *ks_ad, kstring_t *ks_dp, kstring_t *ks_oth);

static int csp_mplp_to_mtx(csp_mplp_t *mplp, jfile_t *fs_ad, jfile_t *fs_dp, jfile_t *fs_oth, size_t idx); 

#if DEVELOP
/* 
* Tags for sparse matrices
* It's useful when the input tags are optional.
*/
typedef size_t csp_mtx_value_t;
typedef int csp_mtx_iter_t;
typedef csp_mtx_value_t (*csp_mtx_value_func_t)(csp_mplp_t*, csp_plp_t*);

static inline csp_mtx_value_t csp_mtx_value_AD(csp_mplp_t *mplp, csp_plp_t *plp) {
    return plp->ad;
}

static inline csp_mtx_value_t csp_mtx_value_DP(csp_mplp_t *mplp, csp_plp_t *plp) {
    return plp->dp;
}

static inline csp_mtx_value_t csp_mtx_value_OTH(csp_mplp_t *mplp, csp_plp_t *plp) {
    return plp->oth;
}

static csp_mtx_iter_t csp_mtx_ntags = 3;
static char *csp_mtx_tags[] = {"AD", "DP", "OTH"};
static csp_mtx_value_func_t csp_mtx_value_funcs[] = { &csp_mtx_value_AD, &csp_mtx_value_DP, &csp_mtx_value_OTH };

static inline char* csp_get_mtx_fn(char *tag) {
    kstring_t ks = KS_INITIALIZE;
    ksprintf(&ks, "cellSNP.tag.%s.mtx", tag);
    char *p = strdup(ks_str(&ks));
    ks_free(&ks);
    return p;
}

static inline csp_mtx_iter_t csp_get_mtx_idx(char *tag) {
    csp_mtx_iter_t i;
    for (i = 0; i < csp_mtx_ntags; i++) {
        if (0 == strcmp(tag, csp_mtx_tags[i])) { return i; }
    } 
    return -1;
}

static inline csp_mtx_value_func_t csp_get_mtx_value_func(csp_mtx_iter_t i) { return csp_mtx_value_funcs[i]; }

typedef struct {
    char *out_fn;
    csp_mtx_value_func_t vfunc;
} csp_mtx_tag_fs;

static inline csp_mtx_tag_fs* csp_mtx_tag_fs_init(void) {
    return (csp_mtx_tag_fs*) calloc(1, sizeof(csp_mtx_tag_fs));
}

static inline void sp_mtx_tag_fs_destroy(csp_mtx_tag_fs *p) {
    if (p) { free(p->out_fn); free(p); }
}
#endif

#endif
