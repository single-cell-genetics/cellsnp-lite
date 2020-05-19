
#ifndef CELL_SNP_UTIL_H
#define CELL_SNP_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "kvec.h"
#include "general_util.h"

typedef struct _gll_settings global_settings;

/* 
* BAM/SAM/CRAM File API 
*/

/*@abstract  Translate chr name to tid of bam/sam/cram.
@param hdr   Pointer of sam_hdr_t structure.
@param name  Chr name.
@param s     Pointer of kstring_t.
@return      Non-negative number if success, -1 for unrecognized reference seq name, -2 for unparsed hdr.

@note        If the translation failed, this function will try "<name> - chr" if name starts with "chr", try
             "chr + <name>" otherwise. 
 */
static inline int csp_sam_hdr_name2id(sam_hdr_t *hdr, const char *name, kstring_t *s) {
    int tid;
    if ((tid = sam_hdr_name2tid(hdr, name)) < 0) {
        if (tid < -1) { return tid; }
        if (strlen(name) > 3 && 'c' == name[0] && 'h' == name[1] && 'r' == name[2]) { return sam_hdr_name2tid(hdr, name + 3); }
        else {
            kputs("chr", s); kputs(name, s);
            tid = sam_hdr_name2tid(hdr, ks_str(s));
            ks_clear(s);
            return tid;
        }
    } else { return tid; }
}

/*@abstract   The two functions below convert raw cigar op/len value to real value.
@param c      Raw cigar op/len value stored in bam1_t, can be an element of cigar array obtained by bam_get_cigar(b) [uint32_t].
@return       An integer [int].
*/
#define get_cigar_op(c) ((c) & BAM_CIGAR_MASK)
#define get_cigar_len(c) ((c) >> BAM_CIGAR_SHIFT)

/*@abstract  Convert index of seq_nt16_str to the index (0-4) of A/C/G/T/G.
@param i     Index in the seq_nt16_str for the base [int8_t].
@return      Index in the 'ACGTN' for the base [int8_t].
 */
#define seq_nt16_idx2int(i) (seq_nt16_int[i])

/*@abstract  Convert a char (A/C/G/T/N) to index of seq_nt16_str/seq_nt16_int.
@param c     The base char [int8_t].
@return      Index in the seq_nt16_str for the base [int8_t].
  */
#define seq_nt16_char2idx(c) (seq_nt16_table[c])

/*@abstract  Convert a char (A/C/G/T/N) to the index (0-4) of A/C/G/T/G.
@param c     The base char [int8_t].
@return      Index in the 'ACGTN' for the base [int8_t].
  */
#define seq_nt16_char2int(c) (seq_nt16_idx2int(seq_nt16_char2idx(c)))

static char csp_nt5_str[] = "ACGTN";

/*@abstract  Convert index in "ACGTN" to a letter.
@param i     Index in "ACGTN"
@return      A letter.
 */
#define seq_nt16_int2char(i) (csp_nt5_str[i])

/*@abstract  Convert 4-bit integer to a letter. seq_nt16_str is declared in htslib/hts.h 
@param i     A 4-bit integer returned by bam_seqi(), standing for the index in the seq_nt16_str [int8_t].
@return      A letter standing for base char [int8_t].
 */
#define bam_seq_idx2base(i) (seq_nt16_str[i])

/*@abstract  Convert raw qual score stored in bam1_t to qual char.
@param i     Raw qual score stored in bam1_t [int8_t].
@return      Qual char [int8_t].
  */
#define bam_seq_qual2char(i) ((i) + 33)

/*@abstract   Get the content of bam aux tag of 'Z'/'H' (string) type.
@param b      Pointer to the bam1_t structure.
@param tag    Pointer to the bam aux tag.
@return       NULL if donot contain the tag or corrupted data, pointer to the content otherwise.

@note   1. To speed up, the caller should guarantee parameters b and tag are valid. 
        2. The data of the pointer returned by this function is part of bam1_t, so do not double free!
 */
static inline char* get_bam_aux_str(bam1_t *b, const char tag[2]) {
    uint8_t *data;
    if (NULL == (data = bam_aux_get(b, tag))) { return NULL; }
    return bam_aux2Z(data);
}

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
static inline csp_bam_fs* csp_bam_fs_init(void) { return (csp_bam_fs*) calloc(1, sizeof(csp_bam_fs)); }

static inline void csp_bam_fs_destroy(csp_bam_fs *p) {
    if (p) {
        if (p->idx) { hts_idx_destroy(p->idx); }
        if (p->hdr) { sam_hdr_destroy(p->hdr); }
        if (p->fp)  { hts_close(p->fp); }
        free(p);
    }
}

/*@abstract  Build the csp_bam_fs structure.  
@param fn    Filename of the bam/sam/cram.
@param ret   Pointer to the state.
             0 if success, negative numbers otherwise:
               -1, fn is NULL.
               -3, init csp_bam_fs structure error.
               -5, open bam/sam/cram error.
               -7, read sam header error.
               -9, load sam index error.
@return      The pointer to the csp_bam_fs structure if success, NULL otherwise.

@note        The pointer returned successfully by csp_bam_fs_build() should be freed
             by csp_bam_fs_destroy() when no longer used.
*/
static inline csp_bam_fs* csp_bam_fs_build(const char *fn, int *ret) {
    csp_bam_fs *p;
    if (NULL == fn) { *ret = -1; return NULL; }
    if (NULL == (p = csp_bam_fs_init())) { *ret = -3; return NULL; }
    if (NULL == (p->fp = hts_open(fn, "rb"))) { *ret = -5; goto fail; }
    if (NULL == (p->hdr = sam_hdr_read(p->fp))) { *ret = -7; goto fail; }
    if (NULL == (p->idx = sam_index_load(p->fp, fn))) { *ret = -9; goto fail; }
    *ret = 0;
    return p;
  fail:
    csp_bam_fs_destroy(p);
    return NULL;		
}

static inline int csp_bam_fs_reset(csp_bam_fs *p, const char *fn) {
    if (NULL == p) { return -1; }
    if (p->idx) { hts_idx_destroy(p->idx); }
    if (p->hdr) { sam_hdr_destroy(p->hdr); }
    if (p->fp)  { hts_close(p->fp); }		
    if (NULL == (p->fp = hts_open(fn, "rb"))) { return -3; }
    if (NULL == (p->hdr = sam_hdr_read(p->fp))) { return -5; }
    if (NULL == (p->idx = sam_index_load(p->fp, fn))) { return -7; }
    return 0;
}

/* 
* SNP List API
*/
/*@abstruct    The basic data structure to store SNP info.
@param chr     Name of chr.
@param pos     0-based coordinate in the reference sequence.
@param ref     Ref base (a letter). 0 means no ref in the input SNP file for the pos.
@param alt     Alt base (a letter). 0 means no alt in the input SNP file for the pos.
 */
typedef struct {
    char *chr;   
    hts_pos_t pos; 
    int8_t ref, alt;
} csp_snp_t;

/*@abstract  Initilize the csp_snp_t structure.
@return      Pointer to the csp_snp_t structure if success, NULL if failure.
@note        The pointer returned successfully by csp_snp_init() should be freed
             by csp_snp_destroy() when no longer used.
 */
static inline csp_snp_t* csp_snp_init(void) { return (csp_snp_t*) calloc(1, sizeof(csp_snp_t)); }

static inline void csp_snp_destroy(csp_snp_t *p) { 
    if (p) { free(p->chr); free(p); } 
}

/*@abstract  A list containing of several pointers to the csp_snp_t structure.
@param a     Array of csp_snp_t* pointers.
@param n     The next pos of the unused element of the array.
@param m     Size of the whole array.

@note        The csp_snplist_t structure should be freed by csp_snplist_destroy() when no longer used.
             csp_snplist_init() function should be called immediately after the structure was created.

@example (A simple example from kvec.h)
    kvec_t(int) array;
    kv_init(array);
    kv_push(int, array, 10); // append
    kv_A(array, 20) = 4;
    kv_destroy(array);
 */
typedef kvec_t(csp_snp_t*) csp_snplist_t;   /* kvec_t from kvec.h */
#define csp_snplist_init(v) kv_init(v)
#define csp_snplist_resize(v, size) kv_resize(csp_snp_t*, v, size)
#define csp_snplist_push(v, x) kv_push(csp_snp_t*, v, x)
#define csp_snplist_A(v, i) kv_A(v, i)
#define csp_snplist_size(v) kv_size(v)
#define csp_snplist_max(v) kv_max(v)
#define csp_snplist_destroy(v) {																		\
    size_t __j;																							\
    for (__j = 0; __j < csp_snplist_size(v); __j++) csp_snp_destroy(csp_snplist_A(v, __j));					\
    kv_destroy(v);																						\
}

/*@abstract    Extract SNP info from bcf/vcf file.
@param fn      Filename of bcf/vcf.
@param pl      Pointer to client data used to store the extracted SNP info.
@param ret     Pointer to store running state. 0 if success, -1 otherwise.
@return        Num of elements successfully added into the snplist.

@note          If length of Ref or Alt is larger than 1, then the SNP would be skipped.
               If length of Ref or Alt is 0, then their values would be infered during pileup.
 */
static size_t get_snplist_from_vcf(const char *fn, csp_snplist_t *pl, int *ret) {
    htsFile *fp = NULL;
    bcf_hdr_t *hdr = NULL;
    bcf1_t *rec = NULL;
    csp_snp_t *ip = NULL;
    size_t l, n = 0;  /* use n instead of kv_size(*pl) in case that pl is not empty at the beginning. */
    int r;
    *ret = -1;
    if (NULL == fn || NULL == pl) { return 0; }
    if (NULL == (fp = hts_open(fn, "rb"))) { fprintf(stderr, "[E::%s] could not open '%s'\n", __func__, fn); return 0; }
    if (NULL == (hdr = bcf_hdr_read(fp))) { fprintf(stderr, "[E::%s] could not read header for '%s'\n", __func__, fn); goto fail; }
    if (NULL == (rec = bcf_init1())) { fprintf(stderr, "[E::%s] could not initialize the bcf structure.\n", __func__); goto fail; }
    while ((r = bcf_read1(fp, hdr, rec)) >= 0) {
        if (NULL == (ip = csp_snp_init())) { 
            fprintf(stderr, "[E::%s] could not initialize the csp_snp_t structure.\n", __func__); 
            goto fail; 
        }
        ip->chr = safe_strdup(bcf_hdr_id2name(hdr, rec->rid));
        if (NULL == ip->chr) {
            fprintf(stderr, "[W::%s] could not get chr name for rid = %d, pos = %ld (both 0-based).\n", __func__, rec->rid, rec->pos); 
            continue;
        } ip->pos = rec->pos;
        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele) {
            if (1 == (l = strlen(rec->d.allele[0]))) { ip->ref = rec->d.allele[0][0]; }
            else if (l > 1) { 
                fprintf(stderr, "[W::%s] skip No.%ld SNP: chr = %s, pos = %ld, ref = %s\n", __func__, n + 1, ip->chr, ip->pos, rec->d.allele[0]); 
                csp_snp_destroy(ip); 
                continue; 
            } // else: do nothing. keep ref = 0, 0 is its init value.
            if (rec->n_allele > 1) {
                if (1 == (l = strlen(rec->d.allele[1]))) { ip->alt = rec->d.allele[1][0]; }
                else if (l > 1) {
                    fprintf(stderr, "[W::%s] skip No.%ld SNP: chr = %s, pos = %ld, alt = %s\n", __func__, n + 1, ip->chr, ip->pos, rec->d.allele[1]); 
                    csp_snp_destroy(ip); 
                    continue; 					
                } // else: do nothing. keep alt = 0, 0 is its init value.
            }
        }
        csp_snplist_push(*pl, ip);
        n++;
    }
    if (-1 == r) { // end of bcf file.
        csp_snplist_resize(*pl, csp_snplist_size(*pl));
    } else { 
        fprintf(stderr, "[E::%s] error when parsing '%s'\n", __func__, fn); 
        goto fail; 
    }
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);
    *ret = 0;
    return n;
  fail:
    if (rec) bcf_destroy1(rec);
    if (hdr) bcf_hdr_destroy(hdr);
    if (fp) hts_close(fp);
    return n;
}
#define get_snplist(fn, pl, ret) get_snplist_from_vcf((fn), (pl), (ret))

/* 
* Thread API 
*/
/*@abstract    The data structure used as thread-func parameter.
@param gs      Pointer to the global_settings structure.
@param n       Pos of next element in the snplist to be used by certain thread.
@param m       Total size of elements to be used by certain thread, must not be changed.
@param ci      Index of the element in the chrom_all array used by certain thread.
@param i       Id of the thread data.
@param ret     Running state of the thread.
@param out_fn  Name of the output file. Data that the pointer points to would be freed elsewhere.
@param out_fm  Mode of the output file. w/a.
@param is_out_zip  If out file needs to be zipped.
 */
typedef struct {
    global_settings *gs;
    size_t m, n;   // for snplist.
    int ci;        // for chromosome(s).
    int i;
    int ret;
    char *out_fn;
    char *out_fm;
    int is_out_zip;
} thread_data;

/*@abstract  Create the thread_data structure.
@return      Pointer to the structure if success, NULL otherwise.
@note        The pointer returned successfully by thdata_init() should be freed
             by thdata_destroy() when no longer used.
 */
static inline thread_data* thdata_init(void) { return (thread_data*) calloc(1, sizeof(thread_data)); }

static inline void thdata_destroy(thread_data *p) { free(p); }

static inline void thdata_print(FILE *fp, thread_data *p) {
    fprintf(fp, "\tm = %ld, n = %ld\n", p->m, p->n);
    fprintf(fp, "\tci = %d, i = %d, ret = %d\n", p->ci, p->i, p->ret);
    fprintf(fp, "\tout_fn = '%s', out_fm = '%s', is_out_zip = %d\n", p->out_fn, p->out_fm, p->is_out_zip);
}

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
static inline csp_pileup_t* csp_pileup_init(void) {
    csp_pileup_t *p = (csp_pileup_t*) malloc(sizeof(csp_pileup_t));
    if (p) {
        if (NULL == (p->b = bam_init1())) { free(p); return NULL; }
    }
    return p;
}

static inline void csp_pileup_destroy(csp_pileup_t *p) { 
    if (p) {
        if (p->b) bam_destroy1(p->b);	
        free(p);
    } 
}

/* reset the csp_pileup_t structure without reallocating memory.
   return 0 if success, -1 otherwise. */
static inline int csp_pileup_reset(csp_pileup_t *p) {
    if (p) {
        if (p->b) { bam_destroy1(p->b); }
        memset(p, 0, sizeof(csp_pileup_t)); 
        if (NULL == (p->b = bam_init1())) { return -1; }
    } // else: do nothing.
    return 0;
}

/* only reset part of the csp_pileup_t as values of parameters in other parts will be immediately overwritten after
   calling this function. It's often called by pileup_read_with_fetch(). */
static inline void csp_pileup_reset_(csp_pileup_t *p) { }

static inline void csp_pileup_print(FILE *fp, csp_pileup_t *p) {
    fprintf(fp, "qpos = %d\n", p->qpos);
    fprintf(fp, "base = %c, qual = %d\n", p->base, p->qual);
    fprintf(fp, "is_refskip = %d, is_del = %d\n", p->is_refskip, p->is_del);
    fprintf(fp, "umi = %s, cb = %s\n", p->umi, p->cb);
    fprintf(fp, "len_aln = %d\n", p->laln);
}

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

static inline csp_umi_unit_t* csp_umi_unit_init(void) {
    csp_umi_unit_t *p = (csp_umi_unit_t*) calloc(1, sizeof(csp_umi_unit_t));
    return p;   /* will set values just after this function is called so no need to set init values here. */
}
#define csp_umi_unit_reset(p)
static inline void csp_umi_unit_destroy(csp_umi_unit_t *p) { free(p); }

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
#define csp_list_uu_init(v) kv_init(*(v))
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
#define csp_map_ug_destroy(h) {																									\
    if (h) {																														\
        csp_map_ug_iter __k;																									\
        for (__k = csp_map_ug_begin(h); __k != csp_map_ug_end(h); __k++) { 														\
            if (csp_map_ug_exist(h, __k)) csp_list_uu_destroy(csp_map_ug_val(h, __k));																\
        }																													\
        kh_destroy(ug, h);																										\
    }																															\
}

SZ_NUMERIC_OP_INIT(cu_d, double)
SZ_NUMERIC_OP_INIT(cu_s, size_t)

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
static inline int get_qual_vector(double qual, double cap_bq, double min_bq, double *rv) {
    double bq = max2(min2(cap_bq, qual), min_bq);
    double p = pow(0.1, bq / 10);
    rv[0] = log(1 - p);              rv[1] = log(0.75 - 2.0 / 3 * p);
    rv[2] = log(0.5 - 1.0 / 3 * p);  rv[3] = log(p);
    return 0;
}

/*@abstract     Internal function to translate qual matrix to vector of GL.
@param qm       Qual matrix: 5-by-4, columns are [1-Q, 3/4-2/3Q, 1/2-1/3Q, Q]. See csp_plp_t::qmat
@param bc       Base count: (5, ). See csp_plp_t::bcount
@param ref_idx  Index of ref base in 'ACGTN'. 
@param alt_idx  Index of alt base in 'ACGTN'.
@param db       doublet_GL. 1/0.
@param gl       Array of GL: loglikelihood for 
                     GL1: L(rr|qual_matrix, base_count), 
                     GL2-GL5: L(ra|..), L(aa|..), L(rr+ra|..), L(ra+aa|..)
                Size must be no less than 5. Be careful that no size check inside the function
@param n        Pointer of num of elements in GL array.
@return         0 if success, -1 otherwise.

@note           In some special cases, ref=A and alt=AG for example, the ref_idx would be equal with alt_idx.
                Should be fixed in future.
 */
static int qual_matrix_to_geno(double qm[][4], size_t *bc, int8_t ref_idx, int8_t alt_idx, int db, double *gl, int *n) {
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

/*@abstract     Infer ref and alt based on the read count of each possible alleles (ACGTN).
@param bc       Pointer of array that stores the read count of each possible alleles, the read count is
                in the order of A/C/G/T/N.
@param ref_idx  Pointer of index of ref in "ACGTN".
@param alt_idx  Pointer of index of alt in "ACGTN".
@return         Void.

@note           It's usually called when the input pos has no ref or alt.
 */
static inline void csp_infer_allele(size_t *bc, int8_t *ref_idx, int8_t *alt_idx) {
    int8_t i, k1, k2;
    size_t m1, m2;
    if (bc[0] < bc[1]) { m1 = bc[1]; m2 = bc[0]; k1 = 1; k2 = 0; }
    else { m1 = bc[0]; m2 = bc[1]; k1 = 0; k2 = 1; }
    for (i = 2; i < 5; i++) {
        if (bc[i] > m1) { m1 = bc[i]; k1 = i; }
        else if (bc[i] > m2) { m2 = bc[i]; k2 = i; }
    }
    *ref_idx = k1; *alt_idx = k2;
}

/*@abstract      The structure that store pileup info of one cell/sample for certain query pos.
@param bcount    Total read count for each base in the order of 'ACGTN'.
@param tcount    Total read count for all bases.
@param qmat      Matrix of qual with 'ACGTN' vs. [1-Q, 3/4-2/3Q, 1/2-1/3Q, Q].
@param gl        Array of GL: loglikelihood for 
                     GL1: L(rr|qual_matrix, base_count), 
                     GL2-GL5: L(ra|..), L(aa|..), L(rr+ra|..), L(ra+aa|..)
                 Note that the right values of its elements would be altered and should not be used anymore, 
                 after calling csp_plp_to_vcf(). 
@param ngl       Num of valid elements in the array gl.
@param h         Pointer of hash table that stores stat info of UMI groups.
 */
typedef struct {
    size_t bcount[5];
    size_t tcount;
    double qmat[5][4];
    double gl[5];
    int ngl;
    csp_map_ug_t *h;
} csp_plp_t;

static inline csp_plp_t* csp_plp_init(void) { return (csp_plp_t*) calloc(1, sizeof(csp_plp_t)); }

static inline void csp_plp_destroy(csp_plp_t *p) { 
    if (p) { 
        if (p->h) { csp_map_ug_destroy(p->h); }
        free(p); 
    }
}

static inline void csp_plp_reset(csp_plp_t *p) {
    if (p) {
        memset(p->bcount, 0, sizeof(p->bcount));
        p->tcount = 0;
        memset(p->qmat, 0, sizeof(p->qmat));
        p->ngl = 0;
        if (p->h) { csp_map_ug_reset(p->h); } 
    }
}

/*@abstract    Print the content to csp_plp_t to stream.
@param fp      Pointer of stream.
@param p       Pointer of csp_plpt_t to be printed.
@param prefix  Pointer of prefix string. Set to "" if no prefix.
@return        Void.
 */
static void csp_plp_print(FILE *fp, csp_plp_t *p, char *prefix) {
    int i, j;
    csp_map_ug_iter u;
    fprintf(fp, "%stotal read count = %ld\n", prefix, p->tcount);
    fprintf(fp, "%sbase count (A/C/G/T/N):", prefix);
    for (i = 0; i < 5; i++) fprintf(fp, " %ld", p->bcount[i]);
    fputc('\n', fp);
    fprintf(fp, "%squal matrix 5x4:\n", prefix);
    for (i = 0; i < 5; i++) {
        fprintf(fp, "%s\t", prefix);
        for (j = 0; j < 4; j++) fprintf(fp, " %.2f", p->qmat[i][j]);
        fputc('\n', fp);
    }
    fprintf(fp, "%snum of geno likelihood = %d\n", prefix, p->ngl);
    if (p->ngl) {
        fprintf(fp, "%sgeno likelihood:", prefix);
        for (i = 0; i < p->ngl; i++) fprintf(fp, " %.2f", p->gl[i]);
        fputc('\n', fp);
    }
    if (p->h) {
        int size = csp_map_ug_size(p->h);
        fprintf(fp, "%ssize of the csp_map_ug = %d\n", prefix, size);
        if (size) {
            fprintf(fp, "%s", prefix);
            for (u = csp_map_ug_begin(p->h); u != csp_map_ug_end(p->h); u++) {
                if (csp_map_ug_exist(p->h, u)) { fprintf(fp, " %s", csp_map_ug_key(p->h, u)); }
            } fputc('\n', fp);
        }
    }
}

/*@abstract     Format the content of csp_plp_t (the pileup stat info) of one cell/sample for certain query pos to 
                string in the output vcf file. 
@param p        Pointer of csp_plp_t structure corresponding to the pos.
@param s        Pointer of kstring_t which would store the formatted string.
@param ref_idx  Index of ref base in 'ACGTN'.
@param alt_idx  Index of alt base in 'ACGTN'.
@return         0 if success, -1 otherwise.

@note           The right values in plp->gl are changed after calling this function, so do not use the values in
                plp->gl since then.
 */
static int csp_plp_to_vcf(csp_plp_t *p, kstring_t *s, int8_t ref_idx, int8_t alt_idx) {
    if (p->tcount <= 0) { kputs(".:.:.:.:.:.", s); return 0; }
    size_t ref_cnt, alt_cnt, dp_cnt, oth_cnt;
    int i, m;
    double tmp = -10 / log(10);
    char *gt[] = {"0/0", "1/0", "1/1"};
    ref_cnt = p->bcount[ref_idx];   alt_cnt = p->bcount[alt_idx];
    dp_cnt = ref_cnt + alt_cnt;     oth_cnt = p->tcount - dp_cnt;
    m = get_idx_of_max(cu_d, p->gl, 3);
    kputs(gt[m], s);
    ksprintf(s, ":%ld:%ld:%ld:", alt_cnt, dp_cnt, oth_cnt);
    for (i = 0; i < p->ngl; i++) { p->gl[i] *= tmp; }
    if (join_arr_to_str(cu_d, p->gl, p->ngl, ',', "%.0f", s) < p->ngl) { return -1; }
    kputc_(':', s);
    if (join_arr_to_str(cu_s, p->bcount, 5, ',', "%ld", s) < 5) { return -1; }  // here first s is name of SZ_NUMERIC_OP while the latter one is pointer of kstring_t.
    return 0;
}

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
#define csp_map_sg_destroy(h) {																								\
    if (h) {																												\
        csp_map_sg_iter __k;																									\
        for (__k = csp_map_sg_begin(h); __k != csp_map_sg_end(h); __k++) { 													\
            if (csp_map_sg_exist(h, __k)) csp_plp_destroy(csp_map_sg_val(h, __k)); 												\
        }																													\
        kh_destroy(sg, h);																								\
    }																														\
}
#define csp_map_sg_reset_val(h) {																							\
    if (h) {																												\
        csp_map_sg_iter __k;																								\
        for (__k = csp_map_sg_begin(h); __k != csp_map_sg_end(h); __k++) {												\
            if (csp_map_sg_exist(h, __k)) csp_plp_reset(csp_map_sg_val(h, __k)); 												\
        }																													\
    }																														\
}

/*@abstract  The structure stores the stat info of all sample groups for certain query pos.
@param ref_idx  Index of ref in "ACGTN". Negative number means not valid value.
@param alt_idx  Index of alt in "ACGTN". Negative number means not valid value.
@param inf_rid  Infered index of ref in "ACGTN". Negative number means not valid value.
@param inf_aid  Infered index of alt in "ACGTN". Negative number means not valid value.
@param bc    Read count of each base summarizing all sample groups for the pos, in the order of 'ACGTN'.
@param tc    Total read count of all bases for the pos.
@param h     HashMap that stores the stat info of all sample groups for the pos.
@param hiter Pointer of array of csp_map_sg_iter. The iter in the array is in the same order of sg names.
@param nsg   Size of csp_map_sg_iter array hiter.
@param pu    Pool of csp_umi_unit_t structures.
@param pl    Pool of csp_list_uu_t structures.
@param su    Pool of UMI strings.
@param s     The kstring_t structure to store the formatted string of csp_mplp_t. The string would be output to vcf.
@param qvec  A container for the qual vector returned by get_qual_vector().
 */
typedef struct {
    int8_t ref_idx, alt_idx, inf_rid, inf_aid;
    size_t bc[5];
    size_t tc;
    csp_map_sg_t *h;
    csp_map_sg_iter *hiter;
    int nsg;
    csp_pool_uu_t *pu;
    csp_pool_ul_t *pl;
    csp_pool_ps_t *su;
    kstring_t s;
    double qvec[4];
} csp_mplp_t;

/*@abstract  Initialize the csp_mplp_t structure.
@return      Pointer to the csp_mplp_t structure if success, NULL otherwise.

@note        1. The kstring_t s is also initialized inside this function.   
             2. The valid pointer returned by this function should be freed by csp_mplp_destroy() function
                   when no longer used.
 */
static inline csp_mplp_t* csp_mplp_init(void) { 
    csp_mplp_t *p = (csp_mplp_t*) calloc(1, sizeof(csp_mplp_t));
    return p;
}

static inline void csp_mplp_destroy(csp_mplp_t *p) { 
    if (p) {
        if (p->h) { csp_map_sg_destroy(p->h); }
        if (p->hiter) { free(p->hiter); }
        if (p->pu) { csp_pool_uu_destroy(p->pu); }
        if (p->pl) { csp_pool_ul_destroy(p->pl); }
        if (p->su) { csp_pool_ps_destroy(p->su); }
        ks_free(&p->s);
        free(p); 
    }
}

static inline void csp_mplp_reset(csp_mplp_t *p) {
    if (p) {
        memset(p->bc, 0, sizeof(p->bc));
        p->tc = 0;
        if (p->h) { csp_map_sg_reset_val(p->h); }
        if (p->pu) { csp_pool_uu_reset(p->pu); }
        if (p->pl) { csp_pool_ul_reset(p->pl); }
        if (p->su) { csp_pool_ps_reset(p->su); }
        memset(p->qvec, 0, sizeof(p->qvec));
        ks_clear(&p->s);
    }
}

/*@abstract    Print the content to csp_mplp_t to stream.
@param fp      Pointer of stream.
@param p       Pointer of csp_mplpt_t to be printed.
@param prefix  Pointer of prefix string. Set to "" if no prefix.
 */
static void csp_mplp_print(FILE *fp, csp_mplp_t *p, char *prefix) {
    int i;
    csp_plp_t *plp;
    kstring_t ks = KS_INITIALIZE;
    kstring_t *s = &ks;
    fprintf(fp, "%sref_idx = %d, alt_idx = %d\n", prefix, p->ref_idx, p->alt_idx);
    fprintf(fp, "%sinf_rid = %d, inf_aid = %d\n", prefix, p->inf_rid, p->inf_aid);
    fprintf(fp, "%stotal base count = %ld\n", prefix, p->tc);
    fprintf(fp, "%sbase count (A/C/G/T/N):", prefix);
    for (i = 0; i < 5; i++) { fprintf(fp, " %ld", p->bc[i]); }
    fputc('\n', fp);
    fprintf(fp, "%snum of sample group = %d\n", prefix, p->nsg);
    if (p->nsg) {
        kputs(prefix, s); kputc('\t', s);
        for (i = 0; i < p->nsg; i++) {
            fprintf(fp, "%sSG-%d = %s:\n", prefix, i, csp_map_sg_key(p->h, p->hiter[i]));
            plp = csp_map_sg_val(p->h, p->hiter[i]);
            csp_plp_print(fp, plp, ks_str(s));
        }
    }
    ks_free(s);
}

static void csp_mplp_print_(FILE *fp, csp_mplp_t *p, char *prefix) {
    int i;
    fprintf(fp, "%sref_idx = %d, alt_idx = %d\n", prefix, p->ref_idx, p->alt_idx);
    fprintf(fp, "%sinf_rid = %d, inf_aid = %d\n", prefix, p->inf_rid, p->inf_aid);
    fprintf(fp, "%stotal base count = %ld\n", prefix, p->tc);
    fprintf(fp, "%sbase count (A/C/G/T/N):", prefix);
    for (i = 0; i < 5; i++) { fprintf(fp, " %ld", p->bc[i]); }
    fputc('\n', fp);
    fprintf(fp, "%snum of sample group = %d\n", prefix, p->nsg);
}

/*@abstract  Set sample group names for the csp_mplp_t structure.
@param p     Pointer to the csp_mplp_t structure.
@parma s     Pointer to the array of names of sample groups.
@param n     Num of sample groups.
@return      0, no error; negative numbers otherwise:
               -1, invalid parameter;
               -3, failed to allocate memory for sg HashMap.
               -5, failed to allocate memory for csp_map_sg_iter array.
               -7, kh_put failed or have repeted sg names;
               -9, invalid sg names, i.e., NULL;
               -11, HashMap error.

@note        1. This function should be called just one time right after csp_mplp_t structure was created
                becuase the sgname wouldn't change once set.
             2. The HashMap (for sgnames) in csp_mplp_t should be empty or NULL.
 */
static int csp_mplp_set_sg(csp_mplp_t *p, char **s, const int n) {
    if (NULL == p || NULL == s || 0 == n) { return -1; }
    int i, r;
    csp_map_sg_iter k;
    if (NULL == p->h && NULL == (p->h = csp_map_sg_init())) { return -3; }
    if (NULL == p->hiter && NULL == (p->hiter = (csp_map_sg_iter*) malloc(sizeof(csp_map_sg_iter) * n))) { return -5; }
    for (i = 0; i < n; i++) {
        if (s[i]) { 
            k = csp_map_sg_put(p->h, s[i], &r); 
            if (r <= 0) { csp_mplp_destroy(p); return -7; } /* r = 0 means repeatd sgnames. */
            else { csp_map_sg_val(p->h, k) = NULL; }
        } else { csp_mplp_destroy(p); return -9; }
    }
    /* Storing iter index for each sg (sample group) name must be done after all sg names have been pushed into 
       the HashMap in case that the internal arrays of HashMap autoly shrink or some else modifications. */
    for (i = 0; i < n; i++) {
        k = csp_map_sg_get(p->h, s[i]);
        if (k == csp_map_sg_end(p->h)) { return -11; }
        else { p->hiter[i] = k; }
    }
    p->nsg = n;
    return 0;
}

/*@abstract  Format the content of csp_mplp_t (the pileup stat info) of certain query pos to string in the output vcf file.
@param mplp  Pointer of the csp_mplp_t structure corresponding to the pos.
@return      0 if success, -1 otherwise.
 */
static int csp_mplp_to_vcf(csp_mplp_t *mplp) {
    kstring_t *s = &mplp->s;
    int i;
    for (i = 0; i < mplp->nsg; i++) {
        kputc_('\t', s);
        if (csp_plp_to_vcf(csp_map_sg_val(mplp->h, mplp->hiter[i]), s, mplp->ref_idx, mplp->alt_idx) < 0) { return -1; }
    } //s->s[--(s->l)] = '\0';    /* s->l could not be negative unless no csp_plp_t(s) are printed to s->s. */
    return 0;
}

/*
* File API
 */
#endif