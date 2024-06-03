// jsam.h - sequence alignment file operations.


#ifndef SZ_JSAM_H
#define SZ_JSAM_H

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"

/* 
* BAM/SAM/CRAM File API 
*/

/*!@func
@abstract    Translate chr name to tid of bam/sam/cram.
@param hdr   Pointer of sam_hdr_t structure.
@param name  Chr name.
@param s     Pointer of kstring_t.
@return      Non-negative number if success, -1 for unrecognized reference seq name, -2 for unparsed hdr.
@note        If the translation failed, this function will try "<name> - chr" if name starts with "chr", try
             "chr + <name>" otherwise. 
 */
int csp_sam_hdr_name2id(sam_hdr_t *hdr, const char *name, kstring_t *s);

/*!@func
@abstract  Convert chrom name to be the same with the one in sam header. 
@return    Pointer to ref/chrom name if success, NULL otherwise.
@note      No need to free the returned char* pointer when success.
*/
const char* csp_fmt_chr_name(const char *name, sam_hdr_t *hdr, kstring_t *s);

/*!@func
@abstract  The two functions below convert raw cigar op/len value to real value.
@param c   Raw cigar op/len value stored in bam1_t, can be an element of cigar array obtained by bam_get_cigar(b) [uint32_t].
@return    An integer [int].
*/
#define get_cigar_op(c) ((c) & BAM_CIGAR_MASK)
#define get_cigar_len(c) ((c) >> BAM_CIGAR_SHIFT)

extern const char csp_nt5_str[5];

/*!@func
@abstract    Base Conversions
@param char  A char (A/C/G/T/N) [int8_t]
@param idx   Index in the seq_nt16_str for the base [int8_t].
@param int   Index in the 'ACGTN' for the base [int8_t].
 */
#define seq_nt16_idx2int(i) (seq_nt16_int[i])
#define seq_nt16_char2idx(c) (seq_nt16_table[c])
#define seq_nt16_char2int(c) (seq_nt16_idx2int(seq_nt16_char2idx(c)))
#define seq_nt16_int2char(i) (csp_nt5_str[i])

/*!@func
@abstract  Convert 4-bit integer to a letter. seq_nt16_str is declared in htslib/hts.h 
@param i   A 4-bit integer returned by bam_seqi(), standing for the index in the seq_nt16_str [int8_t].
@return    A letter standing for base char [int8_t].
 */
#define bam_seq_idx2base(i) (seq_nt16_str[i])

/*!@func
@abstract  Convert raw qual score stored in bam1_t to qual char.
@param i   Raw qual score stored in bam1_t [int8_t].
@return    Qual char [int8_t].
  */
#define bam_seq_qual2char(i) ((i) + 33)

/*!@func
@abstract   Get the content of bam aux tag of 'Z'/'H' (string) type.
@param b    Pointer to the bam1_t structure.
@param tag  Pointer to the bam aux tag.
@return     NULL if donot contain the tag or corrupted data, pointer to the content otherwise.
@note   1. To speed up, the caller should guarantee parameters b and tag are valid. 
        2. The data of the pointer returned by this function is part of bam1_t, so do not double free!
 */
char* get_bam_aux_str(bam1_t *b, const char tag[2]);

#endif