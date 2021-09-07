/* Global configure
 * Author: Xianjie Huang <hxj5@hku.hk>
 */
#ifndef CSP_CONFIG_H
#define CSP_CONFIG_H

#define DEBUG 0
#define VERBOSE 1
/* DEVELOP defined to 1 means some codes for future version of cellsnp will be included. */
#define DEVELOP 0

#define CSP_NAME "cellsnp-lite"
#define CSP_VERSION "1.2.1"
#define CSP_AUTHOR "hxj5<hxj5@hku.hk>"

#define JF_ZIP_TYPE         2            // use bgzip as zip method for JFile
#define CSP_FIT_MULTI_SMP   1            // if let nthread auto fit multi samples. 0, no; 1, yes.
#define TP_MAX_OPEN         1024         // default max number of open files

typedef struct _gll_settings global_settings;

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
// default max depth for one site of one bam file, 0 means highest possible value. used by bam_mplp_set_maxcnt().
#define CSP_PLP_MAX_DEPTH   0
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

#endif
