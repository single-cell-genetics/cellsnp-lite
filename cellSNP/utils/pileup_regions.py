# Utilility functions for pileup SNPs across regions 
# originally in from pileup_utils.py
# Author: Yuanhua Huang
# Date: 21/05/2018

from .pileup_utils import *

## ealier high error in pileup whole genome might come from
## using _read.query_sequence, which has only partially aligned
## pileupread.query_position is based on the full length of the reads
## _read.qqual is based on aligned reads segment

# def pileup_bases(pileupColumn):
#     """ pileup all reads mapped to the genome position.
#     """
#     base_list, read_list, qual_list = [], [], []
#     for pileupread in pileupColumn.pileups:
#         # query position is None if is_del or is_refskip is set.
#         if pileupread.is_del or pileupread.is_refskip:
#             continue
#         #query_POS = pileupread.query_position
#         query_POS = pileupread.query_position
#         _read = pileupread.alignment
#         try:
#             _base = _read.query_sequence[query_POS - 1].upper()
#             _qual = _read.qqual[query_POS - 1]
#         except:
#             print("warnings: a read fails to give _base or _qual.", 
#                   query_POS, len(_read.qqual), len(_read.qual), len(_read.query_sequence))
#             print(_read.qqual)
#             continue
#         #_qual = "J"
#         read_list.append(_read)
#         base_list.append(_base)
#         qual_list.append(_qual)
#     return base_list, qual_list, read_list




def pileup_bases(pileupColumn, real_POS=None):
    """ pileup all reads mapped to the genome position.
    """
    base_list, read_list, qual_list = [], [], []
    for pileupread in pileupColumn.pileups:
        # query position is None if is_del or is_refskip is set.
        if pileupread.is_del or pileupread.is_refskip:
            continue
            
        _read = pileupread.alignment
        if real_POS is not None:
            try:
                idx = _read.positions.index(real_POS-1)
            except:
                continue
            _qual = _read.qqual[idx]
            _base = _read.query_alignment_sequence[idx].upper()
        else:
            query_POS = pileupread.query_position
            _qual = _read.qual[query_POS - 1]
            _base = _read.query_sequence[query_POS - 1].upper()
            
        read_list.append(_read)
        base_list.append(_base)
        qual_list.append(_qual)
    return base_list, qual_list, read_list


def pileup_regions(samFile, barcodes, out_file=None, chrom=None, cell_tag="CR", 
                   UMI_tag="UR", min_COUNT=20, min_MAF=0.1, min_MAPQ=20, 
                   max_FLAG=255, min_LEN=30, doublet_GL=False, verbose=True):
    """Pileup allelic specific expression for a whole chromosome in sam file.
    TODO: 1) multiple sam files, e.g., bulk samples; 2) optional cell barcode
    """
    samFile, chrom = check_pysam_chrom(samFile, chrom)
    if out_file is not None:
        fid = open(out_file, "w")
        fid.writelines(VCF_HEADER + CONTIG)
        fid.writelines("\t".join(VCF_COLUMN + barcodes) + "\n")
    
    POS_CNT = 0
    vcf_lines_all = []
    for pileupcolumn in samFile.pileup(contig=chrom):
        POS_CNT += 1
        if verbose and POS_CNT % 1000000 == 0:
            print("%s: %dM positions processed." %(chrom, POS_CNT/1000000))
        if pileupcolumn.n < min_COUNT:
            continue

        base_list, qual_list, read_list = pileup_bases(pileupcolumn, 
                                                       pileupcolumn.pos)
        RV = filter_reads(read_list, cell_tag, UMI_tag, min_MAPQ, max_FLAG, 
                          min_LEN)
        UMIs_list = RV["UMIs_list"]
        cell_list = RV["cell_list"]
        base_list = [base_list[ii] for ii in RV["idx_list"]]
        qual_list = [qual_list[ii] for ii in RV["idx_list"]]
        read_list = [read_list[ii] for ii in RV["idx_list"]]
        
        if len(cell_list) < min_COUNT:
            continue
        base_merge, base_cells, qual_cells = map_barcodes(base_list, qual_list, 
            cell_list, UMIs_list, barcodes)
        
        vcf_line = get_vcf_line(base_merge, base_cells, qual_cells, 
            pileupcolumn.reference_name, pileupcolumn.pos, min_COUNT, min_MAF,
            doublet_GL = doublet_GL)

        if vcf_line is not None:
            if out_file is None:
                vcf_lines_all.append(vcf_line)
            else:
                fid.writelines(vcf_line)
    
    if out_file is not None:
        fid.close() 
    return vcf_lines_all
