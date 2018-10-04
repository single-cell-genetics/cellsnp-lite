# Utilility functions for pileup SNPs
# Author: Yuanhua Huang
# Date: 22/08/2018

## TODO: samFile.fetch is more efficient, but may gives
## low quality reads, e.g., deletion or refskip

## Note, pileup is not the fastest way, fetch reads and deal 
## with CIGARs will be faster.

import sys
import pysam
import cellSNP
from .base_utils import id_mapping, unique_list

VCF_HEADER = (
    '##fileformat=VCFv4.2\n'
    '##source=cellSNP_v%s\n'
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FILTER=<ID=.,Description="Filter info not available">\n'
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="total counts for ALT and '
    'REF">\n'
    '##INFO=<ID=AD,Number=1,Type=Integer,Description="total counts for ALT">\n'
    '##INFO=<ID=OTH,Number=1,Type=Integer,Description="total counts for other '
    'bases from REF and ALT">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="total counts for ALT and '
    'REF">\n'
    '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="total counts for ALT">\n'
    '##FORMAT=<ID=OTH,Number=1,Type=Integer,Description="total counts for other '
    'bases from REF and ALT">\n'
    '##FORMAT=<ID=ALL,Number=5,Type=Integer,Description="total counts for all '
    'bases in order of A,C,G,T,N">\n' %cellSNP.__version__)

CONTIG = "".join(['##contig=<ID=%s>\n' %x for x in list(range(1,23))+['X', 'Y']])
header_line="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

VCF_COLUMN = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", 
              "INFO", "FORMAT"]

BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
BASE_ZERO = {"A": 0, "C": 0, "G": 0, "T": 0, "N": 0}


def check_pysam_chrom(samFile, chrom=None):
    """Chech if samFile is a file name or pysam object, and if chrom format. 
    """
    if type(samFile) == str:
        samFile = pysam.Samfile(samFile)
    if chrom != None:
        if chrom not in samFile.references:
            if chrom.startswith("chr"):
                chrom = chrom.split("chr")[1]
            else:
                chrom = "chr" + chrom
        if chrom not in samFile.references:
            print("Can't find references %s in samFile" %chrom)
            return None, None
    return samFile, chrom


def fetch_bases(samFile, chrom, POS):
    """ Fetch all reads mapped to the genome position.
    """
    base_list, read_list = [], []
    for _read in samFile.fetch(chrom, POS-1, POS):
        try:
            idx = _read.positions.index(POS-1)
        except:
            continue
        _base = _read.query_alignment_sequence[idx].upper()
        base_list.append(_base)
        read_list.append(_read)
    return base_list, read_list


def pileup_bases(pileupColumn):
    """ pileup all reads mapped to the genome position.
    """
    base_list, read_list = [], []
    for pileupread in pileupColumn.pileups:
        # query position is None if is_del or is_refskip is set.
        if pileupread.is_del or pileupread.is_refskip:
            continue
        query_POS = pileupread.query_position
        _read = pileupread.alignment
        _base = _read.query_sequence[query_POS - 1].upper()
        base_list.append(_base)
        read_list.append(_read)
    return base_list, read_list


def filter_reads(read_list, cell_tag="CR", UMI_tag="UR", min_MAPQ=20, 
                 max_FLAG=255, min_LEN=30):
    """Filter reads and check read tag, e.g., cell and UMI barcodes.
    """
    idx_keep, UMIs_list, cell_list = [], [], []
    for i in range(len(read_list)):
        _read = read_list[i]
        if (_read.mapq < min_MAPQ or _read.flag > max_FLAG or 
            len(_read.positions) < min_LEN): 
            continue
        if cell_tag is not None and _read.has_tag(cell_tag) == False: 
            continue
        if UMI_tag is not None and _read.has_tag(UMI_tag) == False: 
            continue
        if UMI_tag is not None:
            UMIs_list.append(_read.get_tag(UMI_tag))
        if cell_tag is not None:
            cell_list.append(_read.get_tag(cell_tag))
        idx_keep.append(i)
    RV = {}
    RV["idx_list"] = idx_keep
    RV["UMIs_list"] = UMIs_list
    RV["cell_list"] = cell_list
    return RV


def fetch_positions(samFile_list, chroms, positions, REF=None, ALT=None, 
                    barcodes=None, sample_ids=None, out_file=None, 
                    cell_tag="CR", UMI_tag="UR", min_COUNT=20, min_MAF=0.1, 
                    min_MAPQ=20, max_FLAG=255, min_LEN=30, verbose=True):
    """Fetch allelic expression for a list of variants across multiple samples.
    Option 1: one single-cell sam file, a list of barcodes
    Option 2: multiple bulk sam files, multiple sample ids
    No support for multiple sam files and barcodes.
    """
    samFile_list = [check_pysam_chrom(x, chroms[0])[0] for x in samFile_list]
    if out_file is not None:
        fid = open(out_file, "w")
        fid.writelines(VCF_HEADER + CONTIG)
        if barcodes is not None:
            fid.writelines("\t".join(VCF_COLUMN + barcodes) + "\n")
        else:
            fid.writelines("\t".join(VCF_COLUMN + sample_ids) + "\n")
    
    POS_CNT = 0
    vcf_lines_all = []
    for i in range(len(positions)):
        POS_CNT += 1
        if verbose and POS_CNT % 10000 == 0:
            print("%.2fM positions processed." %(POS_CNT/1000000))
        
        base_cells_sample = []
        base_merge_sample = BASE_ZERO.copy()
        for samFile in samFile_list:
            samFile, chrom = check_pysam_chrom(samFile, chroms[i])
            base_list, read_list = fetch_bases(samFile, chrom, positions[i])
            RV = filter_reads(read_list, cell_tag, UMI_tag, min_MAPQ, max_FLAG, 
                              min_LEN)
            UMIs_list = RV["UMIs_list"]
            cell_list = RV["cell_list"]
            base_list = [base_list[ii] for ii in RV["idx_list"]]
            read_list = [read_list[ii] for ii in RV["idx_list"]]
            base_merge, base_cells = map_barcodes(base_list, cell_list, 
                                                  UMIs_list, barcodes)
            
            if barcodes is None:
                for _key in base_merge_sample.keys():
                    base_merge_sample[_key] += base_merge[_key]
                base_cells_sample.append(base_cells[0])
        
        if barcodes is None:
            base_merge, base_cells = base_merge_sample, base_cells_sample
            
        if sum(base_merge.values()) < min_COUNT:
                continue  
        
        if REF is not None and ALT is not None:
            _REF, _ALT = REF[i], ALT[i]
            #only support single nucleotide variants
            if len(_REF) > 1 or len(_ALT) > 1:
                continue
        else:
            _REF, _ALT = None, None
        vcf_line = get_vcf_line(base_merge, base_cells, 
            chrom, positions[i], min_COUNT, min_MAF, _REF, _ALT)

        if vcf_line is not None:
            if out_file is None:
                vcf_lines_all.append(vcf_line)
            else:
                fid.writelines(vcf_line)
    
    if out_file is not None:
        fid.close() 
    return vcf_lines_all


def pileup_regions(samFile, barcodes, out_file=None, chrom=None, cell_tag="CR", 
                   UMI_tag="UR", min_COUNT=20, min_MAF=0.1, min_MAPQ=20, 
                   max_FLAG=255, min_LEN=30, verbose=True):
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

        base_list, read_list = pileup_bases(pileupcolumn)
        RV = filter_reads(read_list, cell_tag, UMI_tag, min_MAPQ, max_FLAG, 
                          min_LEN)
        UMIs_list = RV["UMIs_list"]
        cell_list = RV["cell_list"]
        base_list = [base_list[ii] for ii in RV["idx_list"]]
        read_list = [read_list[ii] for ii in RV["idx_list"]]
        
        if len(cell_list) < min_COUNT:
            continue
        base_merge, base_cells = map_barcodes(base_list, cell_list,
            UMIs_list, barcodes)
        
        vcf_line = get_vcf_line(base_merge, base_cells, 
            pileupcolumn.reference_name, pileupcolumn.pos, min_COUNT, min_MAF)

        if vcf_line is not None:
            if out_file is None:
                vcf_lines_all.append(vcf_line)
            else:
                fid.writelines(vcf_line)
    
    if out_file is not None:
        fid.close() 
    return vcf_lines_all


def map_barcodes(base_list, cell_list, UMIs_list, barcodes):
    """map cell barcodes and pileup bases
    """
    base_merge = BASE_ZERO.copy()
    
    if len(base_list) == 0:
        base_cells = [[0,0,0,0,0]] # need check
        return base_merge, base_cells
    
    # count UMI rather than reads
    if len(UMIs_list) == len(base_list):
        UMIs_uniq, UMIs_idx, tmp = unique_list(UMIs_list)
        cell_list = [cell_list[i] for i in UMIs_idx]
        base_list = [base_list[i] for i in UMIs_idx]
        
    if barcodes is not None and len(cell_list) > 0:
        base_cells = [[0,0,0,0,0] for x in barcodes]
        match_idx = id_mapping(cell_list, barcodes, uniq_ref_only=False, 
                               IDs2_sorted=True)

        for i in range(len(base_list)):
            _idx = match_idx[i]
            _base = base_list[i]
            if _idx is not None:
                base_cells[_idx][BASE_IDX[_base]] += 1
                base_merge[_base] += 1
    else:
        for i in range(len(base_list)):
            base_merge[base_list[i]] += 1
        base_cells = [[base_merge[x] for x in "ACGTN"]]
    return base_merge, base_cells


def get_vcf_line(base_merge, base_cells, chrom, POS, min_COUNT, min_MAF, 
                 REF=None, ALT=None):
    """Convert the counts for all bases into a vcf line
    """
    base_sorted = sorted(base_merge, key=base_merge.__getitem__, reverse=True)
    if REF is None or ALT is None:
        REF = base_sorted[0]
        ALT = base_sorted[1]
            
    min_cnt_2nd = min_MAF * sum(base_merge.values())      
    if (sum(base_merge.values()) < min_COUNT or 
        base_merge[base_sorted[1]] < min_cnt_2nd):
        return None

    FORMAT = "AD:DP:OTH:ALL"
    REF_cnt = base_merge[REF]
    ALT_cnt = base_merge[ALT]
    OTH_cnt = sum(base_merge.values()) - REF_cnt - ALT_cnt
    
    INFO = "AD=%d;DP=%d;OTH=%d" %(ALT_cnt, ALT_cnt+REF_cnt, OTH_cnt)
    
    cells_str = []
    for _base_cell in base_cells:
        if sum(_base_cell) == 0:
            cells_str.append(".:.:.:.")
        else:
            _REF_cnt = _base_cell[BASE_IDX[REF]]
            _ALT_cnt = _base_cell[BASE_IDX[ALT]]
            _OTH_cnt = sum(_base_cell) - _REF_cnt - _ALT_cnt
    
            all_str = ",".join([str(x) for x in _base_cell])
            cnt_lst = [str(_ALT_cnt), str(_ALT_cnt + _REF_cnt), str(_OTH_cnt)]
            out_lst = ":".join(cnt_lst + [all_str])
            cells_str.append(out_lst)
    
    vcf_val = [chrom, str(POS), ".", REF, ALT, ".", "PASS", INFO, FORMAT]
    vcf_line = "\t".join(vcf_val + cells_str) + "\n"
    
    return vcf_line
