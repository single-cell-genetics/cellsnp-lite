# pileup SNPs across the genome with pysam
# Author: Yuanhua Huang
# Date: 22-08-2018

import os
import sys
import gzip
import time
import pysam
import subprocess
import multiprocessing
from optparse import OptionParser, OptionGroup
from .utils.run_utils import merge_vcf
from .utils.pileup_utils import pileup_10X, VCF_HEADER

START_TIME = time.time()

def show_progress(RV=None):
    return RV

def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--outFile", "-o", dest="out_file", default=None,
        help=("Output file path and name for VCF file."))
    parser.add_option("--samFile", "-s", dest="sam_file", default=None,
        help=("An indexed sorted sam file."))
    parser.add_option("--barcodeFile", "-b", dest="barcode_file", default=None,
        help=("A plain file listing all effective cell barcode."))
    
    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")
    group.add_option("--cellTAG", dest="cell_tag", default="CR", 
        help="Tag for cell barcodes, turn off with None [default: %default]")
    group.add_option("--UMItag", dest="UMI_tag", default="UR", 
        help="Tag for UMI, turn off with None [default: %default]")
    group.add_option("--chrom", dest="chrom_all", default=None, 
        help="The chromosomes to use, comma separated [default: 1 to 22]")
    group.add_option("--minCOUNT", type="int", dest="min_COUNT", default=100, 
        help="Minimum aggragated count [default: %default]")
    group.add_option("--minMAF", type="float", dest="min_MAF", default=0.10, 
        help="Minimum minor allele frequency [default: %default]")
    group.add_option("--minMAPQ", type="int", dest="min_MAPQ", default=20, 
        help="Minimum MAPQ for read filtering [default: %default]")
    group.add_option("--maxFLAG", type="int", dest="max_FLAG", default=4096, 
        help="Maximum FLAG for read filtering [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to cellSNP!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    if options.barcode_file is None:
        print("Error: need barcodeFile for cell barcode file.")
        sys.exit(1)
    elif os.path.isfile(options.barcode_file) == False:
        print("Error: No such file\n    -- %s" %options.barcode_file)
        sys.exit(1)
    else:
        fid = open(options.barcode_file, "r")
        barcodes = [x.rstrip().split("-")[0] for x in fid.readlines()]
        fid.close()
        barcodes = sorted(barcodes)
        print("[cellSNP] %d effective cell barcodes are used." %len(barcodes))
        
    if options.sam_file is None:
        print("Error: need samFile for sam file.")
        sys.exit(1)
    elif os.path.isfile(options.sam_file) == False:
        print("Error: No such file\n    -- %s" %options.sam_file)
        sys.exit(1)
    else:
        sam_file = options.sam_file
        
    if options.out_file is None:
        print("Error: need outFile for output file path and name.")
        sys.exit(1)
    elif os.path.isdir(os.path.dirname(options.out_file)) == False:
        print("Error: No such directory for file\n -- %s" %options.out_file)
        sys.exit(1)
    else:
        out_file = options.out_file      
    
    if options.cell_tag.upper() == "NONE":
        cell_tag = None
    else:
        cell_tag = options.cell_tag
    if options.UMI_tag.upper() == "NONE":
        UMI_tag = None
    else:
        UMI_tag = options.UMI_tag
    if options.chrom_all is None:
        chrom_all = [str(x) for x in range(1, 23)]
    else:
        chrom_all = options.chrom_all.split(",")
        
    nproc = options.nproc
    min_MAF = options.min_MAF
    min_MAPQ = options.min_MAPQ
    max_FLAG = options.max_FLAG
    min_COUNT = options.min_COUNT
    
    # pileup in each chrom
    result = []
    if nproc > 1:
        pool = multiprocessing.Pool(processes=nproc)
        for _chrom in chrom_all:
            chr_out_file = out_file + ".temp_%s_" %(_chrom)
            result.append(pool.apply_async(pileup_10X, (sam_file, barcodes, 
                chr_out_file, _chrom, cell_tag, UMI_tag, min_COUNT, min_MAF, 
                min_MAPQ, max_FLAG, True), callback=show_progress))
        pool.close()
        pool.join()
    else:
        for _chrom in chrom_all:
            chr_out_file = out_file + ".temp_%s_" %(_chrom)
            pileup_10X(sam_file, barcodes, chr_out_file, _chrom, cell_tag, 
                       UMI_tag, min_COUNT, min_MAF, min_MAPQ, max_FLAG, True)
            show_progress(1)
    result = [res.get() if nproc > 1 else res for res in result]
    
    print("")
    print("[cellSNP] Whole genome pileupped, now merging all variants ...")
    
    merge_vcf(out_file, chrom_all)
    
    run_time = time.time() - START_TIME
    print("")
    print("[cellSNP] All done: %d min %.1f sec" %(int(run_time / 60), 
                                                  run_time % 60))
    
    # print("")
    # print("%.1f sec" %(time.time() - START_TIME))
        
if __name__ == "__main__":
    main()