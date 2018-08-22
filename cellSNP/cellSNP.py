# pileup SNPs across the genome with pysam
# Author: Yuanhua Huang
# Date: 22-08-2018

import os
import sys
import gzip
import time
import math
import pysam
import subprocess
import multiprocessing
from optparse import OptionParser, OptionGroup
from .utils.pileup_utils import pileup_allele, VCF_HEADER

START_TIME = time.time()

def show_progress(RV=None):
    return RV

def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
        help=("Output directory for VCF files."))
    parser.add_option("--samFile", "-s", dest="sam_file", default=None,
        help=("An indexed sorted sam file."))
    
    group = OptionGroup(parser, "Optional arguments")
    group.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")
    group.add_option("--minCount", type="int", dest="min_count", 
        default=200, help="Minimum aggragated count [default: %default]")
    group.add_option("--minMAF", type="float", dest="min_MAF", 
        default=0.15, help="Minimum minor allele frequency [default: %default]")
    group.add_option("--minMAC", type="float", dest="min_MAC", 
        default=5, help="Minimum minor allele count [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to cellSNP!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    if options.sam_file is None:
        print("Error: need samFile for sam file.")
        sys.exit(1)
    elif os.path.isfile(options.sam_file) == False:
        print("Error: No such file\n    -- %s" %options.sam_file)
        sys.exit(1)
    else:
        sam_file = options.sam_file
        
    if options.out_dir is None:
        print("Error: need outDir for output directory.")
        sys.exit(1)
    elif os.path.isdir(options.out_dir) == False:
        print("Error: No such directory\n    -- %s" %options.out_dir)
        sys.exit(1)
    else:
        out_dir = options.out_dir
    
    nproc = options.nproc
    min_MAF = options.min_MAF
    min_MAC = options.min_MAC
    min_count = options.min_count
    
    chrom_all = [str(x) for x in range(1, 23)]
    # pileup in each chrom
    result = []
    if nproc > 1:
        pool = multiprocessing.Pool(processes=nproc)
        for _chrom in chrom_all:
            out_file = out_dir + "/pileup.%s.vcf" %(_chrom)
            result.append(pool.apply_async(pileup_allele, (sam_file, _chrom, 
                None, None, "CR", min_count, min_MAF,  min_MAC, True, out_file), 
                                           callback=show_progress))
        pool.close()
        pool.join()
    else:
        for _chrom in chrom_all:
            out_file = out_dir + "/pileup.%s.vcf" %(_chrom)
            pileup_allele(sam_file, _chrom, None, None, "CR", min_count, 
                          min_MAF,  min_MAC, True, out_file)
            show_progress(1)
    result = [res.get() if nproc > 1 else res for res in result]
    
    # print("")
    # print("%.1f sec" %(time.time() - START_TIME))
    
    # bashCommand = "gzip -f %s" %(options.out_file.split(".gz")[0])
    # pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    # output = pro.communicate()[0]
        
if __name__ == "__main__":
    main()