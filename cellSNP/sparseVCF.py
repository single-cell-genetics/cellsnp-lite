# pileup SNPs across the genome with pysam's fetch or pileup reads
# Author: Yuanhua Huang
# Date: 08-07-2019

import os
import sys
import time
from optparse import OptionParser
from utils.vcf_utils import VCF_to_sparseMat

START_TIME = time.time()

def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--vcfFile", "-i", dest="vcf_file", default=None,
        help=("The input vcf file to parse."))
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
        help=("Directory for output files [default: $vcfFile/sparseVCF]."))
    parser.add_option("--tags", "-t", dest="out_tags", default="AD,DP,OTH",
        help=("A comma separated tags to keep in sparse matrix " 
              "[default: %default]"))

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to sparseVCF!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    if options.vcf_file is None:
        print("Error: need vcfFile for vcf file.")
        sys.exit(1)
    elif os.path.isfile(options.vcf_file) == False:
        print("Error: No such file\n    -- %s" %options.vcf_file)
        sys.exit(1)
    else:
        vcf_file = options.vcf_file

    if options.out_dir is None:
        print("Warning: no outDir provided, we use $vcfFilePath/sparseVCF.")
        out_dir = os.path.dirname(os.path.abspath(vcf_file)) + "/sparseVCF"
    elif os.path.dirname(options.out_dir) == "":
        out_dir= "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    out_tags = options.out_tags.split(",")
    VCF_to_sparseMat(vcf_file, out_tags, out_dir)
    
    # run_time = time.time() - START_TIME
    # print("[cellSNP] All done: %d min %.1f sec" %(int(run_time / 60), 
    #                                               run_time % 60))
    
        
if __name__ == "__main__":
    main()
