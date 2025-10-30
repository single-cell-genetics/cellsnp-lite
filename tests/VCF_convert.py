# convert cellSNP output from vcf.gz to HDF5 format
# Author: Yuanhua Huang
# Date: 09-06-2019

import os
import sys
from optparse import OptionParser, OptionGroup
from cellSNP.utils.vcf_utils import load_VCF, write_VCF_to_hdf5

def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--vcfFILE", "-i", dest="vcf_file", default=None,
        help=("Input vcf file."))
    parser.add_option("--outFile", "-o", dest="out_file", default=None,
        help=("Output file path and name for HDF5 file."))

    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("Welcome to VCF_convert!\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    if options.vcf_file is None:
        print("Error: need vcf file as input.")
        sys.exit(1)
    else:
        if options.vcf_file.endswith(".gz"):
            vcf_file = options.vcf_file.split(".gz")[0]
        else:
            vcf_file = options.vcf_file

    if options.out_file is None:
        out_file = vcf_file + '.h5'
    elif os.path.dirname(options.out_file) == "":
        out_file = "./" + options.out_file
    else:
        out_file = options.out_file
    if os.path.isdir(os.path.dirname(out_file)) == False:
        print("Error: No such directory for file\n -- %s" %out_file)
        sys.exit(1)  

    ## save to hdf5 file
    vcf_dat = load_VCF(vcf_file + ".gz", load_sample=True, sparse=True)
    write_VCF_to_hdf5(vcf_dat, out_file)


if __name__ == "__main__":
    main()
