# Utilility functions for pooling vcf files
# Author: Yuanhua Huang
# Date: 22/08/2018


import os
import sys
import subprocess

def merge_vcf(out_file, chrom_all):
    """Merge vcf for all chromsomes
    """
    if out_file.endswith(".gz"):
        out_file_use = out_file.split(".gz")[0]
    else:
        out_file_use = out_file
        
    CNT = 0
    fid_out = open(out_file_use, "w")
    for _chrom in chrom_all:
        _file = out_file + ".temp_%s_" %(_chrom)
        with open(_file, "r") as fid_in:
            for line in fid_in:
                if line.startswith("#") and _chrom != chrom_all[0]:
                    continue
                else:
                    CNT += 1
                    fid_out.writelines(line)
        os.remove(_file)
    fid_out.close()
    print("%d lines in final vcf file" %CNT)
    
    bashCommand = "gzip -f %s" %(out_file_use)
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]
            
    return None
    