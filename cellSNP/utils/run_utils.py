# Utilility functions for pooling vcf files
# Author: Yuanhua Huang
# Date: 22/08/2018


import os
import sys
import gzip
import subprocess

def merge_vcf(out_file, out_files):
    """Merge vcf for all chromsomes
    """
    if out_file.endswith(".gz"):
        out_file_use = out_file.split(".gz")[0]
    else:
        out_file_use = out_file
        
    CNT = 0
    fid_out = open(out_file_use, "w")
    for _file in out_files:
        with open(_file, "r") as fid_in:
            for line in fid_in:
                if line.startswith("#") and _file != out_files[0]:
                    continue
                else:
                    CNT += 1
                    fid_out.writelines(line)
        os.remove(_file)
    fid_out.close()
    print("[cellSNP] %d lines in final vcf file" %CNT)
    
    bashCommand = "gzip -f %s" %(out_file_use)
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]
            
    return None
    
    
def parse_vcf_file(vcf_file, SNP_only=True):
    """Parse vcf file
    """
    if vcf_file[-3:] == ".gz" or vcf_file[-4:] == ".bgz":
        infile = gzip.open(vcf_file, "rb")
    else:
        infile = open(vcf_file, "r")
        
    chrom_list, pos_list, ids_list = [], [], []
    REF_list, ALT_list = [], []
    contig_lines = []
    for line in infile:
        line = line.decode('utf-8')
        if line.startswith("#"):
            if line.startswith("##contig="): 
                contig_lines.append(line)
            if line.startswith("#CHROM"):
                list_val = line.rstrip().split("\t")[:8]
        else:
            list_val = line.rstrip().split("\t")[:5] #:8
            if SNP_only:
                if len(list_val[3]) > 1 or len(list_val[4]) > 1:
                    continue
            chrom_list.append(list_val[0])
            pos_list.append(int(list_val[1]))
            ids_list.append(list_val[2])
            REF_list.append(list_val[3])
            ALT_list.append(list_val[4])
    RV = {}
    RV["pos"] = pos_list
    RV["ids"] = ids_list
    RV["REF"] = REF_list
    RV["ALT"] = ALT_list
    RV["chrom"] = chrom_list
    RV["contig_lines"] = contig_lines
    return RV
    