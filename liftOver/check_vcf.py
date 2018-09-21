##

import pysam

def vcf_check_with_fasta(vcf_file, fasta_file, n_lines=1000):
    faFile = pysam.FastaFile(fasta_file)
    if vcf_file[-3:] == ".gz":
        is_gzip = True
        fid_in = gzip.open(vcf_file, "r")
    else:
        is_gzip = False
        fid_in = open(vcf_file, "r")

    ii = -1
    fasta_base, vcf_REF, vcf_ALT = [], [], []
    for line in fid_in:
        if is_gzip:
            line = line.decode('utf-8')
        if line.startswith("#"):
            continue

        ii += 1
        if ii == n_lines:
            break

        line_val = line.rstrip().split()
        vcf_REF.append(line_val[3])
        vcf_ALT.append(line_val[4])
        _base = faFile.fetch(line_val[0], int(line_val[1])-1, 
                             int(line_val[1]))[0]
        fasta_base.append(_base)
    fid_in.close()
    return fasta_base, vcf_REF, vcf_ALT