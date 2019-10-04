## A python wrap of UCSC liftOver function for vcf file
## UCSC liftOver binary and hg19 to hg38 chain file:
## https://genome.ucsc.edu/cgi-bin/hgLiftOver
## http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
## http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

import sys
import gzip
import subprocess
from optparse import OptionParser

LIFTOVER_INFO = '##INFO=<ID=OLD,Number=1,Type=Integer,'
LIFTOVER_INFO += 'Description="position before liftover">\n'

def vcf_to_bed(vcf_file, out_file, chr_in=True):
    if vcf_file[-3:] == ".gz":
        is_gzip = True
        fid_in = gzip.open(vcf_file, "r")
    else:
        is_gzip = False
        fid_in = open(vcf_file, "r")
    
    fid_out = open(out_file, "w")
    for line in fid_in:
        if is_gzip:
            line = line.decode('utf-8')
        if line.startswith("#") == False:
            line_val = line.rstrip().split("\t")[:8]
            if chr_in and line_val[0].startswith("chr") == False:
                line_val[0] = "chr" + line_val[0]
            line_val[2] = str(int(line_val[1]) + 1)
            fid_out.writelines("\t".join(line_val[:3]) + "\n")
    fid_in.close()
    fid_out.close()
    return None

def update_vcf(vcf_file, bed_new, bed_unmap, out_file):       
    ## unmapped lines
    unmap_pos = []
    _fid = open(bed_unmap, "r")
    for line in _fid:
        if not line.startswith("#"):
            _pos_id = "_".join(line.rstrip().split("\t")[:2])
            unmap_pos.append(_pos_id)
    _fid.close()
    
    if vcf_file[-3:] == ".gz":
        is_gzip = True
        fid_in = gzip.open(vcf_file, "r")
    else:
        is_gzip = False
        fid_in = open(vcf_file, "r")
        
    cnt1 = 0 
    idx_unmap = 0
    fid_bed = open(bed_new, "r")
    fid_out = open(out_file, "w")    
    for line in fid_in:
        if is_gzip:
            line = line.decode('utf-8')
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                fid_out.writelines(LIFTOVER_INFO)
            fid_out.writelines(line)
        else:
            line_val = line.rstrip().split("\t")
            if idx_unmap < len(unmap_pos):
                _pos_id = "_".join(line_val[:2])
                if line_val[0].startswith("chr") == False:
                    _pos_id = "chr" + _pos_id
                if _pos_id == unmap_pos[idx_unmap]:
                    idx_unmap += 1
                    continue
            cnt1 += 1
            bed_line = fid_bed.readline()
            line_val[7] = "OLD=" + line_val[1] + ";" + line_val[7]
            line_val[1] = bed_line.rstrip().split("\t")[1]
            fid_out.writelines("\t".join(line_val) + "\n")
    print(cnt1, idx_unmap)
    fid_in.close()
    fid_bed.close()
    fid_out.close()
    return None


def main():
    import warnings
    warnings.filterwarnings('error')

    # parse command line options
    parser = OptionParser()
    parser.add_option("--chainFile", "-c", dest="chain_file", default=None,
        help=("Chain file, full path."))
    parser.add_option("--inFile", "-i", dest="in_file", default=None,
        help=("Input vcf file, full path."))
    parser.add_option("--outFile", "-o", dest="out_file", default=None,
        help=("Output VCF file, full path."))
    parser.add_option("--liftOverPath", "-P", dest="liftOver_path", default=None,
        help=("liftOver_path if it is not in PATH variable."))
    
    (options, args) = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        print("liftOver-vcf: a wrap of UCSC liftOver for VCF file.\n")
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    in_file = options.in_file
    bed_file = options.in_file.split(".vcf")[0] + ".bed"
    new_bed_file = options.out_file.split(".vcf")[0] + ".bed"
    unmap_bed_file = options.out_file.split(".vcf")[0] + ".unmap.bed"
    
    ## generate bed file
    print("converting vcf to bed file ... ")
    vcf_to_bed(in_file, bed_file)
    
    ## UCSC liftOver on bed file
    chain_file = options.chain_file
    if options.liftOver_path is None:
        liftOver = "liftOver"
    else:
        # check if path exists
        liftOver = options.liftOver_path
                
    print("liftOver bed file ... ")
    bashCommand = "%s %s %s %s %s" %(liftOver, bed_file, chain_file, 
                                  new_bed_file, unmap_bed_file)
    #print(bashCommand)
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]

    ## update vcf file
    out_file = options.out_file
    if out_file[-3:] == ".gz":
        out_file = out_file[:-3]
    print("updating vcf file ... ")
    update_vcf(in_file, new_bed_file, unmap_bed_file, out_file)
    
    print("gzip vcf file ... ")
    import shutil
    if shutil.which("bgzip") is not None:
        bashCommand = "bgzip -f %s" %(out_file)
    else:
        bashCommand = "gzip -f %s" %(out_file)
    pro = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    pro.communicate()[0]
    return None

if __name__ == "__main__":
    main()
    