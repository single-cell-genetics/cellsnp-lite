#!/usr/bin/env python
#-*-coding:utf-8-*-
#this script is aimed to convert AD values of bcftools mpileup to cellSNP-style ref & alt matrices
#hxj5<hxj5@hku.hk>

import pysam
import sys
import os
import getopt

def rewrite_mtx(tmp, out, nvar, nsmp, nrec):
    """
    @abstract    Add mtx header to the tmp mtx and rewrite to new mtx
    @param tmp   Path to the tmp mtx which doesn't have mtx header [str]
    @param out   Path to the new mtx [str]
    @param nvar  Number of variants [int]
    @param nsmp  Number of samples [int]
    @param nrec  Number of records [int]
    @return      0 if success, -1 otherwise [int]
    """
    buf_size = 1024 * 1024
    mtx_header = "%%MatrixMarket matrix coordinate integer general\n%\n"
    fp = open(out, "w")
    fp_tmp = open(tmp, "r")
    fp.write(mtx_header)
    fp.write("%d\t%d\t%d\n" % (nvar, nsmp, nrec))
    while True:
        dat = fp_tmp.read(buf_size)
        if dat:
            fp.write(dat)
        else:
            break
    fp_tmp.close()
    fp.close()
    os.remove(tmp)
    return 0

def bcf2mtx(vcff0, vcff1, out_dir):
    """
    @abstract       Convert AD values of bcftools mpileup to cellSNP-style ref & alt matrices
    @param vcff0    Sorted raw input vcf of bcftools mpileup [str]
    @param vcff1    Sorted output vcf of bcftools mpileup [str]
    @param out_dir  Path to output dir for converted matrices [str]
    @return         0 if success, -1 otherwise [int]
    """
    try:
        vf0 = pysam.VariantFile(vcff0)
        vf1 = pysam.VariantFile(vcff1)
    except:
        return -1

    ref_mtx = os.path.join(out_dir, "bcftools.ref.mtx")
    alt_mtx = os.path.join(out_dir, "bcftools.alt.mtx")
    oth_mtx = os.path.join(out_dir, "bcftools.oth.mtx")
    var_file = os.path.join(out_dir, "bcftools.variants.tsv")
    ref_mtx_tmp = ref_mtx + ".tmp"
    alt_mtx_tmp = alt_mtx + ".tmp"
    oth_mtx_tmp = oth_mtx + ".tmp"
    
    ref_fp = open(ref_mtx_tmp, "w")
    alt_fp = open(alt_mtx_tmp, "w")
    oth_fp = open(oth_mtx_tmp, "w")
    var_fp = open(var_file, "w")
    var_fp.write("#CHROM\tPOS\tREF\tALT\n")

    snps0 = {}
    for rec in vf0.fetch():
        snps0.setdefault(rec.contig, []).append((rec.pos, rec.ref, rec.alts[0]))
    vf0.close()
    nsnp = 0
    for c, a in snps0.items():
        nsnp += len(a)

    ret = 0
    nr_ref = nr_alt = nr_oth = nr_var = 0
    i = 0        # index of the next avaliable snp in snp0[chrom]
    nsample = 0
    for rec in vf1.fetch():
        write_any = 0
        chrom, pos = rec.contig, rec.pos   # pos is 1-based.
        if i == 0:
            old_chrom = chrom
        elif old_chrom != chrom:   # new chrom
            i = 0
            old_chrom = chrom
        n = len(snps0.get(chrom, []))
        if n <= 0:         # query snp's chrom not in vcf0
            sys.stderr.write("Error: %s:%d is not in original vcf.\n" % (chrom, pos))
            ret = -1
            break
        while i < n and pos > snps0[chrom][i][0]:
            i += 1
        if i >= n or pos < snps0[chrom][i][0]:       # query snp's pos not in vcf0
            sys.stderr.write("Error: %s:%d is not in original vcf.\n" % (chrom, pos))
            ret = -1
            break
        else:      # find the query snp!
            nsample = len(rec.samples.keys())   # CHECKME!!! each rec has the same number of samples?
            ref, alt = snps0[chrom][i][1:3]
            total_dp = rec.info.get("DP", 0)
            if total_dp is None or total_dp <= 0:
                sys.stderr.write("Warning: skip %s:%d, total DP <= 0\n" % (chrom, pos))
                continue
            #allele_idx = [i for i, a in enumerate(rec.alleles) if a in ("A", "C", "G", "T")]
            #if not allele_idx:
            #    sys.stderr.write("Warning: %s:%d has no valid allele." % (chrom, pos))
            #    continue
            #alleles = [rec.alleles[i] for i in allele_idx]
            try:
                ref_idx = rec.alleles.index(ref)
            except ValueError:
                ref_idx = -1
            try:
                alt_idx = rec.alleles.index(alt)
            except ValueError:
                alt_idx = -1
            try:
                ad_idx = rec.format.keys().index("AD")
            except ValueError:
                sys.stderr.write("Error: %s:%d format has no AD.\n" % (chrom, pos))
                ret = -1
                break
            total_smp_dp = 0
            for smp_i, smp in enumerate(rec.samples.values()):
                smp_v = smp.values()
                assert ad_idx < len(smp_v),  \
                       "Error: %s:%d in sample %d has invalid format." % (chrom, pos, smp_i + 1)
                vref = valt = voth = 0
                vtotal = sum(smp_v[ad_idx])
                total_smp_dp += vtotal
                if vtotal <= 0:
                    continue
                if ref_idx >= 0:
                    vref = smp_v[ad_idx][ref_idx]
                    if vref > 0:
                        ref_fp.write("%d\t%d\t%d\n" % (nr_var + 1, smp_i + 1, vref))
                        nr_ref += 1
                        write_any = 1
                if alt_idx >= 0:
                    valt = smp_v[ad_idx][alt_idx]
                    if valt > 0:
                        alt_fp.write("%d\t%d\t%d\n" % (nr_var + 1, smp_i + 1, valt))
                        nr_alt += 1
                        write_any = 1
                oth_allele_idx = [i for i, a in enumerate(rec.alleles) if a in ("A", "C", "G", "T", "N") and a != ref and a != alt]
                if oth_allele_idx:
                    voth = sum([smp_v[ad_idx][i] for i in oth_allele_idx])
                if voth > 0:
                    oth_fp.write("%d\t%d\t%d\n" % (nr_var + 1, smp_i + 1, voth))
                    nr_oth += 1
                    write_any = 1
            if write_any:
                var_fp.write("\t".join([chrom, str(pos), ref, alt]) + "\n")
                nr_var += 1
            assert total_smp_dp <= total_dp, "%s:%d, total sample DP should be no more than total DP" % (chrom, pos)

    vf1.close()     
    ref_fp.close()
    alt_fp.close()
    oth_fp.close()
    var_fp.close()
    
    if ret < 0:
        return ret

    assert nr_var <= nsnp, "Number of variants in vcff1 should be no more than vcff0."
    rewrite_mtx(ref_mtx_tmp, ref_mtx, nr_var, nsample, nr_ref)
    rewrite_mtx(alt_mtx_tmp, alt_mtx, nr_var, nsample, nr_alt)
    rewrite_mtx(oth_mtx_tmp, oth_mtx, nr_var, nsample, nr_oth)

    return 0

def usage(fp = sys.stderr):
    msg = "\n"
    msg += "Usage: %s [options]\n" % sys.argv[0]
    msg += "\n"                                           \
           "Options:\n"                                   \
           "  --vcf0 FILE    Sorted raw input vcf of bcftools mpileup\n"           \
           "  --vcf1 FILE    Sorted output vcf of bcftools mpileup\n"          \
           "  --outdir DIR   Path to output dir for converted matrices\n"   \
           "  -h, --help     Print this message\n"                  \
           "\n"
    fp.write(msg)

if __name__ == "__main__":    
    if len(sys.argv) < 4:
        usage()
        sys.exit(1)
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["vcf0=", "vcf1=", "outdir=", "help"])
    except getopt.GetoptError as e:
        print(str(e))
        usage()
        sys.exit(2)
    vcf0 = vcf1 = out_dir = None
    for opt, value in opts:
        if opt == "--vcf0": vcf0 = value
        elif opt == "--vcf1": vcf1 = value
        elif opt == "--outdir": out_dir = value
        elif opt in ("-h", "--help"): usage(); sys.exit(3)
        else: assert False, "unhandled option"
    ret = bcf2mtx(vcf0, vcf1, out_dir)
    if ret < 0:
        sys.stderr.write("Error: failed to run bcf2mtx.\n")
        sys.exit(1)

