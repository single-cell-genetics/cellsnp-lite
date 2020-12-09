#!/usr/bin/env python
#-*-coding:utf-8-*-
#this script is aimed to fix ref & alt for cellsnp mode 2 vcf and output new Ref & Alt matrices
#hxj5<hxj5@hku.hk>

# https://github.com/hxj5/cellSNP-project/blob/bm3_yy_110520a/scripts/bm/bcf2mtx.py

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

def csp2mtx(vcff0, vcff1, out_dir):
    """
    @param vcff0    Sorted vcf containing correct ref & alt in target region [str]
    @param vcff1    Sorted cellsnp mode 2 vcf (in target region) to be fixed [str]
    @param out_dir  Path to output dir for matrices [str]
    @return         0 if success, -1 otherwise [int]
    """
    try:
        vf0 = pysam.VariantFile(vcff0)
        vf1 = pysam.VariantFile(vcff1)
    except:
        return -1

    ref_mtx = os.path.join(out_dir, "cellsnp.ref.mtx")
    alt_mtx = os.path.join(out_dir, "cellsnp.alt.mtx")
    oth_mtx = os.path.join(out_dir, "cellsnp.oth.mtx")
    var_file = os.path.join(out_dir, "cellsnp.variants.tsv")
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
    for _, a in snps0.items():
        nsnp += len(a)

    ret = 0
    nr_ref = nr_alt = nr_oth = nr_var = 0
    ns_total = ns_fwd = ns_rev = ns_oth = 0
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

        # find the query snp
        n = len(snps0.get(chrom, []))
        assert n > 0, "Error: %s:%d is not in target region." % (chrom, pos) 
        while i < n and pos > snps0[chrom][i][0]:
            i += 1
        assert i < n and pos == snps0[chrom][i][0], "Error: %s:%d is not in target region." % (chrom, pos)

        # query snp found
        nsample = len(rec.samples.keys())   # CHECKME!!! each rec has the same number of samples?
        ref, alt = snps0[chrom][i][1:3]
        assert ref in ("A", "C", "G", "T", "N"), "Error: No.%d SNP %s:%d in target-region vcf has no valid ref" % (i + 1, chrom, pos)
        ref_idx = ("A", "C", "G", "T", "N").index(ref)
        assert alt in ("A", "C", "G", "T", "N"), "Error: No.%d SNP %s:%d in target-region vcf has no valid alt" % (i + 1, chrom, pos)
        alt_idx = ("A", "C", "G", "T", "N").index(alt)

        ns_total += 1
        if rec.ref == ref and rec.alts and rec.alts[0] == alt: ns_fwd += 1
        elif rec.alts and rec.alts[0] == ref and rec.ref == alt: ns_rev += 1
        else: ns_oth += 1

        # iter each sample
        for smp_i, smp in enumerate(rec.samples.values()):
            smp_all = smp.get("ALL", None)
            assert smp_all, "Error: %s:%d in Sample %d has no field ALL" % (chrom, pos, smp_i + 1)
            if len(smp_all) < 5:
                continue
            vref = valt = voth = 0
            vtotal = sum(smp_all)
            if vtotal <= 0:
                continue
            vref = smp_all[ref_idx]
            if vref > 0:
                ref_fp.write("%d\t%d\t%d\n" % (nr_var + 1, smp_i + 1, vref))
                nr_ref += 1
                write_any = 1
            valt = smp_all[alt_idx]
            if valt > 0:
                alt_fp.write("%d\t%d\t%d\n" % (nr_var + 1, smp_i + 1, valt))
                nr_alt += 1
                write_any = 1
            voth = vtotal - vref - valt
            if voth > 0:
                oth_fp.write("%d\t%d\t%d\n" % (nr_var + 1, smp_i + 1, voth))
                nr_oth += 1
                write_any = 1
        if write_any:
            var_fp.write("\t".join([chrom, str(pos), ref, alt]) + "\n")
            nr_var += 1

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
           "  --vcf0 FILE    Sorted vcf containing correct ref & alt in target region\n"           \
           "  --vcf1 FILE    Sorted cellsnp mode 2 vcf (in target region) to be fixed\n"          \
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
    ret = csp2mtx(vcf0, vcf1, out_dir)
    if ret < 0:
        sys.stderr.write("Error: failed to run csp2mtx.\n")
        sys.exit(1)
