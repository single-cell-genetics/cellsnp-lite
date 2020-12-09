#!/bin/bash
#this script is aimed to compare results of mode 2 and mode 3
#hxj5<hxj5@hku.hk>

root_dir=./portal
if [ ! -d "$root_dir" ]; then
    echo "Error: valid root dir needed!" >&2
    exit 1
fi

root_dir=`cd ${root_dir}; pwd`

bin_analysis=$root_dir/run/bm2_bm3/snp/snp_analysis.sh
dir2=~/projects/csp-bm/result/submit2/bm2_bm3/run/result/cellsnp-lite_ncores8_rep1
dir3=~/projects/csp-bm/result/submit2/bm3_carde/run/result/cellsnp-lite-R_ncores8_rep1
origin_variants=~/data/cellsnp-bm/bm1_demux/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.noindel.nobiallele.nodup.vcf.gz
out_dir=~/projects/csp-bm/result/submit2/bm2_bm3/snp
mkdir -p $out_dir &> /dev/null

qsub -q cgsd -N bm2_bm3_snp -l nodes=1:ppn=10,mem=60gb,walltime=100:00:00 -o $out_dir/run.out.log \
     -e $out_dir/run.err.log -- $bin_analysis --dir2 $dir2 --dir3 $dir3 \
     --variants $origin_variants --ncores 10 --out-dir $out_dir --rootdir $root_dir
