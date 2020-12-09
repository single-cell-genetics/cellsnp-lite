#!/bin/bash
#this script is aimed to analysis results of different apps for mode 3
#hxj5<hxj5@hku.hk>

root_dir=./portal
if [ ! -d "$root_dir" ]; then
    echo "Error: valid root dir needed!" >&2
    exit 1
fi

root_dir=`cd ${root_dir}; pwd`

bin_analysis=$root_dir/run/bm3/mode3_analysis.sh
data_dir=~/projects/csp-bm/result/submit2/bm3_carde
cellsnp_dir=$data_dir/run/result/cellsnp-lite-R_ncores8_rep1
bcftools_dir=$data_dir/run/result/bcftools-R_ncores8_rep1
perf_dir=$data_dir/run/perf
origin_variants=~/data/cellsnp-bm/bm1_demux/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.noindel.nobiallele.nodup.vcf.gz
out_dir=$data_dir/analysis
mkdir -p $out_dir &> /dev/null

qsub -q cgsd -N bm3_carde_a -l nodes=1:ppn=10,mem=60gb,walltime=100:00:00 -o $out_dir/run.out.log \
     -e $out_dir/run.err.log -- $bin_analysis --cellsnp-dir $cellsnp_dir --bcftools-dir $bcftools_dir \
     --variants $origin_variants --perf $perf_dir --ncores 10 --out-dir $out_dir --rootdir $root_dir

# compare cellsnp-lite-R with cellsnp-lite-T by `diff`
# compare bcftools-R with bcftools-T by `diff`
