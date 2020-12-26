#!/bin/bash
#this script is aimed to analysis results of cellsnp-lite mode 2 and bcftools (no region file as inputs)
#hxj5<hxj5@hku.hk>

root_dir=./portal
if [ ! -d "$root_dir" ]; then
    echo "Error: valid root dir needed!" >&2
    exit 1
fi

root_dir=`cd ${root_dir}; pwd`

bin_analysis=$root_dir/run/bm2_bm3/chrom/chrom_analysis.sh
data_dir=~/projects/csp-bm/result/submit2/bm2_bm3/
cellsnp_dir=$data_dir/run/result/cellsnp-lite_ncores8_rep1
bcftools_dir=$data_dir/run/result/bcftools_ncores8_rep1
perf_dir=$data_dir/run/perf
out_dir=$data_dir/chrom
mkdir -p $out_dir &> /dev/null

qsub -q cgsd -N bm2_bm3_chrom -l nodes=1:ppn=10,mem=60gb,walltime=100:00:00 -o $out_dir/run.out.log \
     -e $out_dir/run.err.log -- $bin_analysis --cellsnp-dir $cellsnp_dir --bcftools-dir $bcftools_dir \
     --perf $perf_dir --ncores 10 --out-dir $out_dir --rootdir $root_dir
