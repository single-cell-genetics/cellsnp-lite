#!/bin/bash
#this script is aimed to analysis results of different apps for mode 2
#hxj5<hxj5@hku.hk>

root_dir=./portal
if [ ! -d "$root_dir" ]; then
    echo "Error: valid root dir needed!" >&2
    exit 1
fi

root_dir=`cd ${root_dir}; pwd`

bin_analysis=$root_dir/run/bm2/pr/mode2_analysis.sh
data_dir=~/projects/csp-bm/result/submit2/bm2_soupc
cellsnp_gt=$data_dir/run/result/cellsnp-lite_ncores8_rep1/cellSNP.cells.vcf.gz
freebayes_gt=$data_dir/run/result/freebayes_ncores8_rep1/freebayes.vcf
array_gt=~/data/cellsnp-bm/bm_sc/ERZ368756/HPSI0914i-euts_1.wec.gtarray.HumanCoreExome-12_v1_0.20160912.genotypes.vcf.gz
perf_dir=$data_dir/run/perf

out_dir=$data_dir/analysis_pr
mkdir -p $out_dir &> /dev/null
qsub -q cgsd -N bm2_pr -l nodes=1:ppn=10,mem=60gb,walltime=100:00:00 -o $out_dir/analysis.out.log \
     -e $out_dir/analysis.err.log -- $bin_analysis --cellsnp $cellsnp_gt --freebayes $freebayes_gt \
     --array $array_gt --perf $perf_dir --out-dir $out_dir --rootdir $root_dir

#out_dir=$data_dir/analysis_GP0.9
#mkdir -p $out_dir &> /dev/null
#qsub -q cgsd -N bm2_GP0.9 -l nodes=1:ppn=10,mem=60gb,walltime=100:00:00 -o $out_dir/run.out.log \
#     -e $out_dir/run.err.log -- $bin_analysis --cellsnp $cellsnp_gt --freebayes $freebayes_gt \
#     --array $array_gt --perf $perf_dir --out-dir $out_dir --rootdir $root_dir --minGP 0.9

#out_dir=$data_dir/analysis_GP0.99
#mkdir -p $out_dir &> /dev/null
#qsub -q cgsd -N bm2_GP0.99 -l nodes=1:ppn=10,mem=60gb,walltime=100:00:00 -o $out_dir/run.out.log \
#     -e $out_dir/run.err.log -- $bin_analysis --cellsnp $cellsnp_gt --freebayes $freebayes_gt \
#     --array $array_gt --perf $perf_dir --out-dir $out_dir --rootdir $root_dir --minGP 0.99

