#!/bin/bash 
#this script is aimed to compare different softwares for mode 2 with different cores.
#hxj5<hxj5@hku.hk>

root_dir=./portal
if [ ! -d $root_dir ]; then
    echo "Error: $root_dir does not exist." >&2
    exit 1
fi

root_dir=`cd ${root_dir}; pwd`

bin_mc_mode2=$root_dir/run/bm2/multi_core_mode2.sh
out_dir=~/projects/csp-bm/result/submit2/bm2_soupc/run
mkdir -p $out_dir &> /dev/null
multi_cores=8,16,32
tools=cellSNP,cellsnp-lite,freebayes
data_dir=$root_dir/data/bm2_soupc

# run multi core mode 2
$bin_mc_mode2 \
  --bam $data_dir/soupc.bam \
  --fasta $data_dir/cellranger.hg19.3.0.0.fa \
  -t $tools \
  -p $multi_cores \
  -O $out_dir \
  --rootdir $root_dir
