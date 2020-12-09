#!/bin/bash 
#this script is aimed to compare different softwares for mode 3 with different cores.
#hxj5<hxj5@hku.hk>

root_dir=./portal
if [ ! -d $root_dir ]; then
    echo "Error: $root_dir does not exist." >&2
    exit 1
fi

root_dir=`cd ${root_dir}; pwd`

bin_mc_mode3=$root_dir/run/bm3/multi_core_mode3.sh
out_dir=~/projects/csp-bm/result/submit2/bm3_carde/run
mkdir -p $out_dir &> /dev/null
multi_cores=1,2,4,8
tools=bcftools-R,bcftools-T,cellSNP,cellsnp-lite-R,cellsnp-lite-T
data_dir=$root_dir/data/bm3_carde

# run multi core mode 1
$bin_mc_mode3 \
  --bam $data_dir/bam.lst \
  --snp $data_dir/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.noindel.nobiallele.nodup.vcf.gz \
  --sample $data_dir/sample.lst \
  --fasta $data_dir/cellranger.hg19.3.0.0.fa \
  -t $tools \
  -p $multi_cores \
  -O $out_dir \
  --rootdir $root_dir
