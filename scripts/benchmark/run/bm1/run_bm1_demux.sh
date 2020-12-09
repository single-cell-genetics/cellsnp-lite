#!/bin/bash 
#this script is aimed to compare different softwares for mode 1 with different cores.
#hxj5<hxj5@hku.hk>

root_dir=./portal
if [ ! -d $root_dir ]; then
    echo "Error: $root_dir does not exist." >&2
    exit 1
fi

root_dir=`cd ${root_dir}; pwd`

bin_mc_mode1=$root_dir/run/bm1/multi_core_mode1.sh
out_dir=~/projects/csp-bm/result/submit2/bm1_demux/run
mkdir -p $out_dir &> /dev/null
multi_cores=8,16,32
tools=cellSNP,cellsnp-lite,vartrix
data_dir=$root_dir/data/bm1_demux

# run multi core mode 1
$bin_mc_mode1 \
  --bam $data_dir/demux.B.merged.bam \
  --snp $data_dir/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.noindel.nobiallele.nodup.vcf.gz \
  --barcode $data_dir/demux.B.barcodes.tsv \
  --fasta $data_dir/cellranger.hg19.3.0.0.fa \
  -t $tools \
  -p $multi_cores \
  -O $out_dir \
  --rootdir $root_dir

