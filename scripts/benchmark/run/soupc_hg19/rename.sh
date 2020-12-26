#!/bin/bash
#this script is aimed to rename fastq files from souporcell article to be the names that meet cell-ranger's requirements.
#hxj5<hxj5@hku.hk>

fq_dir=/home/xianjie/data/cellsnp-bm/bm_sc/fq
for fq in `ls $fq_dir/*.fastq.gz`; do
    fn=`basename $fq`
    sample=${fn%%_*}
    tmp_reads=${fn##*_}
    reads=${tmp_reads%%.fastq.gz}
    new_fq=$fq_dir/${sample}_S1_L001_R${reads}_001.fastq.gz
    cmd="mv $fq $new_fq"
    echo "$cmd" && eval $cmd
done
