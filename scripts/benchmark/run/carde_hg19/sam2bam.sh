#!/bin/bash

in_dir=/home/xianjie/data/cellsnp-bm/bm3_fb/aln
out_dir=/home/xianjie/data/cellsnp-bm/bm3_fb/sort
mkdir $out_dir &> /dev/null

for f in `ls $in_dir/*.sam`; do
    sample=`basename ${f%%.*}`
    tmp_bam=$out_dir/${sample}.bam
    bam=$out_dir/${sample}.sort.bam
    echo "=> $sample"
    cmd="samtools view -h -b $f > $tmp_bam"
    echo "$cmd" && eval $cmd
    cmd="samtools sort -O BAM $tmp_bam > $bam"
    echo "$cmd" && eval $cmd
    cmd="rm $tmp_bam"
    echo "$cmd" && eval $cmd
    cmd="samtools index $bam"
    echo "$cmd" && eval $cmd
    echo
done
