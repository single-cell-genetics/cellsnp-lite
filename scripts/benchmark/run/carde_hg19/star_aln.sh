#!/bin/bash

fq_dir=/home/xianjie/data/cellsnp-bm/bm3_fb
genome_dir=/home/xianjie/data/cellsnp-bm/bm3_fb/star
out_dir=/home/xianjie/data/cellsnp-bm/bm3_fb/aln
mkdir $out_dir &> /dev/null
ncores=8
mem=60   # gb

fq_list=`cd ${fq_dir}; ls *.fastq.gz | tr ' ' '\n' | sed 's/_.*//' | sort -u`
for name in $fq_list; do
    fq1=$fq_dir/${name}_1.fastq.gz
    fq2=$fq_dir/${name}_2.fastq.gz
    out_prefix=$out_dir/$name
    out_log=$out_dir/${name}.out.log
    err_log=$out_dir/${name}.err.log
    cmd="qsub -q cgsd -N $name -l nodes=1:ppn=$[ncores + 1],mem=${mem}gb -o $out_log -e $err_log -- /home/xianjie/.anaconda3/envs/CSP/bin/STAR --runThreadN 8 --genomeDir $genome_dir --readFilesCommand zcat --readFilesIn $fq1 $fq2 --outFileNamePrefix ${out_prefix}."
    echo "$cmd" && eval $cmd
done
