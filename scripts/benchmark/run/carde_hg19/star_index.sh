#!/bin/bash

genome_dir=/home/xianjie/data/cellsnp-bm/bm3_fb/star
fasta_file=/home/xianjie/data/cellranger/refdata-cellranger-hg19-3.0.0/fasta/genome.fa
/home/xianjie/.anaconda3/envs/CSP/bin/STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $genome_dir --genomeFastaFiles $fasta_file
