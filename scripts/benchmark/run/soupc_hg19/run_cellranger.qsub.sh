#!/bin/bash
## declare a name for this job
#PBS -N cr_count
# request the queue for this job
#PBS -q cgsd
# request a total of x processors for this job (y nodes and z processors per node)
#PBS -l nodes=1:ppn=31
# request memory
#PBS -l mem=210gb
# specify walltime
#PBS -l walltime=200:00:00
# out log file
#PBS -o run.out.log
# err log file
#PBS -e run.err.log

#this script is aimed to run cellranger count pipeline on 64 souporcell cell line fastq data
#hxj5<hxj5@hku.hk>

bin_cr=~/tools/cellranger-4.0.0/cellranger
work_dir=/home/xianjie/data/cellsnp-bm/bm_sc
fq_dir=$work_dir/fq
out_dir=$work_dir/cr_hg19
sample_file=$work_dir/sub.samples.lst

sample_lst=`cat $sample_file | tr '\n' ',' | sed 's/,$//'`
id=SAMEA4810598_euts_1
fa_dir=~/data/cellranger/refdata-cellranger-hg19-3.0.0
chem=SC3Pv2
ncores=30
mem=200

echo "=> START @`date '+%Y-%m-%d %H:%M:%S'`"
echo "cellranger version: `$bin_cr --version`"
cd $out_dir
cmd="$bin_cr count --id=$id 
                   --transcriptome=$fa_dir 
                   --fastqs=$fq_dir 
                   --chemistry=$chem 
                   --sample=$sample_lst 
                   --localcores=$ncores 
                   --localmem=$mem"
echo "$cmd" && eval $cmd
echo "=> END @`date '+%Y-%m-%d %H:%M:%S'`"
