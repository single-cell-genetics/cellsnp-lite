#!/bin/bash
# declare a name for this job
#PBS -N SRR6860519-genotype
# request the queue for this job
#PBS -q cgsd
# request a total of x processors for this job (y nodes and z processors per node)
#PBS -l nodes=1:ppn=16
# request memory
#PBS -l mem=100gb
# specify walltime
#PBS -l walltime=100:00:00
# out log file
#PBS -o SRR6860519-genotype.out
# err log file
#PBS -e SRR6860519-genotype.err

### Use a scRNA-seq bam as bulk bam for genotyping, with 
### cellsnp-lite and freebayes separately

# path to cellsnp-lite (v1.2.0, packaged by conda)
bin_cellsnp=~/.anaconda3/envs/CSP/bin/cellsnp-lite

# path to freebayes/freebayes-parallel 
# (freebayes v1.3.2-dirty, packaged by conda)
bin_freebayes=~/.anaconda3/envs/CSP/bin/freebayes-parallel

# path to samtools (v1.10, using htslib v1.10.2, packaged by conda)
bin_samtools=~/.anaconda3/envs/CSP/bin/samtools   

out_dir=~/projects/csp-bm/result/srr_020321/genotype_mode2
fa=~/data/cellranger/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa
fai=~/data/cellranger/refdata-cellranger-GRCh38-1.2.0/fasta/genome.fa.fai
ncores=16

# download the bam file of SRR6860519, which is for donor1-rep1 of GSE112013. 
# and then index bam
bam=$out_dir/SRR6860519.bam
wget https://sra-pub-src-1.s3.amazonaws.com/SRR6860519/14515X1.bam.1 -O $bam
$bin_samtools index $bam

# run cellsnp-lite
chroms="`seq 1 22` X Y"
chroms=`echo $chroms | sed 's/ /,/g'`

# with minMAF=0,minCOUNT=20
/usr/bin/time --verbose $bin_cellsnp -s $bam -O $out_dir/cellsnp-lite-maf0 \
  -p $ncores --minMAF 0 --minCOUNT 20 --gzip --genotype --chrom $chroms \
  --cellTAG None --UMItag None --exclFLAG 772 \
  > $out_dir/cellsnp-lite-maf0.out 2> $out_dir/cellsnp-lite-maf0.err

# with minMAF=0,minCOUNT=1,minLEN=0
/usr/bin/time --verbose $bin_cellsnp -s $bam -O $out_dir/cellsnp-lite-maf0-count1 \
  -p $ncores --minMAF 0 --minCOUNT 1 --gzip --genotype --chrom $chroms \
  --cellTAG None --UMItag None --exclFLAG 772 --minLEN 0 \
  > $out_dir/cellsnp-lite-maf0-count1.out 2> $out_dir/cellsnp-lite-maf0-count1.err

# with minMAF=0.1,minCOUNT=20
/usr/bin/time --verbose $bin_cellsnp -s $bam -O $out_dir/cellsnp-lite-maf0.1 \
  -p $ncores --minMAF 0.1 --minCOUNT 20 --gzip --genotype --chrom $chroms \
  --cellTAG None --UMItag None --exclFLAG 772 \
  > $out_dir/cellsnp-lite-maf0.1.out 2> $out_dir/cellsnp-lite-maf0.1.err

# run freebayes

# tools used by freebayes (eg. fasta_generate_regions.py, vcffirstheader, 
# vcfstreamsort, vcfuniq and parallel) should be added into PATH
export PATH=~/.anaconda3/envs/CSP/bin:$PATH

/usr/bin/time --verbose $bin_freebayes <(fasta_generate_regions.py $fai 100000) \
  $ncores -f $fa --genotype-qualities --use-duplicate-reads $bam \
  > $out_dir/freebayes.vcf 2> $out_dir/freebayes.err

