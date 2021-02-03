#!/bin/bash
# declare a name for this job
#PBS -N soupc-pileup-p32
# request the queue for this job
#PBS -q cgsd
# request a total of x processors for this job (y nodes and z processors per node)
#PBS -l nodes=1:ppn=32
# request memory
#PBS -l mem=200gb
# specify walltime
#PBS -l walltime=100:00:00
# out log file
#PBS -o soupc-pileup-p32.out
# err log file
#PBS -e soupc-pileup-p32.err

### Pileup scRNA-seq bam file, using cellSNP, cellsnp-lite and vartrix separately,
### with pre-compiled input vcf
vcf=~/data/cellsnp-bm/bm1_demux/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.noindel.nobiallele.nodup.vcf.gz

ncores=32
bam=~/data/cellsnp-bm/bm_sc/cr_hg19/SAMEA4810598_euts_1/outs/possorted_genome_bam.bam
barcode=~/data/cellsnp-bm/bm_sc/cr_hg19/SAMEA4810598_euts_1/outs/filtered_feature_bc_matrix/barcodes.tsv
fa=~/data/cellranger/refdata-cellranger-hg19-3.0.0/fasta/genome.fa
out_dir=~/projects/csp-bm/result/soupc_020321/pileup_p$ncores

# path to cellSNP (v0.3.2)
bin_cellsnp=~/.anaconda3/envs/CSP/bin/cellSNP

# path to cellsnp-lite (v1.2.0, packaged by conda)
bin_cellsnp_lite=~/.anaconda3/envs/CSP/bin/cellsnp-lite

# path to vartrix (v1.1.16, downloaded from Github
# https://github.com/10XGenomics/vartrix/releases/download/v1.1.16/vartrix_linux)
bin_vartrix=~/bin/vartrix

# run cellSNP
/usr/bin/time --verbose $bin_cellsnp -s $bam -O $out_dir/cellSNP \
  -R $vcf -b $barcode -p $ncores --minMAF 0 --minLEN 0 --minMAPQ 20 \
  --minCOUNT 1 --cellTAG CB --UMItag UB --maxFLAG 4096 \
  > $out_dir/cellSNP.out 2> $out_dir/cellSNP.err

# run cellsnp-lite
/usr/bin/time --verbose $bin_cellsnp_lite -s $bam -O $out_dir/cellsnp-lite \
  -R $vcf -b $barcode -p $ncores --minMAF 0 --minLEN 0 --minMAPQ 20 \
  --minCOUNT 1 --gzip --cellTAG CB --UMItag UB --exclFLAG 772 \
  > $out_dir/cellsnp-lite.out 2> $out_dir/cellsnp-lite.err

# run vartrix
/usr/bin/time --verbose $bin_vartrix --primary-alignments --umi \
  --mapq 20 -b $bam --bam-tag CB -c $barcode --scoring-method coverage \
  --threads $ncores --ref-matrix $out_dir/vartrix.ref.mtx \
  --out-matrix $out_dir/vartrix.alt.mtx --out-variants $out_dir/vartrix.variants.tsv \
  -v $vcf --fasta $fa > $out_dir/vartrix.out 2> $out_dir/vartrix.err

