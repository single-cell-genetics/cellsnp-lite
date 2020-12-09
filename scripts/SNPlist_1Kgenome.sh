#!/bin/sh

DAT_DIR=/hps/nobackup/stegle/users/huangh/donorID/genotypes/genome1K

### download data
mkdir $DAT_DIR/chroms
cd $DAT_DIR/chroms
for chr in `seq 1 22`
do
    echo $chr
    file=ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$file
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$file.tbi
done
file=ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$file.tbi


### only keep SNP with AF>0; remove samples
for chr in `seq 1 22`
do
    echo $chr
    VCF_IN=$DAT_DIR/chroms/ALL.chr$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
    OUT_FILE=$DAT_DIR/chroms/ALL_SNP.chr$chr.vcf.gz
    bcftools view -i 'AF>0.0 & TYPE="snp"' -s "." $VCF_IN -O z -o $OUT_FILE --force-samples #&
done
VCF_IN=$DAT_DIR/chroms/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz
OUT_FILE=$DAT_DIR/chroms/ALL_SNP.chrX.vcf.gz
bcftools view -i 'AF>0.0 & TYPE="snp"' -s "." $VCF_IN -O z -o $OUT_FILE --force-samples #&


### concat chr1 to chr22 and chrX
OUT_FILE=$DAT_DIR/genome1K.phase3.SNP_All.chr1toX.vcf.gz
chrom_file=$DAT_DIR/chroms/ALL_SNP.chr1.vcf.gz
for chr in `seq 2 22` X
do 
    echo $chr    
    chrom_file="$chrom_file $DAT_DIR/chroms/ALL_SNP.chr$chr.vcf.gz"
done
echo $chrom_file
bcftools concat $chrom_file -Oz -o $OUT_FILE


### filter SNP by allele frequence and lite INFO
INPUT=$DAT_DIR/genome1K.phase3.SNP_All.chr1toX.vcf.gz
OUT_FILE=$DAT_DIR/genome1K.phase3.SNP_AF5e4.chr1toX.vcf.gz
bcftools annotate -x ^INFO/AF -i 'AF>0.0005' $INPUT -O z -o $OUT_FILE


### liftOver to hg38
## Read here: https://github.com/huangyh09/cellSNP/tree/master/liftOver
## change $CHAIN_PATH and $PATH_cellSNP to the real path
cd $DAT_DIR
CHAIN=$CHAIN_PATH/hg19ToHg38.over.chain.gz
BIN_DIR=$PATH_cellSNP/cellSNP/liftOver
python $BIN_DIR/liftOver_vcf.py -c $CHAIN -i $DAT_DIR/genome1K.phase3.SNP_AF5e4.chr1toX.vcf.gz \
    -o $DAT_DIR/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz

