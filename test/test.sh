#!/bin/sh

DOWNLOAD=${1:-FALSE}
DAT_DIR=$HOME/links/test_cellSNP
cd $DAT_DIR

### Download data (bam file: 19G)
if [ $DOWNLOAD = "TRUE" ]; then
    wget https://sra-pub-src-1.s3.amazonaws.com/SRR5398236/B.merged.bam.1 \ 
        -O $DAT_DIR/demux.B.merged.bam
    wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560246/suppl/GSM2560246%5Fbarcodes%2Etsv%2Egz \
        -O $DAT_DIR/demux.B.barcodes.tsv.gz
    wget http://ufpr.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.subset.hg19.vcf.gz \
        -O $DAT_DIR/genome1K.subset.hg19.vcf.gz
    
    gzip -d $DAT_DIR/demux.B.barcodes.tsv.gz
    head -400 $DAT_DIR/demux.B.barcodes.tsv > $DAT_DIR/demux.B.barcodes.400.tsv
    samtools index $DAT_DIR/demux.B.merged.bam
fi

BAM=$DAT_DIR/demux.B.merged.bam
BARCODE=$DAT_DIR/demux.B.barcodes.400.tsv
REGION=$DAT_DIR/genome1K.subset.hg19.vcf.gz

### Mode 1: 10x data with SNP list
# OUT_DIR=$DAT_DIR/demux_B
# cellSNP -s $BAM -O $OUT_DIR -R $REGION -p 20 -b $BARCODE #--UMItag None 

### Mode 2: 10x data with SNP list
OUT_DIR=$DAT_DIR/demux_B_chrM
cellSNP -s $BAM -O $OUT_DIR --chrom chrMT -p 1 --UMItag None --cellTAG None # -b $BARCODE #--UMItag None 

