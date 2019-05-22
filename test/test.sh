#!/bin/sh

DOWNLOAD=${1:-FALSE}
DAT_DIR=$HOME/links/test_cellSNP
cd $DAT_DIR

### Download data (bam file: 19G)
if [ $DOWNLOAD = "TRUE" ]; then
    wget https://sra-download.ncbi.nlm.nih.gov/traces/sra47/SRZ/005398/SRR5398236/B.merged.bam \ 
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

OUT_FILE=$DAT_DIR/demux.B.cellGT.v010.vcf.gz
cellSNP -s $BAM -o $OUT_FILE -R $REGION -p 20 -b $BARCODE #--UMItag None 

