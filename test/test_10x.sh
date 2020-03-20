#!/bin/sh

## Please edit the $DAT_DIR to the directory you want

DOWNLOAD=${1:-FALSE}
DAT_DIR=$HOME/test_cellSNP
mkdir $DAT_DIR
cd $DAT_DIR

### Download data (including bam file: 123 Mb)
if [ $DOWNLOAD = "TRUE" ]; then
    cd -P $DAT_DIR
    wget http://ufpr.dl.sourceforge.net/project/cellsnp/testdata/cellSNP_testdata_10x.zip
    unzip cellSNP_testdata_10x.zip
fi

BAM=$DAT_DIR/demux.B.lite.bam
BARCODE=$DAT_DIR/demux.B.barcodes.400.tsv
REGION=$DAT_DIR/genome1K.subset.hg19.vcf.gz

### Mode 1: 10x data with SNP list
OUT_DIR=$DAT_DIR/demux_B_list
cellSNP -s $BAM -O $OUT_DIR -R $REGION -b $BARCODE --minCOUNT 20
## -p 4 #--UMItag None


### Mode 2: 10x data without SNP list
OUT_DIR=$DAT_DIR/demux_B_chrom
cellSNP -s $BAM -O $OUT_DIR -b $BARCODE --minCOUNT 20 --minMAF 0.1 -p 22
## --chrom chrMT # -p 1 --UMItag None --cellTAG None

