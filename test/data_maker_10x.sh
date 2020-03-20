#!/bin/sh

## Please edit the $DAT_DIR to the directory you want

DOWNLOAD=${1:-FALSE}
DAT_DIR=$HOME/test_cellSNP
mkdir $DAT_DIR
cd -P $DAT_DIR

### Download data (bam file: 19G)
if [ $DOWNLOAD = "TRUE" ]; then
    cd -P $DAT_DIR
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

### Subsampling the bam file using synth_pool.py from vireo
## https://github.com/single-cell-genetics/vireo/blob/master/simulate/synth_pool.py
python $HOME/vireo/simulate/synth_pool.py -r $REGION -b $BARCODE -s $BAM -o ./lite_bam -d 0 


## rename file
# cp lite_bam/pooled.sorted.bam demux.B.lite.bam 
# cp lite_bam/pooled.sorted.bam.bai demux.B.lite.bam.bai
