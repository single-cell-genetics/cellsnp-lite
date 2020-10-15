Full parameters
---------------
Here is a list of full parameters for setting (``cellsnp-lite -V`` always give the 
version you are using):

.. code-block:: html

  Usage: cellsnp-lite [options]
  
  Options:
    -s, --samFile STR    Indexed sam/bam file(s), comma separated multiple samples.
                         Mode 1&2: one sam/bam file with single cell.
                         Mode 3: one or multiple bulk sam/bam files,
                         no barcodes needed, but sample ids and regionsVCF.
    -S, --samFileList FILE   A list file containing bam files, each per line, for Mode 3.
    -O, --outDir DIR         Output directory for VCF and sparse matrices.
    -R, --regionsVCF FILE    A vcf file listing all candidate SNPs, for fetch each variants.
                             If None, pileup the genome. Needed for bulk samples.
    -b, --barcodeFile FILE   A plain file listing all effective cell barcode.
    -i, --sampleList FILE    A list file containing sample IDs, each per line.
    -I, --sampleIDs STR      Comma separated sample ids.
    -V, --version            Print software version and exit.
    -h, --help               Show this help message and exit.
  
  Optional arguments:
    --genotype           If use, do genotyping in addition to counting.
    --gzip               If use, the output files will be zipped into BGZF format.
    --printSkipSNPs      If use, the SNPs skipped when loading VCF will be printed.
    -p, --nproc INT      Number of subprocesses [1]
    --chrom STR          The chromosomes to use, comma separated [1 to 22]
    --cellTAG STR        Tag for cell barcodes, turn off with None [CB]
    --UMItag STR         Tag for UMI: UR, Auto, None. For Auto mode, use UR if barcodes is inputted,
                         otherwise use None. None mode means no UMI but read counts [Auto]
    --minCOUNT INT       Minimum aggragated count [20]
    --minMAF FLOAT       Minimum minor allele frequency [0.00]
    --doubletGL          If use, keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5.
  
  Read filtering:
    --inclFLAG STR|INT   Required flags: skip reads with all mask bits unset []
    --exclFLAG STR|INT   Filter flags: skip reads with any mask bits set [UNMAP,SECONDARY,QCFAIL
                         (when use UMI) or UNMAP,SECONDARY,QCFAIL,DUP (otherwise)]
    --minLEN INT         Minimum mapped length for read filtering [30]
    --minMAPQ INT        Minimum MAPQ for read filtering [20]
    --countORPHAN        If use, do not skip anomalous read pairs.
  
  Note that the "--maxFLAG" option is now deprecated, please use "--inclFLAG" or "--exclFLAG" instead.
  You can easily aggregate and convert the flag mask bits to an integer by refering to:
  https://broadinstitute.github.io/picard/explain-flags.html
