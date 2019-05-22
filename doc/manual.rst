Full parameters
---------------
Here is a list of full parameters for setting (``cellSNP -h`` always give the 
version you are using):

.. code-block:: html

  Usage: cellSNP [options]

  Options:
    -h, --help            show this help message and exit
    -s SAM_FILE, --samFile=SAM_FILE
                          Indexed sorted sam file(s), comma separated multiple 
                          samples. Mode 1&2: one sam file with single cell 
                          barcode; Mode 3: one or multiple bulk sam files, no 
                          barcodes needed, but sample ids and regionsVCF.
    -o OUT_FILE, --outFile=OUT_FILE
                          Output file path and name for VCF file.
    -R REGION_FILE, --regionsVCF=REGION_FILE
                          A vcf file listing all candidate SNPs, for fetch each 
                          variants. If None, pileup the genome. Needed for bulk 
                          samples.
    -b BARCODE_FILE, --barcodeFile=BARCODE_FILE
                          A plain file listing all effective cell barcode.
    -I SAMPLE_IDS, --sampleIDs=SAMPLE_IDS
                          Comma separated sample ids. Only use it when you input
                          multiple bulk sam files.

    Optional arguments:
      -p NPROC, --nproc=NPROC
                          Number of subprocesses [default: 1]
      --chrom=CHROM_ALL   The chromosomes to use, comma separated [default: 1 to 22]
      --cellTAG=CELL_TAG  Tag for cell barcodes, turn off with None [default: CB]
      --UMItag=UMI_TAG    Tag for UMI: UR, Auto, None. For Auto mode, use UR if 
                          barcodes is inputted, otherwise use None. None mode 
                          means no UMI but read counts [default: Auto]
      --minCOUNT=MIN_COUNT
                          Minimum aggragated count [default: 20]
      --minMAF=MIN_MAF    Minimum minor allele frequency [default: 0.0]
      --doubletGL         If use, keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5

    Read filtering:
      --minLEN=MIN_LEN    Minimum mapped length for read filtering [default: 30]
      --minMAPQ=MIN_MAPQ  Minimum MAPQ for read filtering [default: 20]
      --maxFLAG=MAX_FLAG  Maximum FLAG for read filtering [default: 255]

