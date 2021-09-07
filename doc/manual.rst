Manual
======

Quick usage
-----------

Once installed, check all arguments by type ``cellsnp-lite -h``. 
There are two modes of cellsnp-lite:

Mode 1: pileup with given SNPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This mode genotypes single cells or bulk sample at a list of given SNPs, which 
could be common SNPs in human population (see `compiled candidate SNPs`_), or
called heterouzygous variants from ``Mode 2b`` on its own.

.. _compiled candidate SNPs: snp_list.html


Mode 1a: droplet-based single cells
+++++++++++++++++++++++++++++++++++

Use both ``-R`` and ``-b`` to pileup droplet-based dataset (e.g., 10x Genomics) with given SNPs.

Require: a single BAM/SAM/CRAM file, e.g., from cellranger, a list of cell barcodes,
a VCF file for common SNPs. This mode is recommended comparing to mode 2, if a
list of common SNP is known, e.g., human (see `Candidate_SNPs`_)

.. _Candidate_SNPs: https://cellsnp-lite.readthedocs.io/en/latest/snp_list.html

.. code-block:: bash

  cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -R $REGION_VCF -p 20 --minMAF 0.1 --minCOUNT 20 --gzip

As shown in the above command line, we recommend filtering SNPs with <20UMIs
or <10% minor alleles for downstream donor deconvolution, by adding
``--minMAF 0.1 --minCOUNT 20``

Besides, special care needs to be taken when filtering PCR duplicates for scRNA-seq data by
including DUP bit in exclFLAG, for the upstream pipeline may mark each extra read sharing
the same CB/UMI pair as PCR duplicate, which will result in most variant data being lost.
Due to the reason above, cellsnp-lite by default uses a non-DUP exclFLAG value to include PCR
duplicates for scRNA-seq data when UMItag is turned on.


Mode 1b: well-based single cells or bulk
++++++++++++++++++++++++++++++++++++++++

Use ``-R`` but not ``-b`` to pileup well-based dataset (e.g., SMART-seq2) with given SNPs.

Require: one or multiple BAM/SAM/CRAM files (bulk or smart-seq), their according
sample ids (optional), and a VCF file for a list of common SNPs. BAM/SAM/CRAM files
can be input in comma separated way (``-s``) or in a list file (``-S``).

.. code-block:: bash

  cellsnp-lite -s $BAM1,$BAM2 -I sample_id1,sample_id2 -O $OUT_DIR -R $REGION_VCF -p 20 --cellTAG None --UMItag None --gzip

  cellsnp-lite -S $BAM_list_file -i sample_list_file -O $OUT_DIR -R $REGION_VCF -p 20 --cellTAG None --UMItag None --gzip

Set filtering thresholds according to the downstream analysis. Please add
``--UMItag None`` if your bam file does not have UMIs, e.g., smart-seq and bulk
RNA-seq.



Mode 2: pileup whole chromosome(s) without given SNPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Recommend filtering SNPs with <100UMIs or <10% minor alleles for saving space and speed up inference
when pileup whole genome: ``--minMAF 0.1 --minCOUNT 100``.

.. note::
   This mode may output false positive SNPs, for example somatic variants or falses caused by
   RNA editing. These false SNPs are probably not consistent in all cells within one individual, hence
   confounding the demultiplexing. Nevertheless, for species, e.g., zebrafish, without a good list of
   common SNPs, this strategy is still worth a good try.


Mode 2a: droplet based single cells without given SNPs
++++++++++++++++++++++++++++++++++++++++++++++++++++++

Don't use ``-R`` but use ``-b`` to pileup whole chromosome(s) without given SNPs 
for droplet-based dataset (e.g., 10x Genomics).

This mode requires inputting a single BAM/SAM/CRAM file with cell barcoded (add ``-b``):

.. code-block:: bash

  # 10x sample with cell barcodes
  cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -p 22 --minMAF 0.1 --minCOUNT 100 --gzip

Add ``--chrom`` if you only want to genotype specific chromosomes, e.g., ``1,2``, or ``chrMT``.

.. note::
   ``Mode 2a`` does joint calling and genotyping, but it is substantially slower than 
   calling first in a bulk manner by ``Mode 2b`` followed by genotyping in ``Mode 1a``. 
   Otherwise, it is handy for small chromosomes, e.g., mitochondrial.


Mode 2b: well-based single cells or bulk without SNPs
+++++++++++++++++++++++++++++++++++++++++++++++++++++

Don't use ``-R`` and ``-b`` to pileup whole chromosome(s) without given SNPs 
for well-based dataset (e.g., SMART-seq2).

This mode requires inputting one or multiple BAM/SAM/CRAM file(s) of bulk or smart-seq.

.. code-block:: bash

  # a bulk sample without cell barcodes and UMI tag
  cellsnp-lite -s $bulkBAM -I Sample0 -O $OUT_DIR -p 22 --minMAF 0.1 --minCOUNT 100 --cellTAG None --UMItag None --gzip

Add ``--chrom`` if you only want to genotype specific chromosomes, e.g., ``1,2``, or ``chrMT``.



Notes
-----

The ``Too many open files`` issue has been fixed (since v1.2.0). The issue is commonly
caused by exceeding the `RLIMIT_NOFILE`_ resource limit (ie. the max number of files allowed
to be opened by system for single process), which is typically 1024. Specifically, in the
case of ``M`` input files and ``N`` threads, cellsnp-lite would open in total about ``M*N`` files.
So the issue would more likely happen when large M or N is given. In order to fix it, cellsnp-lite
would firstly try to increase the limit to the max possible value (which is typically 4096) and
then use a fail-retry strategy to auto detect the most suitable number of threads (which could
be smaller than the original nthreads specified by user).

The command line option ``--maxFLAG`` is now deprecated (since v1.0.0), please use ``--inclFLAG`` and
``--exclFLAG`` instead, which are more flexible for reads filtering. You could refer to
the explain_flags_ page to easily aggregate and convert all flag bits into one integer.
One example is that the default exclFLAG value (without using UMIs) is 1796, which is
calculated by adding four flag bits: UNMAP (4), SECONDARY (256), QCFAIL (512) and DUP (1024).

.. _RLIMIT_NOFILE: https://man7.org/linux/man-pages/man2/getrlimit.2.html
.. _explain_flags: https://broadinstitute.github.io/picard/explain-flags.html

Full parameters
---------------
Here is a list of full parameters for setting (``cellsnp-lite -V`` always give the 
version you are using):

.. code-block:: html
  
  Usage: cellsnp-lite [options]
  
  Options:
    -s, --samFile STR    Indexed sam/bam file(s), comma separated multiple samples.
                         Mode 1a & 2a: one sam/bam file with single cell.
                         Mode 1b & 2b: one or multiple bulk sam/bam files,
                         no barcodes needed, but sample ids and regionsVCF.
    -S, --samFileList FILE   A list file containing bam files, each per line, for Mode 1b & 2b.
    -O, --outDir DIR         Output directory for VCF and sparse matrices.
    -R, --regionsVCF FILE    A vcf file listing all candidate SNPs, for fetch each variants.
                             If None, pileup the genome. Needed for bulk samples.
    -T, --targetsVCF FILE    Similar as -R, but the next position is accessed by streaming rather
                             than indexing/jumping (like -T in samtools/bcftools mpileup).
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
  
  Note that the "--maxFLAG" option is now deprecated, please use "--inclFLAG" or "--exclFLAG"
  instead. You can easily aggregate and convert the flag mask bits to an integer by refering to:
  https://broadinstitute.github.io/picard/explain-flags.html

