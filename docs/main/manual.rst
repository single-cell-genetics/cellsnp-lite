..
   Manual
   ======


Cellsnp-lite: Efficient Genotyping Bi-Allelic SNPs on Single Cells
==================================================================

.. contents:: Contents
   :depth: 2
   :local:


.. _manual Quick Usage:

Quick Usage (for Single-Cell or Bulk RNA-seq Data)
--------------------------------------------------
To install cellsnp-lite, please check :doc:`Installation</main/install>`. 
Once installed, check all arguments by typing ``cellsnp-lite -h``. 

Cellsnp-lite has four modes to support different genotype inputs and sequencing
platforms. They are summarised below:

.. csv-table:: Cellsnp-lite Modes
   :file: /tables/cellsnp-lite_modes.csv
   :header-rows: 1

Note, Mode 2b + 1a is an internal alternative to Mode 2a.

.. note::
   By default, cellsnp-lite would not output the file ``cellSNP.cells.vcf``
   which contains the genotypes in single cells. Please add ``--genotype``
   option for outputing this file.


Mode 1: pileup with given SNPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This mode genotypes single cells or bulk sample at a list of given SNPs, which 
could be common SNPs in human population (see 
:ref:`Candidate SNPs <data List of common SNPs>`), or
called heterouzygous variants (e.g., by cellsnp-lite ``Mode 2b`` on its own).


Mode 1a: droplet-based single cells
+++++++++++++++++++++++++++++++++++
Use both ``-R`` and ``-b`` to pileup droplet-based dataset (e.g., 10x Genomics)
with given SNPs.

Require:

* ``-s``: a single BAM/SAM/CRAM file, e.g., from CellRanger; 
* ``-b``: a list of cell barcodes, e.g., ``barcodes.tsv`` file in the 
  CellRanger directory ``outs/filtered_gene_bc_matrices/``; 
* ``-R``: a VCF file containing a list of SNPs.

This mode is recommended comparing to mode 2, if a list of common SNP is 
known, e.g., for human (see :ref:`Candidate SNPs <data List of common SNPs>`).

.. code-block:: bash

  cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -R $REGION_VCF -p 10 --minMAF 0.1 --minCOUNT 20 --gzip

As shown in the above command line, we recommend filtering SNPs with <20UMIs
or <10% minor alleles for downstream donor deconvolution, by adding
``--minMAF 0.1 --minCOUNT 20``


.. _manual Quick Usage Mode 1b:

Mode 1b: well-based single cells or bulk
++++++++++++++++++++++++++++++++++++++++
Use ``-R`` but not ``-b`` to pileup bulk or well-based dataset 
(e.g., SMART-seq2) with given SNPs.

Require:

* ``-s`` or ``-S``: one or multiple BAM/SAM/CRAM files (bulk or smart-seq), 
  specified either in comma separated way (``-s``) or in a list file (``-S``).
* ``-R``: a VCF file containing a list of SNPs.
* (Optional) ``-I`` or ``-i``: sample IDs of the BAM/SAM/CRAM files.

.. code-block:: bash

  # Set filtering thresholds according to the downstream analysis.
  cellsnp-lite -s $BAM1,$BAM2 -I SAMPLE_ID1,SAMPLE_ID2 -O $OUT_DIR -R $REGION_VCF -p 10 --cellTAG None --UMItag None --gzip

  cellsnp-lite -S $BAM_LIST_FILE -i SAMPLE_LIST_FILE -O $OUT_DIR -R $REGION_VCF -p 10 --cellTAG None --UMItag None --gzip

**Set filtering thresholds according to the downstream analysis.** Please add
``--UMItag None`` if your bam file does not have UMIs, e.g., smart-seq and bulk
RNA-seq.



Mode 2: pileup whole chromosome(s) without given SNPs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This mode genotypes single cells or bulk sample on whole chromosomes, without
given SNPs. 

Recommend filtering SNPs with <100UMIs or <10% minor alleles for saving space 
and speed up inference when pileup whole genome: 
``--minMAF 0.1 --minCOUNT 100``.

.. note::
   For mode2, by default it runs on chr1 to 22 on human. For mouse, you need 
   to specify it to 1,2,...,19 (replace the ellipsis).

.. warning::
   This mode may output false positive SNPs, for example somatic variants or 
   falses caused by RNA editing. 
   These false SNPs are probably not consistent in all cells within one 
   individual, hence could confound the downstream tasks such as 
   demultiplexing.
   Nevertheless, for species, e.g., zebrafish, without a good list of common
   SNPs, this strategy is still worth a good try.


Mode 2a: droplet based single cells without given SNPs
++++++++++++++++++++++++++++++++++++++++++++++++++++++
Don't use ``-R`` but use ``-b`` to pileup whole chromosome(s) without given 
SNPs for droplet-based dataset (e.g., 10x Genomics).

Require:

* ``-s``: a single BAM/SAM/CRAM file, e.g., from CellRanger; 
* ``-b``: a list of cell barcodes, e.g., ``barcodes.tsv`` file in the 
  CellRanger directory ``outs/filtered_gene_bc_matrices/``; 

.. code-block:: bash

  # 10x sample with cell barcodes
  cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 100 --gzip

Add ``--chrom`` if you only want to genotype specific chromosomes, 
e.g., ``1,2``, or ``chrMT``.

.. note::
   ``Mode 2a`` does joint calling and genotyping, but it is substantially 
   slower than calling first in a bulk manner by ``Mode 2b`` followed by 
   genotyping in ``Mode 1a``. 
   Otherwise, it is handy for small chromosomes, e.g., mitochondrial.


.. _manual Quick Usage Mode 2b:

Mode 2b: well-based single cells or bulk without SNPs
+++++++++++++++++++++++++++++++++++++++++++++++++++++
Don't use ``-R`` and ``-b`` to pileup whole chromosome(s) without given SNPs 
for bulk or well-based dataset (e.g., SMART-seq2).

Require: 

* ``-s`` or ``-S``: one or multiple BAM/SAM/CRAM files (bulk or smart-seq), 
  specified either in comma separated way (``-s``) or in a list file (``-S``).
* (Optional) ``-I`` or ``-i``: sample IDs of the BAM/SAM/CRAM files.

.. code-block:: bash

  # a bulk sample without cell barcodes and UMI tag
  cellsnp-lite -s $bulkBAM -I Sample0 -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 100 --cellTAG None --UMItag None --gzip

  # SMART-seq2 single cells
  cellsnp-lite -S $BAM_LIST_FILE -i SAMPLE_LIST_FILE -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 100 --cellTAG None --UMItag None --gzip

  # 10x scRNA-seq sample in a pseudo-bulk manner
  cellsnp-lite -s $BAM -O $OUT_DIR -p 10 --minMAF 0.1 --minCOUNT 20 --cellTAG None --UMItag UB --gzip

Add ``--chrom`` if you only want to genotype specific chromosomes, e.g., 
``1,2``, or ``chrMT``.


Advanced Usage
--------------
Cellsnp-lite supports data from various sequencing platforms, including
RNA-seq, DNA-seq, ATAC-seq, either in single-cell or bulk.

The default options of *cellsnp-lite* is set for 10x scRNA-seq data, i.e.,
``--cellTAG`` is set to ``CB`` and ``--UMItag`` is set to ``UB``.
However, it is very flexible to make *cellsnp-lite* to support data from other
platforms by changing a few options, 
mainly ``-b``, ``-i``, or ``-I`` to specify whether the data is in 
single-cell or bulk, 
and ``--cellTAG``, ``--UMItag`` to turn on or off cell and UMI tags.

Below shows some advanced usage of *cellsnp-lite* that incorporate different 
combinations of options.


.. _manual Advanced Usage Other Omics:

Processing other omics data
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Different omics data follow distinct experimental protocols, which leads
to the diversity of the output data format. 
For example, 10x 3' or 5' scRNA-seq data has both cell and UMI tags in the
BAM file, while 10x scDNA-seq and scATAC-seq data only have cell tag but not
UMI tag.

If the omics data has cell tag, set it in ``--cellTAG``, e.g., 
``--cellTAG CB`` for the ``CB`` tag for 10x scRNA-seq data.
Otherwise, please turn it off with ``--cellTAG None``.

If the omics data has UMI tag, set it in ``--UMItag``, e.g.,
``--UMItag UB`` for the ``UB`` tag for 10x scRNA-seq data.
Otherwise, please turn it off with ``--UMItag None``. 

We list options for some common omics data:

.. csv-table:: Cellsnp-lite Options for Various Omics
   :file: /tables/cellsnp-lite_options_for_various_omics.csv
   :header-rows: 1

If your data is not from platforms above, please choose proper ``--cellTAG``
and ``--UMItag`` values, e.g., by following the experimental protocols or 
by checking the BAM records with ``samtools view``.


Running in a pseudo-bulk manner
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cellsnp-lite supports bulk data, including bulk RNA-seq, DNA-seq, and
ATAC-seq, in :ref:`Mode 1b <manual Quick Usage Mode 1b>` or 
:ref:`Mode 2b <manual Quick Usage Mode 2b>`.

In certain scenarios, you may want to genotype in a psedo-bulk manner on your
single-cell data.
Then you may specify a single sample name in ``-I`` (e.g., ``-I Sample0``), 
not ``-b``. Note that you need to turn off the cell tag with 
``--cellTAG None``. 
As to the UMI tag (``--UMItag``), please choose a proper value based on your
data.

**Genotype 10x scRNA-seq data in a pseudo-bulk manner**

To genotype 10x scRNA-seq data in a pseudo-bulk manner with cellsnp-lite 
mode 1b (or mode 2b), it is recommended to subset the BAM file first, by 
extracting the alignment records with valid cell barcodes only. 
Here the valid cell barcodes are typically the cell barcodes stored in the 
cellranger output ``folder filtered_gene_bc_matrices``, which are the cells 
with high-quality sequencing data.

See also: issue #100.


Full Parameters
---------------
Please type ``cellsnp-lite`` or ``cellsnp-lite -h`` to see the list of full 
parameters.

General options
~~~~~~~~~~~~~~~
``-s, --samFile STR`` 
    Indexed BAM/CRAM file(s), comma separated multiple samples. 

``-S, --samFileList FILE`` 
    A file listing BAM/CRAM files, each per line.

``-O, --outDir DIR`` 
    Output directory for VCF and sparse matrices.

``-R, --regionsVCF FILE`` 
    A vcf file listing all candidate SNPs, for fetch each variants.

``-T, --targetsVCF FILE``
    Similar as ``-R``, but the next position is accessed by streaming rather
    than indexing/jumping (like ``-T`` in samtools/bcftools mpileup).

``-b, --barcodeFile FILE`` 
    A plain file listing all effective cell barcodes, e.g., the 
    ``barcodes.tsv`` file in the CellRanger directory
    ``outs/filtered_gene_bc_matrices/``.

    The barcodes in the ``-b`` file should match exactly with the
    string in the cell tag (``--cellTAG``), including the suffix 
    (e.g., ``-1``) if applicable.
    Otherwise, no UMIs or reads would be pileup and the output would be 
    empty.

``-i, --sampleList FILE`` 
    A list file containing sample IDs, each per line.

``-I, --sampleIDs STR``
    Comma separated sample IDs, e.g., ``"Sample_0"`` for single sample, or 
    ``"Sample_1,Sample_2,...,Sample_N"`` for multiple sample IDs.

``-V, --version``
    Print software version and exit.

``-h, --help``
    Show this help message and exit.


.. _manual Full Parameters Optional Arguments:

Optional arguments
~~~~~~~~~~~~~~~~~~
``--chrom STR``
    The chromosomes to use, comma separated. 
    Default is ``1 to 22`` (for human).

    You can also pileup specific chromosomes, e.g., ``1,2``, or ``chrMT``.

    If you want to pileup all chromosomes in **mouse data**, 
    you need to specify it to ``1,2,...,19`` (replace the ellipsis).

    **Chromosome names and order**

    For chromosome names: 
    the chromosome names specified by this option should match the ``@SQ``
    records in the SAM/BAM header, especially for mitochondrial chromosome,
    which has multiple names, such as ``chrM`` and ``chrMT``.
    You may check the ``@SQ`` records with ``samtools view -h``.

    Notably, *cellsnp-lite* would internally remove the "chr" 
    prefix (if available) for both BAM and VCF records after loading them.
    Therefore, users do not need to tweak the chromosome names in the two 
    files if they only differ in the "chr" prefix.

    Users do not need to sort the chromosomes as their order in both files 
    do not matter, as long as the BAM records have been sorted by coordinates,
    e.g,. with ``samtools sort``, and there is an BAM index (.bai) file.

``--cellTAG STR``
    Tag for cell barcodes, turn off with ``None``. 
    Default is ``CB``.

    .. note::
       Generally, you need to set this option to ``None`` if the input reads
       do not have cell barcodes, e.g., for 10x scDNA-seq or scATAC-seq data.
       Otherwise, no UMIs or reads would be pileup and the output would be
       empty.

``--UMItag STR``
    Tag for UMI: one of ``UB``, ``Auto``, ``None``. 
    Default is ``Auto``.

    For ``Auto`` mode, use ``UB`` if barcodes (``-b``) are inputted,
    otherwise use ``None``.
    The ``None`` mode means no UMI but read counts.

    .. note::
       For data without UMI, such as bulk RNA-seq, 10x scDNA-seq, 
       10x scATAC-seq, SMART-seq2 etc, please set ``--UMItag None``.
       Otherwise, all pileup counts will be zero.

``--minCOUNT INT``
    Minimum aggregated UMI or read count. 
    Default is ``20``.

    SNPs whose aggregated UMI (if ``--UMItag`` is not ``None``) or read 
    (otherwise) count is smaller than this value would be filtered and
    not outputted.

``--minMAF FLOAT``
    Minimum minor allele frequency. 
    Default is ``0.00``.

    The parameter ``minMAF`` is minimum minor allele frequency, which is 
    the minimum frequency of the allele with second highest read or UMI count 
    for a given SNP site. 

    This parameter can be used for SNP filtering. 
    See issue #77, #90, #93 for detailed discussions.

``-p, --nproc INT``
    Number of threads to use.
    Default is ``1``.

``-f, --refseq FILE``
    Faidx indexed reference sequence file. 
    If set, the real (genomic) reference allele (``REF``) extracted from 
    this file would be used for Mode 2 or for the missing REFs in the input 
    VCF for Mode 1.

    Without this option, cellsnp-lite mode 2 would take the allele with the 
    highest count as ``REF`` and the second highest as ``ALT``, 
    with little input information about the actual (genomic) reference. 
    This is different from mode 1, which uses the ``REF`` and ``ALT`` alleles
    specified in the input VCF.

    See also: issue #28.

``--genotype``
    If use, do genotyping in addition to counting.

    By default, cellsnp-lite would not output the file ``cellSNP.cells.vcf``
    which contains the genotypes (e.g., "0/0", "1/0", "1/1") in single cells. 
    Please add this option for outputing the file.

``--gzip``
    If use, the output VCF files will be zipped into ``BGZF`` format.
    Otherwise, the output VCF files would be plain files.

    Briefly, ``BGZF`` format is compatible with ``gzip``, while it is required
    for some popular HTS tools for indexing, e.g., ``bgzip``. 
    Please see details at https://www.htslib.org/doc/bgzip.html#BGZF_FORMAT.

``--printSkipSNPs``
    If use, the SNPs skipped when loading VCF will be printed. 
    This option is only used by developers for debug.

``--doubletGL``
    If use, keep doublet GT likelihood, i.e., GT=0.5 and GT=1.5. 
    This option will be marked as deprecated.


.. _manual Full Parameters Read Filtering:

Read filtering
~~~~~~~~~~~~~~
``--inclFLAG STR|INT``
    Required flags: skip reads with all mask bits unset.
    Default is ``""``.

``--exclFLAG STR|INT``
    Filter flags: skip reads with any mask bits set.
    Default is ``UNMAP,SECONDARY,QCFAIL`` (when use UMI) or 
    ``UNMAP,SECONDARY,QCFAIL,DUP`` (otherwise).

    You can easily aggregate and convert the flag mask bits to an integer at
    https://broadinstitute.github.io/picard/explain-flags.html

    .. note::
       Special care needs to be taken when filtering PCR duplicates for 
       10x scRNA-seq data by including ``DUP`` bit in ``--exclFLAG``, 
       for the upstream pipeline may mark each extra read sharing the same 
       CB/UMI pair as PCR duplicate, 
       which will result in most variant data being lost.
       Due to the reason above, cellsnp-lite by default uses a non-DUP 
       ``--exclFLAG`` value to include PCR duplicates for 10x scRNA-seq data 
       when ``--UMItag`` is turned on.

``--minLEN INT``
    Minimum mapped length for read filtering. 
    Default is ``30``.

    The mapped length is the number of reference positions that a read aligns 
    to, i.e., only count positions whose CIGAR operation is one of
    ``BAM_CMATCH``, ``BAM_CEQUAL``, ``BAM_CDIFF``.
    
    See also: `pysam::get_reference_positions() <https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_reference_positions>`.

``--minMAPQ INT``
    Minimum MAPQ for read filtering.
    Default is ``20``.

    MAPQ stands for mapping quality.

``--maxPILEUP INT``
    Deprecated. Please use ``--maxDEPTH``.
    
    .. note::
      This option was first introduced in cellsnp-lite v1.2.3, for setting
      a maximum pileup read count at a position per input file.
      It was designed to count those filtered reads as well, to be 
      distinguished from ``--maxDEPTH``, which was designed to exclude those
      filtered reads during counting.

      However, it seems the effect of ``--maxPILEUP`` deviates from the 
      original intention of designing it. 
      It has the same effect as ``--max-depth`` in ``bcftools mpileup``
      (and also ``--maxDEPTH`` in cellsnp-lite),
      which is expected to exclude filtered reads,
      since we used this ``--maxPILEUP`` value in ``bam_mplp_set_maxcnt()``.
      Therefore, we would like to mark this option as deprecated,
      and recommend using ``--maxDEPTH`` instead.

``--maxDEPTH INT``
    At a position, read maximally *INT* reads per input file,
    to avoid excessive memory usage.
    Default is ``0``.
    Note, ``0`` means highest possible value (currently ``INT_MAX``).

    It is expected to mimic the ``--max-depth`` in ``bcftools mpileup``.

``--countORPHAN``
    If use, do not skip anomalous read pairs.


Input
-----
Below are some details of the *cellsnp-lite* inputs.
Note that not all files listed below are required for *cellsnp-lite*.
Please look into section :ref:`Quick Usage <manual Quick Usage>` to check the 
required inputs for each mode of *cellsnp-lite*.


``Sequence alignments``
    BAM/CRAM file(s), specified via ``-s`` or ``-S``.

    Note that these files should be indexed, e.g., with ``samtools index``.

``A list of SNPs``
    VCF file, specified via ``-R`` or ``-T``.

    Note that this file is required for Mode 1, but not Mode 2.
    You may use either a list of genotyped SNPs (e.g., from bulk data), or
    common SNPs in population (we have pre-compiled a list of 7.4 million 
    common variants (AF>5%) for human, see 
    :ref:`List of common SNPs<data List of Common SNPs>` for details).

    **When genotypes for each individual is avaiable for demultiplexing.**

    You may use ``bcftools merge`` to make a combined VCF for all donors.

    **When the input VCF contains missing alleles.**
    
    Usually, the VCF should contain a list of heterozygous SNPs with valid
    ``REF`` and ``ALT`` alleles (i.e., ``REF`` and ``ALT`` should be one
    of ``'A'``, ``'C'``, ``'G'``, ``'T'``, and different from each other).
    In some special scenarios, the input ``REF`` or ``ALT`` could be empty,
    then *cellsnp-lite* can assign specific alleles to them.

    If the ``REF`` field in VCF is not provided, then *cellsnp-lite* will 
    extract the ``REF`` allele from the reference genome sequence
    automatically (FASTA file specified via ``-f``).
    If the ``ALT`` field in VCF is not provided, then *cellsnp-lite* will
    assign the allele (other than ``REF``) with the highest UMI/read counts as 
    the ``ALT``.

``A list of cell barcodes``
    Plain or gzip file, specified via ``-b``.

    One cell barcode per line in the file. 
    This file is required for genotyping single cells in data containing 
    cell tags, e.g., 10x scRNA-seq data.

``A list of sample IDs``
   Either a string specifying one or multiple sample IDs separated by comma 
   (``-I``), or a file listing sample IDs, each per line (``-i``).

   The sample ID(s) are required for genotyping in bulk data (single sample) or
   single cells in data without cell tags, e.g., SMART-seq2 data.


Output
------

Brief Introduction
~~~~~~~~~~~~~~~~~~
Cellsnp-lite outputs at least 5 files listed below 
(assuming ``--gzip`` option was used):

``cellSNP.base.vcf.gz``
    A VCF file listing genotyped SNPs and aggregated ``AD`` & ``DP`` 
    infomation (without ``GT``).

``cellSNP.samples.tsv``
    A TSV file listing cell barcodes or sample IDs.

``cellSNP.tag.AD.mtx``
    A file in "Matrix Market exchange formats", containing the allele depths 
    of the alternative (``ALT``) alleles.

``cellSNP.tag.DP.mtx``
    A file in "Matrix Market exchange formats", containing the sum of 
    allele depths of the reference and alternative alleles (``REF`` + ``ALT``).

``cellSNP.tag.OTH.mtx``
    A file in "Matrix Market exchange formats", containing the sum of 
    allele depths of all the alleles other than ``REF`` and ``ALT``.

Note, an additional VCF file ``cellSNP.cells.vcf.gz`` would be outputed 
if ``--genotype`` option was specified. 
This file contains genotyped SNPs and 
``AD`` & ``DP`` & genotype (``GT``) information for each cell or sample.


About REF and ALT alleles
~~~~~~~~~~~~~~~~~~~~~~~~~
The final output ``REF`` and ``ALT`` alleles are stored in the VCF files
``cellSNP.base.vcf.gz`` and ``cellSNP.cells.vcf.gz`` (if ``--genotype``
is used).

.. note::
   Cellsnp-lite was designed for bi-allelic SNPs.
   In its Mode 1, ``REF`` and ``ALT`` alleles are specified by user
   while in mode 2, ``REF`` and ``ALT`` are inferred from data as the alleles
   with highest and second highest read(UMI) counts.
   Therefore, in Mode 1, the ``REF`` or ``ALT`` in the reference VCF could be
   different from the major or minor allele inferred from data.
   For example, the ``ALT`` in VCF could be ``REF`` in the data.


Mode 1
++++++
In Mode 1, the ``REF`` and ``ALT`` alleles are expected to be specified in the
input VCF file (``-R`` or ``-T``).

1. When both ``REF`` and ``ALT`` are specified in input (most common scenario)
    The two alleles will be outputed as it is.

2. When ``REF`` is specified and ``ALT`` is missing in input
    The ``REF`` will be outputed as it is and the allele (other than ``REF``) 
    with the highest UMI/read counts will be assigned as the ``ALT``.

3. When ``REF`` is missing and ``ALT`` is specified in input
    When ``-f`` is used, the real genomic reference will be extracted from 
    FASTA file as ``REF`` and the allele (other than ``REF``) with the highest
    UMI/read count will be assigned as ``ALT``.

    Otherwise, *cellsnp-lite* would take the allele with the highest count 
    as ``REF`` and the second highest as ``ALT``.

    Note, the infered (output) ``ALT`` could be different from the 
    input ``ALT``.
    
4. When both ``REF`` and ``ALT`` are missing in input
    The same with point 3.


Mode 2
++++++
In Mode 2, the ``REF`` and ``ALT`` alleles are expected to be detected from
data.

When ``-f`` is used, the real genomic reference will be extracted from
FASTA file as ``REF`` and the allele (other than ``REF``) with the highest
UMI/read count will be assigned as ``ALT``.
Otherwise, *cellsnp-lite* would take the allele with the highest count
as ``REF`` and the second highest as ``ALT``.

