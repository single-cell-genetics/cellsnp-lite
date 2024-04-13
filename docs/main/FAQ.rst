..
   FAQ
   ===

..
   What is cellsnp-lite?
   What is the input of cellsnp-lite?
   What is the output of cellsnp-lite?
   How does the SNP-filtering options work?


FAQ
===

.. contents:: Contents
   :depth: 2
   :local:


Top Questions
-------------

..
   Troubleshooting
   ---------------

Why the output files are empty, no SNPs genotyped?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If there is no read captured from the bam file, there can be multiple reasons:

Improper ``--UMItag`` setting
    if your data has UMI tag, please use the right value 
    (default value for ``--UMItag`` is ``UB``).
    Otherwise, please set ``--UMItag None``, e.g., for SMART-seq2,
    bulk RNA-seq, 10x scDNA-seq, scATAC-seq data.
    See also issue #26, #44.

Improper ``--cellTAG`` setting
    if your data has cell tag, please use the right value 
    (default value for ``--cellTAG`` is ``CB``).
    Otherwise, please set ``--cellTAG None``, e.g., for SMART-seq2,
    bulk RNA-seq data.

Unmatched cell barcodes
    please check whether the cell barcode (string) in the cell tag of BAM file
    can match exactly with the barcodes in the input file (``-b``), 
    including the possible suffix (e.g., ``-1``).
    See also issue #44.

Unmatched chromosome names
    the input chromosome names, either in the input VCF file or specified
    by ``--chrom`` option, should match the ``@SQ``
    records in the SAM/BAM header, especially for mitochondrial chromosome,
    which has multiple names, such as ``chrM`` and ``chrMT``.
    Please use the right chromosome names, you may check the ``@SQ`` 
    records with ``samtools view -h``.

    Notably, *cellsnp-lite* would internally remove the "chr"
    prefix (if available) of all chromosome names, including the names
    specified by ``--chrom`` option and the ones in the input BAM and
    VCF records.
    Therefore, users do not need to tweak the chromosome names in the option
    and the two files if they only differ in the "chr" prefix.


Input and Output
----------------

Can cellsnp-lite be used for other omics, such as scDNA-seq or scATAC-seq data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Yes. 
Cellsnp-lite supports multiple omics data from various sequencing 
platforms. 
You may use proper options for different omics, mainly the ``--cellTAG``
for cell tag and ``--UMItag`` for UMI tag, e.g., please set ``--UMItag None``
for 10x scDNA-seq or scATAC-seq data.

See :ref:`Processing other omics data <manual Advanced Usage Other Omics>` 
for details and more options.

See also: issue #26, #44.


Can cellsnp-lite be used for spatial transcriptomics data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Yes.
Efforts have been made to explore using *cellsnp-lite* to calculate the
B allele frequency (BAF) as one input signal for Numbat_ to calculate the
CNV probability per spot in spatial transcriptomics
(`Arora et al, 2023`_).


Can cellsnp-lite be used for mouse data?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Yes.

For Mode 2, by default it runs on chr1 to 22 on human. 
For mouse, you need to specify it to ``1,2,...,19 (replace the ellipsis)``.

On the other hand, if it is lab straint with genotype available, e.g., 
from the mouse genome project (ftp://ftp-mouse.sanger.ac.uk/current_snps), 
you can also run Mode 1.

See also: issue #3.


Why file "cellSNP.cells.vcf.gz" is missing in the output folder?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``--genotype`` option should be added for *cellsnp-lite* to output the 
file ``cellSNP.cells.vcf.gz``.

See also: issue #4, #22.


Running Mode
------------

Mode 2a is too slow on large dataset, what should I do?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Mode 2a is more suitable for small datasets. 
For large datasets, you may try ``Mode 2b + Mode 1a``. 
Mode 2a does joint calling and genotyping, but it is substantially slower 
than calling first in a bulk manner by Mode 2b followed by genotyping in 
Mode 1a. 
To speed up, you may try ``--minMAF 0.1 --minCOUNT 100`` options in both modes.

See also: issue #70.


My dataset is not large, why it takes so long time?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To speedup, you may

* Check whether the cell barcodes are "filtered", i.e., from 
  ``filtered_gene_bc_matrices`` instead of from ``raw_gene_bc_matrices`` 
  in the cellranger output folder (update ``-b``);
* Try to use the SNP list from ``AF5e2`` VCF file instead of ``AF5e4`` in 
  this `human SNP list`_ folder (update ``-R``);
* Use more threads or cores (update ``-p``).

See also: issue #78.


Pileup and Genotype Results
---------------------------

Why the reference alleles (REF) are different from FASTA file (in Mode 2)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cellsnp-lite Mode 2 takes the allele with the highest count as ``REF`` and 
the second highest as ``ALT`` by default. 
Therefore, neither allele is necessarily identical to the actual (genomic)
reference in Mode 2.
This is different from Mode 1, which uses the ``REF`` and ``ALT`` alleles 
specified in the input VCF. 

However, since v1.2.2, *cellsnp-lite* has the ``-f`` or ``--refseq`` option
to extract the real (genomic) reference allele from FASTA file as ``REF``,
and assign the allele (other than ``REF``) with the highest UMI/read counts 
as the ``ALT``.

See also: issue #28.


How to check the real reads cover an SNP in a cell?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can extract the reads by *samtools* and then view them in *IGV*.

To extract reads covering a SNP and output to a BAM file 
(**assuming the SNP position is chr1:100000**):

.. code-block:: bash

  samtools view -h -b  "input_BAM"  chr1:100000  >  "output_BAM"
  samtoos index "output_BAM"

If you only want to extract SNP reads in specific cell 
(**assuming cell barcode is XXX-1 and cell tag is CB**):

.. code-block:: bash

  samtools view -h -b  -d CB:XXX-1  "input_BAM"  chr1:100000  >  "output_BAM"
  samtoos index "output_BAM"

Then you can load the output BAM file above into *IGV* to view the reads.

See also: issue #107.


Why the output allele counts are different from IGV view?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IGV would filter some reads by default, which could lead to the difference
in allele counts with cellsnp-lite output.
The allele counts should be the same if given the same read filtering settings.

You may refer to the question on this page
``Which reads will be filtered by cellsnp-lite?`` and 
``Preferences -> Alignments`` for read filtering settings of *cellsnp-lite*
and IGV, respectively.

See also: issue #95.


Why the output read counts are different from samtools view?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The inconsistency of read counts between *samtools* and *cellsnp-lite* is 
probably due to the different filtering settings of the two tools, 
e.g., by default, *cellsnp-lite* will filter some low-quality reads 
(please check ``--exclFLAG`` option) while samtools do not. 
To make the filtering settings the same, you can use ``-F`` option in 
*samtools view*.

See also: issue #107.


Why the output read counts are different from cellSNP (python version)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The difference in read counts is probably because the two methods used 
different read filtering settings, especially in Mode 2.

In Mode 2, cellSNP (actually the dependency ``pysam.pileup()``) has a default 
limitation that the ``max_depth`` (i.e., max pileup-ed read count) 
is ``8000``, 
However, cellsnp-lite does not have this ``max_depth`` limitation by default, 
it will pileup as many reads as possible. 
You may try using the same read filtering settings for both cellsnp-lite and
cellSNP, to make their read counts highly concordant in Mode 2.

See also: issue #33.


About Implementation
--------------------

Which SNPs will be filtered by cellsnp-lite?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cellsnp-lite has a few options for SNP filtering.
By default, SNP will be filtered if

* its aggregated UMI (if ``--UMItag`` is not ``None``) or read (otherwise) 
  count is <20 (``--minCOUNT``);
* its minor allele frequency (the frequency of the allele with second highest
  read or UMI count) is <0 (``--minMAF``).

See :ref:`Optional Arguments <manual Full Parameters Optional Arguments>`
in manual for details and more options.


Which reads will be filtered by cellsnp-lite?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cellsnp-lite has a few options for read filtering. 
By default, read will be filtered if

* it does not contain target cell tag (if set in ``--cellTAG``) or 
  its cell barcode is not in the input barcode list (``-b``);
* it does not contain target UMI tag (if set in ``--UMItag``);
* any mask bits is set in SAM FLAG: 
  ``UNMAP``, ``SECONDARY``, ``QCFAIL`` (when use UMI)
  or ``UNMAP``, ``SECONDARY``, ``QCFAIL``, ``DUP`` (otherwise).
* its mapped length is <30 (``--minLEN``);
* its mapping quality MAPQ is <20 (``--minMAPQ``);
* total pileup read count per input file is >INT_MAX (``--maxDEPTH``);
* it is not mapped in proper pairs (``--countORPHAN``).

See :ref:`Read Filtering <manual Full Parameters Read Filtering>` 
in manual for details and more options.

See also: issue #25.


What is the computational efficiency of cellsnp-lite?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In theory, the computational complexity (i.e., running time) of cellSNP-lite 
is ``O(n)`` for number of variants and ``O(n*log(n))`` for number of cells, 
which means it is more sensitive to the cell counts.

For `human SNP list`_, we suggest using the version with ``AF5e2`` 
(i.e., AF>5%, 7.4M SNPs), instead of ``AF5e4`` (i.e., AF>0.05%, 36.6M SNPs).


Is it fine to have multiple primary alignments in BAM file?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
IMPO, the "multi-primary" strategy, in which multiple alignments with the 
best score are labeled as primary, should be fine for downstream tasks 
if the fraction of the "extra" primary alignments is low.

Generally, we recommend to use "single-primary" strategy for genotyping,
in which only one alignment with best alignment score is labelled as primary
and the rest as secondary.

See detailed discussion in issue #39.


Command Line Options
--------------------

What is ``--minMAF``, can I post-filter SNPs based on ``minMAF``?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Cellsnp-lite was designed for bi-allelic SNPs. 
In its Mode 1, ``REF`` and ``ALT`` alleles are specified by user 
while in mode 2, ``REF`` and ``ALT`` are inferred from data as the alleles
with highest and second highest read(UMI) counts. 
Therefore, in Mode 1, the ``REF`` or ``ALT`` in the reference VCF could be
different from the major or minor allele inferred from data. 
For example, the ``ALT`` in VCF could be ``REF`` in the data.

In cellsnp cmdline (for both Mode 1 and 2), ``MAF`` is always caculated as 
the fraction read(UMI)_count_of_minor_allele / read(UMI)_count_of_all_alleles,
where the minor allele is the allele with second highest read(UMI) count 
inferred from data. 
See also issue #77.

Therefore, in Mode 1, post-filtering SNPs based on the minimum allele 
frequency of the ``REF`` and ``ALT`` alleles in VCF file could be different 
from filtering SNPs with ``--minMAF`` in the cellsnp cmdline, 
for a small subset of SNPs whose major allele (with highest read/UMI count) 
or minor allele (second highest) is neither ``REF`` or ``ALT`` allele 
but one of the ``OTH`` alleles. 
See also issue #90.

The number of SNPs whose major or minor allele is one of the ``OTH`` alleles 
is expected to be quite small (in Mode 1), given the input reference VCF is 
reliable (e.g., with common SNPs compiled from 1000 genome project), hence 
should have limited influence on downstream donor deconvolution.

See also: issue #77, #90, #93.


Downstream Demultiplexing with Vireo
------------------------------------

Cellsnp-lite output not working with Vireo using mitochondrial SNPs?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The large mitochondrial read counts in cellsnp-lite output makes it more 
likely for vireo to reach local optima so that the parameters of donors become
the same and hence vireo cannot assign the cells to certain donor.

Besides, vireo is designed for nuclear SNVs. 
For mito SNVs, you may want to try this `vireo Mito tutorial`_, 
which was used by MQuad_.
Note that the duplicate reads should probably be removed beforehand, 
if there are no UMIs in your data.

See also: issue #33.


How can I use genotypes for each individual when avaiable?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You may use ``bcftools merge`` to make a combined VCF for all donors.

See also: issue #21, #100, #106.


.. _Arora et al, 2023: https://doi.org/10.1038/s41467-023-40271-4 
.. _human SNP list: https://sourceforge.net/projects/cellsnp/files/SNPlist/
.. _MQuad: https://github.com/single-cell-genetics/MQuad
.. _Numbat: https://github.com/kharchenkolab/numbat
.. _vireo Mito tutorial: https://vireosnp.readthedocs.io/en/latest/vireoSNP_clones.html

