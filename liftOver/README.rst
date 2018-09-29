====================================
A wrap of UCSC liftOver for vcf file
====================================

The liftOver function from UCSC_ is a very useful tool, though it only supports 
BED_ format. An alternative software, CrossMap_, allows multiple formats, 
including VCF_, but occasionally fails. Therefore, I wrapped some functions for 
the UCSC liftOver software to support input/output in VCF format.

This Python wrap function contains three steps, 
1. convert the position in VCF into BED format
2. use UCSC liftOver on this input BED file to generated the new BED file
3. convert to the new BED file into VCF file.

Note, a tricky thing here is that the UCSC liftOver ouputs two files for mapped 
and unmapped positions respectively in the same order as the input VCF (I 
believe). Then these two files will be parsed for generate the new VCF file.

Also, for you liftOver between genome builds, you need the chain file to link 
two builds. You can find most chain file in UCSC 

.. _UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/
.. _BED: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _VCF: http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/
.. _CrossMap: http://crossmap.sourceforge.net/

Usage
=====

Use Python to execute the liftOver_vcf.py_ as the following comman line:

.. code-block:: bash

  python liftOver_vcf.py -c hg19ToHg38.over.chain.gz -i my_input.hg19.vcf.gz -o my_output.hg38.vcf.gz


.. _liftOver_vcf.py: https://github.com/huangyh09/cellSNP/blob/master/liftOver/liftOver_vcf.py

Resources
=========
**Software**

* UCSC liftOver online version: https://genome.ucsc.edu/cgi-bin/hgLiftOver
* UCSC liftOver local software (we need here): https://genome-store.ucsc.edu

**Chain files**

* UCSC chain files (all): http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/
* Common files: 
  `hg38ToHg19.over.chain.gz <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz>`_, 
  `hg19ToHg38.over.chain.gz <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz>`_
  
