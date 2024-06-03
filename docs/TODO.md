## TODO List

### Input
-	Wrap codes into functions which open BAM and their index files.


### Output
- Output VCF header according to input BAM header.
- Output optionally qual values/letters to mtx file.


### Implementation
- Add fetch and pileup sub-commands?
  * merge `csp_fetch()` and `csp_pileup()` since the two functions share
  most codes?

- Support calling germline SNPs for multiple bam files?
  * in bulk mode

- Update `-T` method
  * use `qsort` & linear search instead of `regidx_t` of htslib;
  * note the bug of qsort in lower version of `glibc`, 
    refer to https://sourceware.org/bugzilla/show_bug.cgi?id=11655
  * as the linear searching assume that the BAM has been sortted by start pos,
    then how to deal with the partly aligned reads and one read in paired 
    reads aligned to query chrom.

- Try multi-process (process pool) for multi input samples.

- Pileup SNPs within specific genomic regions
  * this feature should be useful when there are only a few target regions.
  * it can improve the overall running efficiency, particularly for mode 2 
    (or mode 1-T), e.g., by spliting all chromosomes equally based on their 
    lengths to make full advantage of multi-threading, as now only one thread 
    is used when only one chrom is inputted (e.g., chrom M).


### Data Structure
- Save memory for snp list
  * use fixed 32-bit index for chrom string.
  * use shared chrom string for all snps in that chrom, like the way in the 
    `develop` branch.

- Use optional sparse matrices tags with the help of function pointers.

- Wrap up the steps in the whole workflow into proper data structures.
  * e.g., the snp iteration, the read extraction from input files 
    (especially the redundant `hts_itr_xxx`), mplp statistics, temporary file 
    creation & output, etc. 

- Update `mplp_t::su` and `plp_t::hug`.

- Separate `htsFile` from `csp_bam_fs` as it cannot be shared among threads.

- Try using `multi_iter` fetching method of BAM/SAM/CRAM for multi regions 
  (SNPs) if it can in theory speed cellsnp up.

- Improve the `jfile_t` structure, for example, adding @p `is_error`.

- Improve the `JMEMPOOL` structure, for example, adding @p `base_init_f`.


### Reads operations
- Max per base quality
  * Add `--max-base-qual` option to set max per base quality for genotyping;
  * in htslib, per base qualities are stored in the `Phred` scale with no 
    `+33` offset.
  * `get_qual_vector()` in `mplp.c`.

- UMI collapsing
  * consistency correction could be done in UMI groups with the help of 
    @p `pu` & @p `pl` inside `mplp` structure.
  * update `map_ug_t` first!

- Deal with the problem that some UMIs have the letter 'N'.


### General running operations
- Improve logging
  * add time str to show_progress (percentage of processed SNPs) so that 
    we can know how long it has been from the last (output file) update.


### Docs
-	Add tutorials 
  * add `Post-hoc Analysis` notebooks such as “how to post-filtering SNPs”. 
  * to use files outside the root source directory of sphinx (readthedocs),
    refer to the 
    [conf.py](https://github.com/single-cell-genetics/vireo/blob/master/doc/conf.py)
    in vireo repo for copying notebooks during building time. 
  * see also: issue 90 and 108.


### Tests
- Write test scripts for some key functions.


## Discussion

### Cmdline options
Cellsnp-lite now only considers biallele (related to -f/--refseq).
It may have unexpected results affected by `OTH` (other) alleles.
We may only consider the `REF` and `ALT` alleles (either specified or 
inferred), but not `OTH` alleles, when calculating `MAF`.

**Case 1**
In mode 1, the output could contain Homozygous SNV even with `--minMAF 0.3`.
For example, assuming one input SNV has `REF/ALT - 'A'/'C'`, 
while the two real alleles are `'C'/'G'` with AF 0.6/0.4, 
then this SNV would pass `--minMAF 0.3` and the genotype is `1/1`
(as `REF` is `'A'`, `ALT` is `'C'`), 
while the real genotype should be `1/2` (as two alt alleles `'C','G'`).

However, SNVs of this kind are not so many in practice?
In a recent case, only 158 out of 133k SNVs.

**Case 2**
A SNP with `AD=9;DP=9;OTH=1;` will not be filtered when 
`minMAF=0.1; minCOUNT=10`.


### Multi-threading by spliting input BAM files
Processing well-based data with multiple threads could lead to the
`too many files open` issue when the number of input cells (BAM files) is 
large.
One possible solution to address the issue is to split input cells instead of
splitting SNPs/Chromosomes when using multi-threading.
However, its implementation is non-trivial because

- there will be too large intermediate results (and IO operations will be 
  heavy) as whether one SNP is filtered or not is not determined until all 
  files are processed; 
- users could merge small BAM files to reduce the number of files, and add
  tags, e.g., `RG` or `CB`, to distinguish the read source.