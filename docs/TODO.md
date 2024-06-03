
## TODO List

### Input


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
  * in htslib, per base qualities are stored in the Phred scale with no +33 offset.
  * `get_qual_vector()` in `mplp.c`.

- UMI collapsing
  * consistency correction could be done in UMI groups with the help of 
    @p `pu` & @p `pl` inside `mplp` structure.
  * update `map_ug_t` first!

- Deal with the problem that some UMIs have the letter 'N'.


### Docs

### Tests
- Write test scripts for some key functions.


## Discussion
Cellsnp-lite now only considers biallele (related to -f/--refseq).
It may have unexpected results affected by `OTH` (other) alleles.
We may only consider the `REF` and `ALT` alleles (either specified or 
inferred), but not `OTH` alleles, when calculating `MAF`.

#### Case 1
In mode 1, the output could contain Homozygous SNV even with `--minMAF 0.3`.
For example, assuming one input SNV has `REF/ALT - 'A'/'C'`, 
while the two real alleles are `'C'/'G'` with AF 0.6/0.4, 
then this SNV would pass `--minMAF 0.3` and the genotype is `1/1`
(as `REF` is `'A'`, `ALT` is `'C'`), 
while the real genotype should be `1/2` (as two alt alleles `'C','G'`).

However, SNVs of this kind are not so many in practice?
In a recent case, only 158 out of 133k SNVs.

#### Case 2
A SNP with `AD=9;DP=9;OTH=1;` will not be filtered when 
`minMAF=0.1; minCOUNT=10`.