# TODO List

## Input
### BAM Files
#### Wrap codes into functions which open BAM and their index files.



## Output
### Matrix Files
#### Output optionally qual values/letters to mtx file.

### VCF Files
#### Output VCF header according to input BAM header.



## Implementation
### Commands and Options
#### Add fetch and pileup sub-commands.
Merge `csp_fetch()` and `csp_pileup()` since the two functions 
share most codes?


### Data Structure
#### Save memory for snp list.
- Use fixed 32-bit index for chrom string.
- Use shared chrom string for all snps in that chrom, like the way in the 
  `develop` branch.

#### Use optional sparse matrices tags with the help of function pointers.

#### Wrap up the steps in the whole workflow into proper data structures.
- e.g., the snp iteration, the read extraction from input files 
  (especially the redundant `hts_itr_xxx`), mplp statistics, temporary file 
  creation & output, etc. 

#### Update `mplp_t::su` and `plp_t::hug`.

#### Separate `htsFile` from `csp_bam_fs` as it cannot be shared among threads.

#### Try using `multi_iter` fetching method 
Try using `multi_iter` fetching method of BAM/SAM/CRAM for multi regions 
(SNPs) if it can in theory speed cellsnp up.

#### Improve the `jfile_t` structure, for example, adding @p `is_error`.

#### Improve the `JMEMPOOL` structure, for example, adding @p `base_init_f`.


### Efficiency
#### Reduce the chance that one thread gets stuck.
How to reduce the chance that one thread gets stuck for too long time 
due to super high coverage of certain SNPs?

One strategy is to shuffle input SNPs first and then push them into threads.

Additionally, try splitting input SNPs into smaller batches (for -R option) 
and push them into thread pools.
Notably, pushing SNPs into a thread for pileup requires extra time overhead, 
e.g., re-opening BAM file(s) and switching threads,
hence the number of SNPs in each batch should be reasonable such that the
extra time overhead of pushing SNPs is acceptable.
Assuming there are N input SNPs and M threads, one proposal is to assign
sqrt(N/M) SNPs to each batch.

#### Update `-T` method.
- Use `qsort` & linear search instead of `regidx_t` of htslib;
- Note the bug of qsort in lower version of `glibc`, 
  refer to https://sourceware.org/bugzilla/show_bug.cgi?id=11655
- As the linear searching assume that the BAM has been sortted by start pos,
  then how to deal with the partly aligned reads and one read in paired 
  reads aligned to query chrom.


### Feature Requests
#### Support calling germline SNPs for multiple bam files.
for bulk mode.

#### Pileup SNPs within specific genomic regions.
- This feature should be useful when there are only a few target regions.
- It can improve the overall running efficiency, particularly for mode 2 
  (or mode 1-T), e.g., by spliting all chromosomes equally based on their 
  lengths to make full advantage of multi-threading, as now only one thread 
  is used when only one chrom is inputted (e.g., chrom M).

#### Try multi-process (process pool) for multi input samples.
Processing well-based data with multiple threads could lead to the
`too many files open` issue when the number of input cells (BAM files) is 
large.

**Multi-threading by spliting input BAM files**

One possible solution to address the issue is to split input cells instead of
splitting SNPs/Chromosomes when using multi-threading.
However, its implementation is non-trivial because

- there will be too large intermediate results (and IO operations will be 
  heavy) as whether one SNP is filtered or not is not determined until all 
  files are processed; 
- users could merge small BAM files to reduce the number of files, and add
  tags, e.g., `RG` or `CB`, to distinguish the read source.

**Alternatively, use multi-processing**

Then the `too many files open` issue could be fixed as each process 
is independent, whereas the overall memory usage could be increased.


### General
#### Improve logging.
- Add time str to show_progress (percentage of processed SNPs) so that 
  we can know how long it has been from the last (output file) update.


### Reads Operations
#### Max per base quality.
- Add `--max-base-qual` option to set max per base quality for genotyping;
- In htslib, per base qualities are stored in the `Phred` scale with no 
  `+33` offset.
- Refer to `get_qual_vector()` in `mplp.c`.

#### UMI collapsing.
- Consistency correction could be done in UMI groups with the help of 
  @p `pu` & @p `pl` inside `mplp` structure.
- Update `map_ug_t` first!

#### Deal with the problem that some UMIs have the letter 'N'.



## Docs
### Tutorials 
#### Add `Post-hoc Analysis` notebooks such as “how to post-filtering SNPs”. 
To use files outside the root source directory of sphinx (readthedocs),
refer to the 
[conf.py](https://github.com/single-cell-genetics/vireo/blob/master/doc/conf.py)
in the `vireo` repo for copying notebooks during building time. 
See also: issue 90 and 108.



## Scripts
### Tests
#### Write test scripts for some key functions.



## Discussion
### MAF
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
`minMAF=0.1; minCOUNT=10`, as the `OTH` is also counted when calculating
aggregated counts (related to `minCOUNT`).