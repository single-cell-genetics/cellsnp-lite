
This directory stores scripts and data for processing the output
of [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
in a post-hoc manner.

We initialized such scripts and data in 
[scripts/utils](https://github.com/single-cell-genetics/cellsnp-lite/tree/master/scripts/utils) 
dir and then set up this new dir specifically for post-hoc analysis.
To keep the file URL stable, we will keep the scripts and data in 
`scripts/utils` unchanged, and put new files here.

The `subset_with_minAD_issue108.ipynb` is a jupyter notebook giving an example
of how to load cellsnp-lite output, subset SNPs, and the save the subset data.
Specifically, this example shows how to filter SNPs whose aggregated AD is 
less than `minAD` (issue 108).

