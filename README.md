# cardips-ipsc-eqtl

CARDIPS iPSC eQTL study.

## Dependencies

You'll need to install this project's Python package `ciepy` included in this
repo using either `python setup.py install` or `python setup.py develop`.

## Directories

### `ciepy`

Project-specific Python package.

### `misc`

Miscellaneous files that are not output from notebooks.

### `notebook`

Jupyter notebooks for running analyses. Each notebook has corresponding
output directories in `output` and `private_output`.

### `output`

This directory contains one directory per Jupyter notebook. Output files here
will be shared on Figshare. Here is a description of a selection of output
directories and files.

#### `eqtl_processing`

This directory contains subdirectories for each eQTL analysis I've run. Subdirs
of the form `eqt[number]` are results from using all samples. `eqtl01` is for
the initial eQTL analysis, `eqtl02` is the secondary analysis for all genes
that were significant in `eqtl01` using the lead variant as a covariate, etc.
Each subdir has the following files:

* `lead_variants.tsv` - This file has the lead variants for all genes including
  variants that equally significant.
* `lead_variants_single.tsv` - This file has a single lead variant for each
  gene. I choose randomly if there are ties.
* `pvalues.tsv` - This file has the permutation p-value for each gene.
* `qvalues.tsv` - This file has the permutation p-value, the q-value, and
  whether the gene is significant.
* `gene_variant_pairs.tsv` - This file has all significant variant-gene pairs.
* `sig_lead_snvs_single.tsv` - This file has the most significant SNV per gene
  for genes that were significant. I choose randomly to break ties and drop
  genes that don't have a significant SNV.
* `independent_lead_snvs.tsv` - This file is created by LD pruning the variants
  in `sig_lead_snvs_single.tsv`.
* `independent_lead_snvs.bed` - Bed file for `independent_lead_snvs.bed`.

### `private_output`

This directory contains one directory per Jupyter notebook. Output files here
contain sensitive information or are too large to be shared on Figshare.

### `scripts`

Contains scripts that are run from the command line.
