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
  variants that are equally significant.
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

#### `input_data`

The files in this directory are created from data only available on the Frazer
lab server. These files serve as input files throughout the project. Some of
the files are:

* `cnvs.tsv` - CNVs discovered between reprogrammed stem cells and blood cells.
* `mbased_major_allele_freq.tsv` - Major allele frequency estimate (i.e. the
percent of transcripts estimated to come from the dominantly-expressed allele)
from MBASED. If a value is missing, that means MBASED didn't have enough
variants or coverage to use the gene in that sample.
* `mbased_p_val_ase.tsv` - ASE p-value from MBASED. If a value is missing, that
means MBASED didn't have enough variants or coverage to use the gene in that
sample.
* `mbased_p_val_het.tsv` - ASE heterogeneous p-value from MBASED. This means
that the gene seems to have some sort of transcript regulation where some het
variants are imbalanced and others aren't. If a value is missing, that means
MBASED didn't have enough variants or coverage to use the gene in that sample.
* `rnaseq_metadata.tsv` - Metadata for RNA-seq samples.
* `rsem_tpm.tsv` - RSEM TPM values for all samples used in this study.
* `star_logs.tsv` - Contents of STAR `Log.final.out` files for all samples used
in this study.
* `subject_metadata.tsv` - Metadata for each subject (person) in this study.
* `wgs_metadata.tsv` - Metadata for each WGS samples used in this study.

### `private_output`

This directory contains one directory per Jupyter notebook. Output files here
contain sensitive information or are too large to be shared on Figshare.

### `scripts`

Contains scripts that are run from the command line.
