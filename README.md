# cardips-ipsc-eqtl

This repository contains the code for the CARDIPS iPSC eQTL study. This code
starts after most of the computationally intensive/time-consuming data
processing steps like alignment, variant calling, expression estimate, etc. The
repository mainly consists of Jupyter notebooks although there is also a
command line script for running EMMAX and a few small data files.

## Dependencies

If you want to run some of the code here, you'll need to install this project's
Python package `ciepy` included in this repo using either `python setup.py
install` or `python setup.py develop`. The notebooks in this repo will only run
as-is on the Frazer lab cluster though with some minor modifications (mainly
specifying some paths to data) it should be possible to run them elsewhere.

## Directories

### `ciepy`

Project-specific Python package needed to run some of the code here.

### `misc`

Miscellaneous files that are not output from notebooks.

### `notebook`

Jupyter notebooks for running analyses. Each notebook has corresponding
output directories in `output` and `private_output`.

### `output`

This directory contains one directory per Jupyter notebook. Output files here
will be shared on Figshare. Here is a description of some of the output
directories and files.

#### `eqtl_processing`

This directory contains subdirectories for each eQTL analysis I've run. Subdirs
of the form `eqt[number]` are results from using all samples. `eqtl01` is for
the initial eQTL analysis, `eqtl02` is the secondary analysis for all genes
that were significant in `eqtl01` using the lead variant as a covariate, etc.

Each subdir has the following files:

* `lead_variants.tsv` - This file has the lead variants for *all* genes
  including variants that are equally significant. Note that this file includes
  genes that were not significant. The `perm_sig` column indicates whether a
  gene was an eGene (i.e. had a significant eQTL).
* `lead_variants_single.tsv` - This file has a single lead variant for each
  gene. I choose randomly if there are ties. This file also contains *all*
  genes, not just significant genes.
* `pvalues.tsv` - This file has the permutation p-value for each gene.
* `qvalues.tsv` - This file has the permutation p-value, the q-value, and
  whether the gene is significant.
* `gene_variant_pairs.tsv` - This file has all significant variant-gene
  associations only for eGenes.
* `sig_lead_snvs_single.tsv` - This file has the most significant SNV per gene
  for eGenes. I choose randomly to break ties and drop genes that don't have a
  significant SNV.
* `independent_lead_snvs.tsv` - This file is created by LD pruning the variants
  in `sig_lead_snvs_single.tsv`.
* `independent_lead_snvs.bed` - Bed file for `independent_lead_snvs.bed`.
* `all_results_sorted.tsv.gz` - This file contains the output for every
  variant-gene pair tested by EMMAX. The file is sorted by position and tabix
  indexed.
* `top_results_sorted.tsv.gz` - This file contains the most significant
  association for each variant tested by EMMAX. The file is sorted by position
  and tabix indexed.

I also ran a few different variations of the eQTL analysis. `unrelated_eqtls01`
contains results from running the eQTL analysis using only the 131 unrelated
subjects. `no_peer_no_std_norm01` contains results from running the analysis
using log TPM expression values without transformation and including batch,
sex, and donor age as covariates.

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
