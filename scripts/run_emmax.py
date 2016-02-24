import argparse
import gzip
import os
import random
random.seed(20150605)
import shutil
import subprocess
import sys

import pandas as pd

def _cov(covf, covariates, samples):
    """
    Make .cov file for EMMAX. The first column is individual ID and the
    subsequent files are covariates.
    
    Parameters
    ----------
    covf : str
        Path for output cov file.

    covariates : str
        Path to input covariates file. This file will be checked and written in
        the correct order for EMMAX.

    samples : str
        Pandas Series with samples (in correctly sorted order) as values.

    """
    c = pd.read_table(covariates, header=None, index_col=0)
    assert len(set(c.index) & set(samples.values)) == samples.shape[0], \
            'Covariates and samples files do not match.'
    c.ix[samples.values].to_csv(covf, sep='\t', header=False)

def _ind(indf, samples):
    """
    Make .ind file for EMMAX. The first column is individual ID and the second
    column is the position of the sample in VCF file. Since I make the VCF based
    on the input ind file, I just use the order of the samples to write the temp
    ind file.
    
    Parameters
    ----------
    indf : str
        Path for output ind file.

    samples : str
        Pandas Series with samples (in correctly sorted order) as values.

    phenos : pandas.DataFrame
        DataFrame with phenotype values.

    permut : list
        A list of indices to permute the phenotype data before writing it.

    """
    se = pd.Series(range(1, samples.shape[0] + 1), index=samples.values)
    se.to_csv(indf, header=False, sep='\t')

def _phe(phenof, phenos, permut=False):
    """
    Make .phe file for EMMAX. The first column is individual ID and the second
    column is the phenotype value.
    
    Parameters
    ----------
    phenof : str
        Path for output phe file.

    phenos : pandas.Series
        Pandas Series with phenotype as values and correctly sorted sample IDs
        as index.

    permut : bool
        If true, permute the expression values.

    """
    if permut:
        p = random.sample(range(phenos.shape[0]), phenos.shape[0])
        pd.Series(phenos.values[p], index=phenos.index).to_csv(
            phenof, header=False, sep='\t')
    else:
        phenos.to_csv(phenof, header=False, sep='\t')

def _reml(eigf, remlf, phef, indf, kinf, pEmmax_path, covf=None):
    """
    Make reml file for EMMAX.
    
    Parameters
    ----------
    eigf : str
        Path for output eigR file.

    remlf : str
        Path for output reml file.

    """
    c = ('{} reml --phenof {} --kinf {} --indf {} --out-eigf {} '
         '--out-remlf {}'.format(
             pEmmax_path,
             phef,
             kinf,
             indf,
             eigf,
             remlf))
    if covf:
        c += ' --covf {}'.format(covf)
    subprocess.check_call(c, shell=True)

def run_emmax(
    gene_id, 
    vcfs, 
    regions, 
    phenotypes, 
    samples, 
    kin, 
    tempdir, 
    outdir,
    covariates=None, 
    num_lesser_pval=15, 
    min_permut=1000, 
    max_permut=10000,
    maf=0.05,
    bcftools_path='bcftools',
    pEmmax_path='pEmmax',
    verbose=False,
):
    """
    Run EMMAX for a single gene given a ped file and permuted ped files.
    
    This function will run EMMAX and...
    
    Parameters
    ----------
    gene_id : str
        Gencode gene ID for gene to test.
        
    regions : list
        List of strings of the form 'chr1:100-200'. Biallelic variants in these
        regions will be tested.

    phenotypes : str
        Path to file with phenotype information. Columns should be individual
        IDs that are in the same order as the ind and cov files. Index should be
        phenotype names (i.e. gene IDs).
    
    samples : str
        Path to file defining samples (one per line) to use for EMMAX. 
        
    kin : str
        Path to kinship matrix file.
        
    tempdir : str
        Path to directory where temp directory should be made.
        
    outdir : str
        Path to directory where results will be saved.
        
    covariates : str 
        Covariates file for use with EMMAX. This should be a tab separated file
        with NO header. The first column is sample names. Subsequent columns are
        covariates. Categorical covariates should be binarized.  I think all
        covariates should be numerical (so male-female should be coded as 1-2
        etc.).
        
    num_lesser_pval : int
        Number of p-values less than minimum gene p-value that must be observed
        to stop testing permutations.
        
    min_permut : int
        The minimum number of permutations to test.
        
    max_permut : int
        The maximum number of permutations to test.

    maf : float
        Minor allele frequency cut-off. Variants below this frequency will not
        be tested for association.

    """
    if verbose:
        import datetime
    
    tempdir = os.path.join(tempdir, gene_id)
    try:
        os.makedirs(tempdir)
    except OSError:
        pass

    curdir = os.path.realpath(os.curdir)
    os.chdir(tempdir)
    
    # I believe the EMMAX input files should all be sorted by sample in the same
    # order, so I'll define a sample ordering here and use it throughout.
    sorted_samples = pd.read_table(samples, header=None, squeeze=True)
    sorted_samples.sort_values(inplace=True)

    # Make VCF file. This VCF file will only have biallelic variants in the
    # regions of interest.
    vcf = _make_emmax_vcf(vcfs, gene_id, tempdir, regions, sorted_samples, maf,
                          bcftools_path)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('VCF file created at {}.\n'.format(res))
        sys.stdout.flush()

    # Make ind file.
    indf = os.path.join(tempdir, '{}.ind'.format(gene_id))
    _ind(indf, sorted_samples)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('ind file created at {}.\n'.format(res))
        sys.stdout.flush()

    # Make cov file.
    if covariates:
        covf = os.path.join(tempdir, '{}.cov'.format(gene_id))
        _cov(covf, covariates, sorted_samples)
        if verbose:
            res = str(datetime.datetime.now())
            sys.stdout.write('cov file created at {}.\n'.format(res))
            sys.stdout.flush()
    else:
        covf = None

    # Make phe file.
    phenos = pd.read_table(phenotypes, index_col=0)
    phenos = phenos.ix[gene_id, sorted_samples.values]
    phenof = os.path.join(tempdir, '{}.phe'.format(gene_id))
    _phe(phenof, phenos, permut=False)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('phe file created at {}.\n'.format(res))
        sys.stdout.flush()
    
    # Make reml file.
    eigf = os.path.join(tempdir, '{}.eigR'.format(gene_id))
    remlf = os.path.join(outdir, '{}.reml'.format(gene_id))
    _reml(eigf, remlf, phenof, indf, kin, pEmmax_path, covf=covf)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('reml file created at {}.\n'.format(res))
        sys.stdout.flush()

    # Run association.
    real_res = _emmax(vcf, phenof, indf, eigf, remlf, maf, pEmmax_path)
    out = os.path.join(outdir, '{}.tsv'.format(gene_id))
    real_res.to_csv(out, sep='\t', index=False, na_rep='NA')
    _emmax_cleanup(gene_id)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('First association completed at {}.\n'.format(res))
        sys.stdout.flush()
    min_pval = real_res.PVALUE.min()

    # Do the previous steps for permutations.
    pvalues = []
    min_pvalues = []
    reml_info = []
    i = 0
    num_lesser_pvals = 0
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('Permutations started at {}.\n'.format(res))
        sys.stdout.flush()
    while i < max_permut:
        # Run EMMAX for this permutation.
        phenof = os.path.join(tempdir, '{}.phe'.format(gene_id))
        _phe(phenof, phenos, permut=True)
        eigf = os.path.join(tempdir, '{}.eigR'.format(gene_id))
        remlf = os.path.join(tempdir, '{}.reml'.format(gene_id))
        _reml(eigf, remlf, phenof, indf, kin, pEmmax_path, covf=covf)
        out = os.path.join(tempdir, '{}_permutation_{}.tsv'.format(
            gene_id, str(i + 1).zfill(len(str(max_permut)))))
        res = _emmax(vcf, phenof, indf, eigf, remlf, maf, pEmmax_path)

        # Read results.
        res.index = ('chr' + res['#CHROM'].astype(str) + ':' + 
                     res.BEG.astype(str))
        pvalues.append(res.PVALUE.values)
        m = res.PVALUE.min()
        min_pvalues.append(m)
        reml_info.append(pd.read_table(remlf, header=None, index_col=0,
                                       squeeze=True))
        if m < min_pval:
            num_lesser_pvals += 1
        if (num_lesser_pvals >= num_lesser_pval) and (i + 1 >= min_permut):
            break
        i += 1
        
        _emmax_cleanup(gene_id)
        if verbose:
            res = str(datetime.datetime.now())
            sys.stdout.write('Permutation {} completed at {}.\n'.format(i, res))
        sys.stdout.flush()

    # Remove VCF file.
    os.remove('{}.vcf.gz'.format(gene_id))
    os.remove('{}.vcf.gz.tbi'.format(gene_id))
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('VCF file removed at {}.\n'.format(res))
        sys.stdout.flush()
   
    pvalues = pd.DataFrame(pvalues).T
    pvalues.to_csv(os.path.join(outdir, 'permuted_pvalues.tsv'), index=None, 
                   sep='\t', header=None)
    min_pvalues = pd.Series(min_pvalues, index=None)
    min_pvalues.to_csv(os.path.join(outdir, 'minimum_pvalues.tsv'), sep='\t',
                       index=None)
    reml_info = pd.DataFrame(reml_info)
    reml_info.to_csv(os.path.join(outdir, 'permuted_reml.tsv'), sep='\t',
                     index=None)
    
    shutil.rmtree(tempdir)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('Finished everything at {}.\n'.format(res))
        sys.stdout.flush()

def _emmax(vcf, phe, ind, eig, reml, maf, pEmmax_path):
    """
    Execute EMMAX command.

    Parameters
    ----------
    vcf : str
        Path to VCF file that contains only the variants to test.
        
    """
    c = ('{} assoc --vcf {} --phenof {} --field GT --indf {} --eigf {} '
         '--remlf {} --out-assocf /dev/stdout --minMAF {} --maxMAF 1 --maxMAC '
         '1000000000 --minRSQ 0 --minCallRate 0.5 --minMAC 3'.format(
             pEmmax_path,
             vcf,
             phe,
             ind,
             eig,
             reml,
             maf,
         ))
    res = subprocess.check_output(c, shell=True)
    lines = res.split('\n')[0:-1]
    lines = [x.split('\t') for x in lines]
    df = pd.DataFrame(lines[1:], columns=lines[0])
    df['PVALUE'] = pd.to_numeric(df['PVALUE'], errors='coerce')
    return df

def _emmax_cleanup(prefix):
    """
    Delete extra EMMAX files that we don't need.
    """
    to_delete = ['eigR', 'phe', 'reml']
    for suffix in to_delete:
        try:
            os.remove('{}.{}'.format(prefix, suffix))
        except OSError:
            continue

def _make_emmax_vcf(
    vcfs, 
    gene_id, 
    tempdir, 
    regions, 
    samples, 
    maf,
    bcftools_path,
):
    """
    Make VCF file for EMMAX. Rather than reading into the larger VCF file, we'll
    make a small VCF with just the samples and variants to test for this gene.
    
    Parameters
    ----------
    vcf : str
        List of paths to VCF files with samples and variants. Must be compressed
        and indexed.

    gene_id : str
        Gene ID used for naming files.

    regions : list
        List of regions to get variants from.

    samples : str
        Pandas Series with samples (in correctly sorted order) as values.

    maf : float
        Minor allele frequency cut off. Any variants with MAF lower than this
        will be filtered out.

    """
    if len(vcfs) == 1:
        fn = os.path.join(tempdir, '{}.vcf.gz'.format(gene_id))
        c = ('{} view {} -q {}:minor -m2 -M2 -r {} -s {} -O u | '
             '{} filter -m x -O z > {}'.format(
                 bcftools_path,
                 vcfs[0],
                 maf,
                 ','.join(regions),
                 ','.join(samples.values),
                 bcftools_path,
                 fn))
        subprocess.check_call(c, shell=True)
        c = ('{} index --tbi {}'.format(bcftools_path, fn))
        subprocess.check_call(c, shell=True)
    else:
        temp_vcfs = []
        for i,v in enumerate(vcfs):
            fn = os.path.join(tempdir, '{}.vcf.gz'.format(i))
            temp_vcfs.append(fn)
            c = ('{} view {} -q {}:minor -m2 -M2 -r {} -s {} -O u | '
                 '{} filter -m x -O z > {}'.format(
                     bcftools_path,
                     v,
                     maf,
                     ','.join(regions),
                     ','.join(samples.values),
                     bcftools_path,
                     fn))
            subprocess.check_call(c, shell=True)
            c = ('{} index --tbi {}'.format(bcftools_path, fn))
            subprocess.check_call(c, shell=True)
        fn = os.path.join(tempdir, '{}.vcf.gz'.format(gene_id))
        c = '{} concat -Oz -a {} > {}'.format(
            bcftools_path, ' '.join(temp_vcfs), fn)
        subprocess.check_call(c, shell=True)
        c = ('{} index --tbi {}'.format(bcftools_path, fn))
        subprocess.check_call(c, shell=True)
        for v in temp_vcfs:
            os.remove(v)
            os.remove('{}.tbi'.format(v))

    return fn
    
def main():
    import glob

    num_lesser_pval = 15
    min_permut = 1000
    max_permut = 10000
    maf = 0.05
    tempdir = '/dev/shm'
    bcftools_path = 'bcftools'
    pEmmax_path = 'pEmmax'
    parser = argparse.ArgumentParser(description=(
        'This script runs EMMAX for a single gene. EMMAX is run for the "real" '
        'data and permuted data. The testing procedure is based on the one '
        'from the 2015 GTEx project paper (10.1126/science.1262110).'))
    parser.add_argument('gene_id', help=('Gene ID for gene to test. This '
                                         'should be a row in the phenotypes '
                                         'file.'))
    parser.add_argument('vcfs', help=('Compressed, indexed VCF files with '
                                      'variants to test. VCFs in regions will'
                                      'be extracted and concatenated. Files '
                                      'should be separated by commas.'))
    parser.add_argument('regions', help=(
        'List of regions of the form chr3:100-200 or (3:100-200 depending on '
        'your VCF. Multiple regions are '
        'separated by commas: chr3:100-200,chr10:400-500. Biallelic variants '
        'in these regions will be tested.'))
    parser.add_argument('phenotypes', 
                        help=('TSV file with gene expression (or other) '
                              'values. Index '
                              'should be gene IDs and columns should be '
                              'sample IDs that match the VCF file.'))
    parser.add_argument('samples', help=(
        'File with samples (one per line) to use from VCF file. The sample IDs '
        'should match those in the VCF file.'))
    parser.add_argument('kin', help=('Kinship matrix file (.kinf).'))
    parser.add_argument('outdir', help=('Directory to store final results.'))
    parser.add_argument('-c', metavar='covariates', help=(
        'Covariates file for use with EMMAX. This should be a tab separated '
        'file with NO header. The first column is sample names. Subsequent '
        'columns are covariates. Categorical covariates should be binarized. '
        'I think all covariates should be numerical (so male-female should be '
        'coded as 1-2 etc.).'))
    parser.add_argument('-n', metavar='num_lesser_pval',
                        default=num_lesser_pval, type=int, help=(
                            'Minimum number of gene-level permutation p-values '
                            'less than "real" p-value required to stop testing '
                            'permutations. Default: {:,}.'.format(
                                num_lesser_pval)))
    parser.add_argument('-i', metavar='min_permut',
                        default=min_permut, type=int, help=(
                            'Minimum number of permutations to perform. '
                            'Default: {:,}.'.format(min_permut)))
    parser.add_argument('-a', metavar='max_permut',
                        default=max_permut, type=int, help=(
                            'Maximum number of permutations to perform. '
                            'Default: {:,}.'.format(max_permut)))
    parser.add_argument('-m', metavar='maf',
                        default=maf, type=float, help=(
                            'Minor allele frequency cut-off. Variants below '
                            'this frequency will not be tested. Default: '
                            '{}.'.format(maf)))
    parser.add_argument('-t', metavar='tempdir', help=(
        'Temporary directory. Default: {}.'.format(tempdir)),
                        default=tempdir)
    parser.add_argument('-b', metavar='bcftools_path', help=(
        'Path to bcftools. Default: {}.'.format(bcftools_path)),
                        default=bcftools_path)
    parser.add_argument('-e', metavar='pEmmax_path', help=(
        'Path to pEmmax. Default: {}.'.format(pEmmax_path)),
        default=pEmmax_path)
    parser.add_argument('--verbose', help=(
        'Print log information to stdout.'), action='store_true')
    args = parser.parse_args()
    gene_id = args.gene_id
    vcfs = args.vcfs.split(',')
    regions = args.regions.split(',')
    phenotypes = args.phenotypes
    samples = args.samples
    kin = args.kin
    outdir = args.outdir
    covariates = args.c
    num_lesser_pval = args.n
    min_permut = args.i
    max_permut = args.a
    maf = args.m
    tempdir = args.t
    bcftools_path = args.b
    pEmmax_path = args.e
    verbose = args.verbose

    run_emmax(
        gene_id, 
        vcfs, 
        regions, 
        phenotypes, 
        samples, 
        kin, 
        tempdir, 
        outdir,
        covariates=covariates,
        num_lesser_pval=num_lesser_pval,
        min_permut=min_permut, 
        max_permut=max_permut,
        maf=maf,
        bcftools_path=bcftools_path,
        pEmmax_path=pEmmax_path,
        verbose=verbose,
    )

if __name__ == '__main__':
    main()
