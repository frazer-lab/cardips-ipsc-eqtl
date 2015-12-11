import argparse
import gzip
import os
import shutil
import subprocess
import sys

import pandas as pd

def _ind(indf, ind):
    """
    Make .ind file for EMMAX. The first column is individual ID and the second
    column is the position of the sample in VCF file. Since I make the VCF based
    on the input ind file, I just use the order of the samples to write the temp
    ind file.
    
    Parameters
    ----------
    inf : str
        Path for output ind file.

    ind : str
        Input ind file.

    phenos : pandas.DataFrame
        DataFrame with phenotype values.

    permut : list
        A list of indices to permute the phenotype data before writing it.

    """
    se = pd.read_table(ind, index_col=0, squeeze=True, header=None)
    se = pd.Series(range(1, se.shape[0] + 1), index=se.index)
    se.to_csv(indf, header=False, sep='\t')

def _phe(phenof, gene_id, phenos, permut=None):
    """
    Make .phe file for EMMAX. The first column is individual ID and the second
    column is the phenotype value.
    
    Parameters
    ----------
    phenof : str
        Path for output phe file.

    gene_id : str
        ID for gene (must be contained in phenos.index).

    phenos : pandas.DataFrame
        DataFrame with phenotype values.

    permut : list
        A list of indices to permute the phenotype data before writing it.

    """
    se = phenos.ix[gene_id]
    if permut:
        se = pd.Series(se.values[permut], index=se.index)
    se.to_csv(phenof, header=False, sep='\t')

def _reml(eigf, remlf, phe, ind, kin, pEmmax_path, cov=None):
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
             phe,
             kin,
             ind,
             eigf,
             remlf))
    if cov:
        c += ' --covf {}'.format(cov)
    subprocess.check_call(c, shell=True)

def run_emmax(
    gene_id, 
    vcf, 
    regions, 
    phenotypes, 
    ind, 
    kin, 
    tempdir, 
    outdir,
    cov=None, 
    num_lesser_pval=15, 
    min_permut=1000, 
    max_permut=10000,
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
    
    ind : str
        Path to ind file for EMMAX. Order should match order of cov phenotypes
        files.
        
    kin : str
        Path to kinship matrix file.
        
    tempdir : str
        Path to directory where temp directory should be made.
        
    outdir : str
        Path to directory where results will be saved.
        
    cov : str 
        Path to covariates file for EMMAX.
        
    num_lesser_pval : int
        Number of p-values less than minimum gene p-value that must be observed
        to stop testing permutations.
        
    min_permut : int
        The minimum number of permutations to test.
        
    max_permut : int
        The maximum number of permutations to test.

    """
    import random
    random.seed(20150605)
    if verbose:
        import datetime

    tempdir = os.path.join(tempdir, gene_id)
    try:
        os.makedirs(tempdir)
    except OSError:
        pass

    curdir = os.path.realpath(os.curdir)
    os.chdir(tempdir)

    # Make VCF file. This VCF file will only have biallelic variants in the
    # regions of interest.
    vcf = _make_emmax_vcf(vcf, gene_id, tempdir, regions, ind, bcftools_path)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('VCF file created at {}.\n'.format(res))
        sys.stdout.flush()

    # Make ind file.
    indf = os.path.join(tempdir, '{}.ind'.format(gene_id))
    _ind(indf, ind)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('ind file created at {}.\n'.format(res))
        sys.stdout.flush()

    # Make phe file.
    order = pd.read_table(ind, index_col=0, header=None)
    phenos = pd.read_table(phenotypes, index_col=0)[order.index]
    phenof = os.path.join(tempdir, '{}.phe'.format(gene_id))
    _phe(phenof, gene_id, phenos)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('phe file created at {}.\n'.format(res))
        sys.stdout.flush()
    
    # Make reml file.
    eigf = os.path.join(tempdir, '{}.eigR'.format(gene_id))
    remlf = os.path.join(outdir, '{}.reml'.format(gene_id))
    _reml(eigf, remlf, phenof, indf, kin, pEmmax_path, cov=cov)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('reml file created at {}.\n'.format(res))
        sys.stdout.flush()

    # Run association.
    out = os.path.join(outdir, '{}.tsv'.format(gene_id))
    _emmax(out, vcf, phenof, indf, eigf, remlf, pEmmax_path)
    _emmax_cleanup(gene_id)
    if verbose:
        res = str(datetime.datetime.now())
        sys.stdout.write('First association completed at {}.\n'.format(res))
        sys.stdout.flush()
    real_res = pd.read_table(out)
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
        p = random.sample(range(phenos.shape[1]), phenos.shape[1])
        phenof = os.path.join(tempdir, '{}.phe'.format(gene_id))
        _phe(phenof, gene_id, phenos, permut=p)
        eigf = os.path.join(tempdir, '{}.eigR'.format(gene_id))
        remlf = os.path.join(tempdir, '{}.reml'.format(gene_id))
        _reml(eigf, remlf, phenof, indf, kin, pEmmax_path, cov=cov)
        out = os.path.join(tempdir, '{}_permutation_{}.tsv'.format(
            gene_id, str(i + 1).zfill(len(str(max_permut)))))
        _emmax(out, vcf, phenof, indf, eigf, remlf, pEmmax_path)

        # Read results.
        res = pd.read_table(out)
        res.index = ('chr' + res['#CHROM'].astype(str) + ':' + 
                     res.BEG.astype(str))
        pvalues.append(res.PVALUE)
        m = res.PVALUE.min()
        min_pvalues.append(m)
        reml_info.append(pd.read_table(remlf, header=None, index_col=0,
                                       squeeze=True))
        if m < min_pval:
            num_lesser_pvals += 1
        if (num_lesser_pvals >= num_lesser_pval) and (i + 1 >= min_permut):
            break
        i += 1
        os.remove(out)
        
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
        sys.stdout.write('VCF file removed at {}.\n'.format(i, res))
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
        sys.stdout.write('Finished at {}.\n'.format(i, res))
        sys.stdout.flush()

def _emmax(out, vcf, phe, ind, eig, reml, pEmmax_path):
    """
    Execute EMMAX command.

    Parameters
    ----------
    out : str
        Output file?
    
    vcf : str
        Path to VCF file that contains only the variants to test.
        
    """
    c = ('{} assoc --vcf {} --phenof {} --field GT --indf {} --eigf {} '
         '--remlf {} --out-assocf {} --minMAF 0.1 --maxMAF 1 --maxMAC '
         '1000000000 --minRSQ 0 --minCallRate 0.5 --minMAC 3'.format(
             pEmmax_path,
             vcf,
             phe,
             ind,
             eig,
             reml,
             out
         ))
    subprocess.check_call(c, shell=True)

def _emmax_cleanup(prefix):
    """
    Delete extra EMMAX files that we don't need.
    """
    to_delete = ['eigR', 'phe', 'reml']
    for suffix in to_delete:
        fn = '{}.{}'.format(prefix, suffix)
        if os.path.exists(fn):
            os.remove(fn)

def _make_emmax_vcf(vcf, gene_id, tempdir, regions, ind, bcftools_path):
    import ciepy

    ind = pd.read_table(ind, index_col=0, header=None)
    fn = os.path.join(tempdir, '{}.vcf.gz'.format(gene_id))
    c = ('{} view {} -q 0.05:minor -m2 -M2 -r {} -s {} -O u | '
         '{} filter -m x -O z > {}'.format(
             bcftools_path,
             vcf,
             ','.join(regions),
             ','.join(ind.index),
             bcftools_path,
             fn))
    subprocess.check_call(c, shell=True)

    c = ('{} index --tbi {}'.format(bcftools_path, fn))
    subprocess.check_call(c, shell=True)
    return fn
    
def main():
    import glob

    num_lesser_pval = 15
    min_permut = 1000
    max_permut = 10000
    tempdir = '/dev/shm'
    bcftools_path = 'bcftools'
    pEmmax_path = 'pEmmax'
    parser = argparse.ArgumentParser(description=(
        'This script runs EMMAX for a single gene. EMMAX is run for the "real" '
        'data and permuted data. The testing procedure is based on the one '
        'from the 2015 GTEx project paper (10.1126/science.1262110).'))
    parser.add_argument('gene_id', help=('Gene ID for gene to test. This '
                                         'should be a column in the phenotypes '
                                         'file.'))
    parser.add_argument('vcf', help=('VCF file with variants.'))
    parser.add_argument('regions', help=(
        'List of regions of the form chr3:100-200. Multiple regions are '
        'separated by commas: chr3:100-200,chr10:400-500. Biallelic variants '
        'in these regions will be tested.'))
    parser.add_argument('phenotypes', 
                        help=('TSV file with gene expression values. Index '
                              'should be gene IDs and columns should be '
                              'sample IDs that match the ind file.'))
    parser.add_argument('ind', help=('ind file with positions of each sample '
                                     'in the VCF header.'))
    parser.add_argument('kin', help=('Kinship matrix file (.kinf).'))
    parser.add_argument('outdir', help=('Directory to store final results.'))
    parser.add_argument('-c', metavar='cov', 
                        help=('Covariates file for EMMAX.'))
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
    vcf = args.vcf
    regions = args.regions.split(',')
    phenotypes = args.phenotypes
    ind = args.ind
    kin = args.kin
    outdir = args.outdir
    cov = args.c
    num_lesser_pval = args.n
    min_permut = args.i
    max_permut = args.a
    tempdir = args.t
    bcftools_path = args.b
    pEmmax_path = args.e
    verbose = args.verbose

    run_emmax(
        gene_id, 
        vcf, 
        regions, 
        phenotypes, 
        ind, 
        kin, 
        tempdir, 
        outdir,
        cov=cov, 
        num_lesser_pval=num_lesser_pval,
        min_permut=min_permut, 
        max_permut=max_permut,
        bcftools_path=bcftools_path,
        pEmmax_path=pEmmax_path,
        verbose=verbose,
    )

if __name__ == '__main__':
    main()
