import argparse
import gzip
import os
import shutil
import subprocess
import sys

import pandas as pd
import projectpy as ppy

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

def _reml(eigf, remlf, phe, ind, kin, cov=None):
    """
    Make reml file for EMMAX.
    
    Parameters
    ----------
    eigf : str
        Path for output eigR file.

    remlf : str
        Path for output reml file.

    """
    c = ('{}/pEmmax reml --phenof {} --kinf {} --indf {} --out-eigf {} '
         '--out-remlf {}'.format(
             os.path.split(ppy.epacts)[0],
             phe,
             kin,
             ind,
             eigf,
             remlf))
    if cov:
        c += ' --covf {}'.format(cov)
    subprocess.check_call(c, shell=True)

def run_emmax(gene_id, vcf, regions, phenotypes, ind, kin, tempdir, outdir,
              cov=None, num_lesser_pval=15, min_permut=1000, max_permut=10000):
    """
    Run EMMAX for a single gene given a ped file and permuted ped files.
    
    This function will run EMMAX and...
    
    Parameters
    ----------
    gene_id : str
        Gencode gene ID for gene to test.
        
    regions : list
        List of strings of the form 'chr1:100-200'. Biallelic SNVs in these
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

    tempdir = os.path.join(tempdir, gene_id)
    ppy.makedir(tempdir)

    curdir = os.path.realpath(os.curdir)
    os.chdir(tempdir)

    # Make VCF file. This VCF file will only have biallelic SNVs in the regions
    # of interest.
    vcf = _make_emmax_vcf(vcf, gene_id, tempdir, regions)

    # Make phe file.
    phenos = pd.read_table(phenotypes, index_col=0)
    phenof = os.path.join(tempdir, '{}.phe'.format(gene_id))
    _phe(phenof, gene_id, phenos)
    
    # Make reml file.
    eigf = os.path.join(tempdir, '{}.eigR'.format(gene_id))
    remlf = os.path.join(tempdir, '{}.reml'.format(gene_id))
    _reml(eigf, remlf, phenof, ind, kin, cov=cov)

    # Run association.
    out = os.path.join(outdir, '{}.tsv'.format(gene_id))
    _emmax(out, vcf, phenof, ind, eigf, remlf)
    _emmax_cleanup(gene_id)
    real_res = pd.read_table(out)
    min_pval = real_res.PVALUE.min()

    # Do the previous steps for permutations.
    pvalues = []
    min_pvalues = []
    i = 0
    num_lesser_pvals = 0
    while i < max_permut:
        # Run EMMAX for this permutation.
        p = random.sample(range(phenos.shape[1]), phenos.shape[1])
        phenof = os.path.join(tempdir, '{}.phe'.format(gene_id))
        _phe(phenof, gene_id, phenos, permut=p)
        eigf = os.path.join(tempdir, '{}.eigR'.format(gene_id))
        remlf = os.path.join(tempdir, '{}.reml'.format(gene_id))
        _reml(eigf, remlf, phenof, ind, kin, cov=cov)
        out = os.path.join(tempdir, '{}_permutation_{}.tsv'.format(
            gene_id, str(i + 1).zfill(len(str(max_permut)))))
        _emmax(out, vcf, phenof, ind, eigf, remlf)
        _emmax_cleanup(gene_id)

        # Read results.
        res = pd.read_table(out)
        res.index = ('chr' + res['#CHROM'].astype(str) + ':' + 
                     res.BEG.astype(str))
        pvalues.append(res.PVALUE)
        m = res.PVALUE.min()
        min_pvalues.append(m)
        if m < min_pval:
            num_lesser_pvals += 1
        if (num_lesser_pvals >= num_lesser_pval) and (i + 1 >= min_permut):
            break
        i += 1
        if i % 50 == 0:
            sys.stderr.write('Finished {} permutations'.format(i))
        os.remove(out)
        

    # Remove VCF file.
    os.remove('{}.vcf.gz'.format(gene_id))
    os.remove('{}.vcf.gz.tbi'.format(gene_id))
   
    pvalues = pd.DataFrame(pvalues).T
    pvalues.to_csv(os.path.join(outdir, 'permuted_pvalues.tsv'), index=None, 
                   sep='\t', header=None)
    min_pvalues = pd.Series(min_pvalues, index=None)
    min_pvalues.to_csv(os.path.join(outdir, 'minimum_pvalues.tsv'), sep='\t',
                       index=None)
    
    shutil.rmtree(tempdir)

def _emmax(out, vcf, phe, ind, eig, reml):
    """
    Execute EMMAX command.

    Parameters
    ----------
    out : str
        Output file?
    
    vcf : str
        Path to VCF file that contains only the SNVs to test.
        
    """
    c = ('{}/pEmmax assoc --vcf {} --phenof {} --field GT --indf {} --eigf {} '
         '--remlf {} --out-assocf {} --minMAF 0.1 --maxMAF 1 --maxMAC '
         '1000000000 --minRSQ 0 --minCallRate 0.5 --minMAC 3'.format(
             os.path.split(ppy.epacts)[0],
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

def _make_emmax_vcf(vcf, gene_id, tempdir, regions):
    import ciepy as cpy

    fn = os.path.join(tempdir, '{}.vcf.gz'.format(gene_id))
    c = ('{} view {} -q 0.1:minor -m2 -M2 -v snps -r {} | '
         '{} annotate --rename-chrs {} -O z > {}'.format(
             ppy.bcftools,
             vcf,
             ','.join(regions),
             ppy.bcftools,
             os.path.join(cpy.root, 'data', 'chromosome_conversion.txt'),
             fn))
    subprocess.check_call(c, shell=True)

    c = ('{} index --tbi {}'.format(ppy.bcftools, fn))
    subprocess.check_call(c, shell=True)
    return fn
    
def _delete_extra_files(prefix):
    """
    Delete extra EMMAX files that we don't need.
    """
    to_delete = ['cov', 'eigR', 'epacts.conf', 'epacts.gz.tbi', 'epacts.mh.pdf',
                 'epacts.OK', 'epacts.qq.pdf', 'epacts.R', 'epacts.top5000',
                 'ind', 'Makefile', 'phe', 'reml']
    for suffix in to_delete:
        fn = '{}.{}'.format(prefix, suffix)
        if os.path.exists(fn):
            os.remove(fn)

def main():
    import glob

    num_lesser_pval = 15
    min_permut = 1000
    max_permut = 10000
    tempdir = '/dev/shm'
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
        'separated by commas: chr3:100-200,chr10:400-500. Biallelic SNVs in '
        'these regions will be tested.'))
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

    run_emmax(gene_id, vcf, regions, phenotypes, ind, kin, tempdir, outdir,
              cov=cov, num_lesser_pval=num_lesser_pval,
              min_permut=min_permut, max_permut=max_permut)

if __name__ == '__main__':
    main()
