import argparse
import gzip
import os
import shutil
import subprocess

import pandas as pd
import projectpy as ppy

def run_emmax(gene_id, vcf, regions, ped, kinship_matrix, tempdir, outdir,
              permuted_peds, covariates=[], num_lesser_pval=15, min_permut=1000,
              max_permut=10000):
    """
    Run EMMAX for a single gene given a ped file and permuted ped files.
    
    This function will run EMMAX and 
    
    Parameters
    ----------
    gene_id : str
        Gencode gene ID for gene to test.
        
    regions : list
        List of strings of the form 'chr1:100-200'. Biallelic SNVs in these
        regions will be tested.
    
    ped : str
        Path to PED file with actual genotypes (i.e. not permuted).
        
    kinship_matrix : str
        Path to kinship matrix file.
        
    tempdir : str
        Path to directory where temp directory should be made.
        
    outdir : str
        Path to directory where results will be saved.
        
    permuted_peds : list
        List of strings of paths to permuted PED files.

    covariates : list
        List of covariates that correspond to columns in PED file. These will be
        used when running EMMAX.
        
    num_lesser_pval : int
        Number of p-values less than minimum gene p-value that must be observed
        to stop testing permutations.
        
    min_permut : int
        The minimum number of permutations to test.
        
    max_permut : int
        The maximum number of permutations to test.

    """
    tempdir = os.path.join(tempdir, gene_id)
    ppy.makedir(tempdir)

    curdir = os.path.realpath(os.curdir)
    os.chdir(tempdir)
    
    # Make VCF file.
    vcf = _make_emmax_vcf(vcf, gene_id, tempdir, regions)
    
    # Run EMMAX for real data.
    _emmax(gene_id, ped, kinship_matrix, vcf, gene_id, covariates=covariates)
    out = '{}.epacts.gz'.format(gene_id)
    real_res = read_emmax_output(out)
    min_pval = real_res.PVALUE.min()
    
    # Run EMMAX for permuted data.
    ped_df = pd.read_table(ped)
    names = []
    pvalues = []
    min_pvalues = []
    i = 0
    num_lesser_pvals = 0
    while i < len(permuted_peds) and i < max_permut:
        fn = permuted_peds[i]
        prefix = os.path.splitext(os.path.split(fn)[1])[0]
        names.append(prefix)
        tdf = pd.read_table(fn)
        tdf[gene_id] = ped_df[gene_id].values
        ped_fn = os.path.join(tempdir, '{}.ped'.format(prefix))
        tdf.to_csv(ped_fn, index=False, sep='\t')
        _emmax(gene_id, ped_fn, kinship_matrix, vcf, prefix)
        out = '{}.epacts.gz'.format(prefix)
        res = read_emmax_output(out)
        res.index = 'chr' + res.CHROM.astype(str) + ':' + res.BEG.astype(str)
        pvalues.append(res.PVALUE)
        m = res.PVALUE.min()
        min_pvalues.append(m)
        if m < min_pval:
            num_lesser_pvals += 1
        res.PVALUE.to_csv('{}_pvalues.tsv'.format(prefix), sep='\t')
        c = 'rm {0}.epacts.gz {0}.ped'.format(prefix)
        subprocess.check_call(c, shell=True)
        if (num_lesser_pvals >= num_lesser_pval) and (i + 1 >= min_permut):
            break
        i += 1
        if i % 50 == 0:
            print('Finished {}'.format(i))
        
    pvalues = pd.DataFrame(pvalues, index=names).T
    pvalues.to_csv('permuted_pvalues.tsv', sep='\t')
    min_pvalues = pd.Series(min_pvalues, index=names)
    min_pvalues.to_csv('minimum_pvalues.tsv', sep='\t')
    
    # Remove VCF file.
    c = 'rm {0}.vcf.gz {0}.vcf.gz.tbi'.format(gene_id)
    subprocess.check_call(c, shell=True)
    
    # Copy output to outdir.
    outdir = os.path.join(outdir, gene_id)
    ppy.makedir(outdir)
    shutil.move('{}.epacts.gz'.format(gene_id), outdir)
    shutil.move('permuted_pvalues.tsv', outdir)
    shutil.move('minimum_pvalues.tsv', outdir)
    shutil.rmtree(tempdir)
    
    os.chdir(curdir)

def read_emmax_output(fn):
    """
    Read gzipped EMMAX output file and return as dataframe.
    """
    with gzip.open(fn) as f:
        lines = [x.strip().split('\t') for x in f.readlines()]
    lines[0][0] = lines[0][0][1:]
    res = pd.DataFrame(lines[1:], columns=lines[0])
    res = res.convert_objects(convert_numeric=True)
    return res
    
def _emmax(gene_id, ped, kinship_matrix, vcf, prefix, covariates=[]):
    """
    Execute EMMAX command.

    Parameters
    ----------
    gene_id : str
        Gencode gene ID for gene to test.
    
    ped : str
        Path to PED file with actual genotypes (i.e. not permuted).
        
    kinship_matrix : str
        Path to kinship matrix file.
        
    vcf : str
        Path to VCF file with SNVs to test.
        
    prefix : str
        Prefix for naming output files.

    covariates : list
        List of covariates that correspond to columns in PED file. These will be
        used when running EMMAX.
        
    """
    covs = ' '.join(['--cov {}'.format(c) for c in covariates])
    c = ('{} single --vcf {} --ped {} --min-maf 0.1 --kin {} --pheno {} '
         '{} --test q.emmax --out {} --run 4'.format(
        ppy.epacts,
        vcf,
        ped,
        kinship_matrix,
        gene_id,
        covs,
        prefix))
    subprocess.check_call(c, shell=True)
    _delete_extra_files(prefix)
    
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
        'This script runs EMMAX using EPACTS for a single gene. EMMAX is run '
        'for the "real" data and permuted data. Permuted PED files with all '
        'covariates except the gene expression are required as input. The '
        'testing procedure is based on the one from the 2015 GTEx project '
        'paper (10.1126/science.1262110).'))
    parser.add_argument('gene_id', help=('Gene ID for gene to test. This '
                                         'should be a column in the PED file.'))
    parser.add_argument('vcf', help=('VCF file with variants.'))
    parser.add_argument('regions', help=(
        'List of regions of the form chr3:100-200. Multiple regions are '
        'separated by commas: chr3:100-200,chr10:400-500. Biallelic SNVs in '
        'these regions will be tested.'))
    parser.add_argument('ped', help=('PED file with covariates and gene '
                                     'expression values.'))
    parser.add_argument('kinf', help=('Kinship matrix file (.kinf).'))
    parser.add_argument('outdir', help=('Directory to store final results.'))
    parser.add_argument('permuted_peds', help=(
        'Directory with permuted ped files. All files with .ped suffix are '
        'assumed to be permuted ped files. These files will be sorted by '
        'filename and used in that order.'))
    parser.add_argument('-c', metavar='covariates', help=(
        'Covariates to include when running EMAX. Multiple covariates should '
        'be separated by commas: SEX,AGE'), default=[])
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
    ped = args.ped
    kinf = args.kinf
    outdir = args.outdir
    permuted_peds = os.path.realpath(args.permuted_peds)
    covariates = args.c.split(',')
    num_lesser_pval = args.n
    min_permut = args.i
    max_permut = args.a
    tempdir = args.t

    permuted_peds = glob.glob(os.path.join(permuted_peds, '*.ped'))
    assert len(permuted_peds) >= min_permut, \
            ('Not enough permuted PED files to perform minimum number of '
             'permutations.')
    assert len(permuted_peds) >= max_permut, \
            ('Not enough permuted PED files to perform maximum number of '
             'permutations.')

    run_emmax(gene_id, vcf, regions, ped, kinf, tempdir, outdir, permuted_peds,
              covariates=covariates, num_lesser_pval=num_lesser_pval,
              min_permut=min_permut, max_permut=max_permut)

if __name__ == '__main__':
    main()
