import gzip 
import os
from os.path import join
import subprocess

import pandas as pd

# File locations
# Rather than specifying paths all over, I'll just put the paths to commonly
# used files and directories here and use the package to open them. 
def _get_project_dir():
    return os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-3])

root = _get_project_dir()

## functions

def makedir(p):
    """Make a directory if it doesn't already exist"""
    try:
        os.makedirs(p)
    except OSError:
        pass

def submit_job(fn):
    subprocess.check_call('ssh cdeboever@flc.ucsd.edu \'qsub {}\''.format(fn),
                          shell=True)
        
#def read_emmax_output(fn):
#    with gzip.open(fn) as f:
#        lines = [x.strip().split('\t') for x in f.readlines()]
#    lines[0][0] = lines[0][0][1:]
#    res = pd.DataFrame(lines[1:], columns=lines[0])
#    res = res.convert_objects(convert_numeric=True)
#    return res

def read_emmax_output(fn):
    res = pd.read_table(fn, index_col=None)
    res.columns = [c.replace('#', '') for c in res.columns]
    return res

def calc_bed_enrichment_from_url(url, variants, variants_window):
    """Calculate enrichment for bed file from a URL for variants
    vs. variants_window"""
    bt = pbt.BedTool(cpb.general.read_gzipped_text_url(url), from_string=True)
    bt = bt.sort()
    bt = bt.merge()
    res = variants.intersect(bt, sorted=True, wo=True)
    eqtl_in_peak = len(res)
    eqtl_out_peak = len(variants) - eqtl_in_peak

    res = variants_window.intersect(bt, sorted=True, wo=True)
    not_eqtl_in_peak = 0
    for r in res:
        not_eqtl_in_peak += int(r.fields[-1])
    not_eqtl_in_peak -= eqtl_in_peak
    
    total = 0
    for r in variants_window:
        total += r.length
    not_eqtl_out_peak = total - not_eqtl_in_peak - eqtl_in_peak - eqtl_out_peak
    
    oddsratio, p = scipy.stats.fisher_exact([[eqtl_in_peak, eqtl_out_peak],
                                             [not_eqtl_in_peak, not_eqtl_out_peak]])
    return url, oddsratio, p

def clean_axis(ax):
    "Remove spines and ticks from axis"
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def _comma_func(x, pos):
    """
    Formatter function takes tick label and tick position.
    """
    s = '{:0,d}'.format(int(x))
    return s

# Use: ax.yaxis.set_major_formatter(ds.comma_format)
import matplotlib.ticker as tkr
comma_format = tkr.FuncFormatter(_comma_func)
