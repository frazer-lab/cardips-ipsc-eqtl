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
        
def read_emmax_output(fn):
    with gzip.open(fn) as f:
        lines = [x.strip().split('\t') for x in f.readlines()]
    lines[0][0] = lines[0][0][1:]
    res = pd.DataFrame(lines[1:], columns=lines[0])
    res = res.convert_objects(convert_numeric=True)
    return res
