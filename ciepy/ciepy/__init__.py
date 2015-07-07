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
