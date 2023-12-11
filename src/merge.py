#!/usr/bin/env python3

import os
import subprocess

def merge_root_files(filein, fileout):
    print('Generating:', fileout)
    try:
        subprocess.run(['hadd', '-k', '-O', '-j', fileout] + filein).check_returncode()
    except:
        subprocess.run(['rm', '-f', fileout])
        raise

def get_mtime(path):
    try:
        return os.stat(path).st_mtime
    except Exception:
        return 0.0

def update_merged_root_file(filein, fileout):
    print('Checking:', fileout)
    mtime_in = max(get_mtime(f) for f in filein)
    mtime_out = get_mtime(fileout)
    if mtime_in < mtime_out: return
    merge_root_files(filein, fileout)
