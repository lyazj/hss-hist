#!/usr/bin/env python3

import os
import subprocess

def merge_root_files(filein, fileout):
    print('Generating:', fileout)
    import uproot
    filein = [file for file in filein if not print('Adding', file) and uproot.open(file)['Events'].num_entries]
    try:
        subprocess.run(['rm', '-f', fileout])
        subprocess.run(['mkdir', '-p', os.path.dirname(fileout)])
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
    dirnames = sorted(set(os.path.dirname(f) for f in filein))
    if len(dirnames) != 1: raise RuntimeError('expect exact 1 directory, got %s' % dirnames)
    mtime_in = max(get_mtime(dirnames[0]), max(get_mtime(f) for f in filein))
    mtime_out = get_mtime(fileout)
    if mtime_in < mtime_out: return
    merge_root_files(filein, fileout)
