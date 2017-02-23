"""Tester script for microstate analysis module and associated functions.
"""
from __future__ import print_function
from builtins import str
import numpy as np

from utils import generate_cluster_dictionary, generate_path_pdb
from microstate_analysis import MicrostateAnalysis, ConformerHbonds, HydrogenBondedPath
from progressbar import Bar, Percentage, ProgressBar, ETA




data_dir = "/Users/kamranhaider/ms_dat_testcase/"
substate_list = ["f1", "f2", "f4"]
#snapshot_list = ["f100", "f200", "f400"]
snapshot_list = ["f100"]
sample_freq = 1
for snapshot in snapshot_list:
    print(snapshot)
    msdat = data_dir + snapshot + "_ms.dat"
    head3lst = data_dir + snapshot + "_head3.lst"
    # initialie MicrostateAnalysis object
    msa = MicrostateAnalysis(msdat, head3lst)
    msa.generate_byte_indices(sample_frequency=sample_freq)
    msa.parse_records()
    with open(snapshot + "_energies.txt", "w") as f:
        for e in msa.energies:
            f.write(str(e) + "\n")
