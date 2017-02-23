"""Tester script for microstate analysis module and associated functions.
"""
from __future__ import print_function
from __future__ import division
from builtins import range
from past.utils import old_div
import numpy as np

from utils import generate_cluster_dictionary, generate_path_pdb, print_progress_bar
from microstate_analysis import MicrostateAnalysis, ConformerHbonds, HydrogenBondedPath

def process_snapshots(paths, snapshot, sample_freq=100):
    """
    Parameters
    ----------
    paths : list
        A list of tuples, each of which contains start and end residue corresponding to
        a given snapshot
    snapshot : string
        Name of the snapshot, used as a prefix to all major input files
    """
    #paths = [HydrogenBondedPath(test[0], test[1]) for test in test_residues]
    substate = snapshot[0:2]
    clusters = generate_cluster_dictionary(substate)
    # generate filenames
    hbdat = data_dir + snapshot + "_hb.dat"
    msdat = data_dir + snapshot + "_ms.dat"
    head3lst = data_dir + snapshot + "_head3.lst"
    step2out = data_dir + snapshot + "_step2_out.pdb"
    hbtxt = snapshot + "_hb.txt"
    
    # initialie MicrostateAnalysis object
    msa = MicrostateAnalysis(msdat, head3lst)
    msa.generate_byte_indices(sample_frequency=sample_freq)
    msa.parse_records()
    # initialize conformer H-bond matrix
    confhb = ConformerHbonds(hbdat)
    confhb.generate_hbmatrix_from_hbdat()

    progress_counter = 0
    max_items = msa.trajectory.shape[0]
    #print(len(msa.residue_list))
    print(msa.trajectory.shape)
    #print(len(msa.conformer_data.keys()))
    print(msa.conformer_data[confhb.hb_matrix.shape[0]])
    #print_progress_bar(progress_counter, max_items)
    test_conf_res =  [msa.conformer_data[i][0] for i in range(10, 18)]
    res_index = msa.residue_list.index('ARGA0019')
    print(test_conf_res)
    print(res_index)
    conf_of_res = msa.trajectory[:, 530]
    print(conf_of_res)
    print([msa.conformer_data[k][0] for k in conf_of_res])


 

data_dir = "/Users/kamranhaider/ms_dat_testcase/"
substate_list = ["f1", "f2", "f4"]
#snapshot_list = ["f200", "f400"]
snapshot_list = ["f100"]
sample_freq = 100000
for snapshot in snapshot_list:
    print(snapshot)
    print("start end path path_length total_microstates occupancy possible_paths")
    cluster_dict = generate_cluster_dictionary(snapshot[0:2])
    source_residues = cluster_dict["EXT_INT"]
    target_residues = cluster_dict["EXT_EXP"]
    paths = [HydrogenBondedPath(source_residue, target_residue) for source_residue in source_residues
                                                for target_residue in target_residues ]

    process_snapshots(paths, snapshot, sample_freq)
