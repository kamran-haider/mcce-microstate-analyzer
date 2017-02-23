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
    print_progress_bar(progress_counter, max_items)
    for index in range(msa.trajectory.shape[0]):
        current_conformers = msa.trajectory[index, :]
        current_hb_network = confhb.generate_microstate_hb_network(current_conformers)
        msa.residue_hb_matrix += (confhb.hb_matrix[current_conformers, :][:, current_conformers] * msa.state_counts[index])
        for p in paths:
            # start and end are correct ids
            start = current_conformers[msa.residue_list.index(p.start_residue)]
            end = current_conformers[msa.residue_list.index(p.end_residue)]
            #print p.start_residue, msa.conformer_data[start]
            #print p.end_residue, msa.conformer_data[end]
            p.search_for_path(current_hb_network, index, start, end)
        progress_counter += 1
        print_progress_bar(progress_counter, max_items)

    msa.write_hbtxt(prefix=snapshot + "_")


 

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
