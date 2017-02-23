"""Tester script for microstate analysis module and associated functions.
"""
from __future__ import print_function
from __future__ import division
from past.utils import old_div
import numpy as np

from utils import generate_cluster_dictionary, generate_path_pdb
from microstate_analysis import MicrostateAnalysis, ConformerHbonds, HydrogenBondedPath
from progressbar import Bar, Percentage, ProgressBar, ETA

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
    #msa.parse_records()
    # initialize conformer H-bond matrix
    confhb = ConformerHbonds(hbdat)
    confhb.generate_hbmatrix_from_hbdat()
    
    print(confhb)
    wat_a = "S0864"
    wat_b = "R4956"
    wat_a_conf_ids = []
    wat_b_conf_ids = []
    #HOH01S0864_004->HOH01R4956_054
    for k in list(msa.conformer_data.keys()):
        if msa.conformer_data[k][0][5:10] == wat_a:
            wat_a_conf_ids.append(k)
        if msa.conformer_data[k][0][5:10] == wat_b:
            wat_b_conf_ids.append(k)

    print(len(wat_a_conf_ids), len(wat_b_conf_ids))
    #hb_network = confhb.generate_microstate_hb_network()
    hb_count = 0
    tot_possible_hb = 0
    for i in wat_a_conf_ids:
        for j in wat_b_conf_ids:
            print(msa.conformer_data[i][0], msa.conformer_data[j][0], confhb.hb_matrix[i][j])
            tot_possible_hb += 1 
            if confhb.hb_matrix[i][j] == 1:
                hb_count += 1
    print(hb_count, tot_possible_hb, old_div(float(hb_count),float(tot_possible_hb)))
    """
    pbar = ProgressBar(widgets=[Percentage(), Bar(), ETA()], 
                        maxval=msa.trajectory.shape[0]).start()

    for index in range(msa.trajectory.shape[0]):
        current_conformers = msa.trajectory[index, :]
        current_hb_network = confhb.generate_microstate_hb_network(current_conformers)
        msa.residue_hb_matrix += (confhb.hb_matrix[current_conformers, :][:, current_conformers] * msa.state_counts[index])
        for p in paths:
            start = current_conformers[msa.residue_list.index(p.start_residue)]
            end = current_conformers[msa.residue_list.index(p.end_residue)]
            #print p.start_residue, msa.conformer_data[start]
            #print p.end_residue, msa.conformer_data[end]
            p.search_for_path(current_hb_network, index, start, end)
        pbar.update(index)
    pbar.finish()
    msa.write_hbtxt(prefix=snapshot + "_")
    for p in paths:
        print p.start_residue, p.end_residue,
        shortest_paths, records_with_path = p.filter_paths(msa.conformer_data)
        if len(shortest_paths) != 0:
            shortest_path_conformers = [msa.conformer_data[conf][0] for conf in shortest_paths[0]]    
            n_microstates_path = np.sum(msa.state_counts[records_with_path])
            occ = float(n_microstates_path)/float(msa.total_microstates)
            e_microstates = np.mean(msa.energies)
            e_microstates_path = np.mean(msa.energies[records_with_path])
            print "->".join(shortest_path_conformers), len(shortest_path_conformers),
            print n_microstates_path, msa.total_microstates, occ,
            # obtain energies and occupancies
            print e_microstates, e_microstates_path
            generate_path_pdb(shortest_path_conformers, step2out, prefix=snapshot + "_" + p.start_residue + "_" + p.end_residue)
        else:
            #print "No paths found!, skipping wirting PDB."
            print
            pass

    """
data_dir = "/Users/kamranhaider/ms_dat_testcase/"
substate_list = ["f1", "f2", "f4"]
#snapshot_list = ["f100", "f200", "f400"]
snapshot_list = ["f100"]
sample_freq = 1000
for snapshot in snapshot_list:
    print(snapshot)
    print("start end path path_length microstates_path total_microstates occupancy avg_e_states avg_e_states_path")
    cluster_dict = generate_cluster_dictionary(snapshot[0:2])
    source_residues = cluster_dict["EXT_INT"]
    target_residues = cluster_dict["EXT_EXP"]
    paths = [HydrogenBondedPath(source_residue, target_residue) for source_residue in source_residues
                                                for target_residue in target_residues ]

    process_snapshots(paths, snapshot, sample_freq)
