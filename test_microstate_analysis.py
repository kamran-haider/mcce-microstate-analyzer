"""Tester script for microstate analysis module and associated functions.
"""
from utils import generate_cluster_dictionary, generate_path_pdb
from microstate_analysis import MicrostateAnalysis, ConformerHbonds, HydrogenBondedPath
from progressbar import Bar, Percentage, ProgressBar, ETA


data_dir = "/Users/kamranhaider/ms_dat_testcase/"
substate_list = ["f1", "f2", "f4"]
snapshot_list = ["f100"]#, "f200", "f400"]
sample_freq = 100000
test_residues = [('THRA0100', 'ASPA0485')]

for snapshot in snapshot_list:
    paths = [HydrogenBondedPath(test[0], test[1]) for test in test_residues]
    substate = snapshot[0:2]
    print "substate: ", substate
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
    #msa.write_hbtxt(prefix=snapshot + "_")
    for p in paths:
        #print len(p.microstates_with_path) + len(p.microstates_without_path)
        shortest_path, microstates_with_path = p.filter_paths()
        if shortest_path is not None and microstates_with_path is not None:
            print "writing pdb file for path."
            path_conformers = [msa.conformer_data[conf][0] for conf in shortest_path]
            print path_conformers
            generate_path_pdb(path_conformers, step2out, prefix=p.start_residue + "_" + p.end_residue)
        else:
            print "No paths found!, skipping wirting PDB."