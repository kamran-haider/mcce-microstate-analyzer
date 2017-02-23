
"""
Moduule containing useful classes for the analysis of eqiulibrium distribution of 
conformatioons and protomers, stored as microstate data files.  
"""
from __future__ import print_function
from __future__ import division
from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
import os
import struct

import numpy as np
import networkx as nx
from utils import function_timer, print_progress_bar, is_grotthuss

class MicrostateAnalysis(object):
    """A class representing microstate data and various analysis methods"""
    def __init__(self, ms_data_file, head3lst_file):
        self.ms_data_file = ms_data_file
        self.byte_indices = None
        self.total_microstates = 0
        self.total_records = 0

        res_list = []
        with open(self.ms_data_file, "rb") as md:
            bytes_n_res = md.read(4)
            n_res = struct.unpack('i', bytes_n_res)[0]
            for i in range(n_res):
                resname = str(md.read(8), 'utf-8')
                res_list.append(resname)
            self.n_res = n_res
        self.residue_list = res_list
        self.residue_hb_matrix = np.zeros((n_res, n_res), dtype=float)

        #with open(ms_gold_file, "r") as ms_gold:
        #    self.n_res = len([res.strip() for res in ms_gold.readlines() if len(res) != 0])
        conf_data = {}
        with open(head3lst_file, "r") as h3:
            for line in h3.readlines()[1:]:
                data = line.split()
                conf_id = int(data[0])
                conf_data[conf_id] = [data[1], data[3]]
        self.conformer_data = conf_data

    def generate_byte_indices(self, sample_frequency=100, filename=None):
        """

        Parameters
        ----------
        n_res : TYPE
            Description
        sample_frequency : int, optional
            Description
        filename : None, optional
            Description
        
        Returns
        -------
        rec_indices : list
            A list of the starting bytes for each record in the microstate data file
        """
        if filename is None:
            filename = self.ms_data_file

        start_byte = 4 + (8 * self.n_res)
        bytes_per_record = (self.n_res * 2) + 20
        file_size = os.path.getsize(filename)
        # n_records = (file_size - start_byte) / bytes_per_record
        rec_indices = list(range(start_byte, file_size, sample_frequency * bytes_per_record))
        self.total_records = len(rec_indices)
        self.byte_indices = rec_indices
        
    @function_timer
    def parse_records(self):
        """
        
        """
        trajectory = np.zeros([self.total_records, self.n_res], dtype=int)
        state_counts = np.zeros([self.total_records], dtype=int)
        energies = np.zeros([self.total_records], dtype=float)
        #print trajectory
        progress_counter = 0
        print_progress_bar(progress_counter, self.total_records)
        with open(self.ms_data_file, "rb") as ms:
            for index, record in enumerate(self.byte_indices):
                ms.seek(record)
                bytes_conf_ids = ms.read(2 * self.n_res)
                bytes_energies_1 = ms.read(8)
                ms.seek(ms.tell() + 8)
                energy = struct.unpack("d", bytes_energies_1)[0]
                bytes_state_count = ms.read(4)
                trajectory[index, :] = np.asarray(struct.unpack(str(self.n_res) + "H", bytes_conf_ids))
                #print(struct.unpack(str(self.n_res) + "H", bytes_conf_ids)[-2:])
                state_count = struct.unpack("i", bytes_state_count)[0]
                self.total_microstates += state_count
                state_counts[index] += state_count
                energies[index] += energy
                progress_counter += 1
                print_progress_bar(progress_counter, self.total_records)
        self.trajectory = trajectory
        self.state_counts = state_counts
        self.energies = energies

    def write_hbtxt(self, prefix= ""):
        THRESHOLD_TO_WRITE = 0.001
        hbtxt_file = prefix + "hb.txt"
        print("Writing data to %s" % hbtxt_file)
        with open(hbtxt_file, "w") as t:
            for i in range(len(self.residue_list)):
                for j in range(len(self.residue_list)):
                    if old_div(self.residue_hb_matrix[i, j], self.total_microstates) >= THRESHOLD_TO_WRITE:
                        t.write(
                            "%s\t%s\t%.3f\n" %
                            (self.residue_list[i],
                             self.residue_list[j],
                             old_div(self.residue_hb_matrix[i, j], self.total_microstates)))


class ConformerHbonds(object):
    """A class to represent conformer hydrogen bond matrix and associated methods.
    """
    def __init__(self, hbdat_file):
        self.hbdat_file = hbdat_file

    @function_timer
    def generate_hbmatrix_from_hbdat(self, hbdat_file=None):
        """
        Parameters
        ----------
        hbdat_file : None, optional
            Description
        """
        if hbdat_file is None:
            hbdat_file = self.hbdat_file

        with open(hbdat_file, "rb") as hbhandle:
            bytes_n_conf = hbhandle.read(4)
            n_conf = struct.unpack('i', bytes_n_conf)[0]
            #print "Total conformers: ", n_conf
            bytes_h_matrix = hbhandle.read(n_conf * n_conf)
            hb_array = struct.unpack(str(n_conf * n_conf) + "B", bytes_h_matrix)
            self.hb_matrix = np.array(hb_array).reshape(n_conf, n_conf)
            
    def generate_microstate_hb_network(self, conformer_ids=None, labels=None):
        """
        Parameters
        ----------
        conformer_ids : TYPE
            Description
        """
        if conformer_ids is None:
            G = nx.to_networkx_graph(self.hb_matrix)
            return G
        else:
            G = nx.to_networkx_graph(self.hb_matrix[conformer_ids, :][:, conformer_ids])
            if labels is None:
                nx.relabel_nodes(G, dict(list(zip(G.nodes(), conformer_ids))), copy=False)
            else:
                assert len(labels) == len(conformer_ids), "Please make sure number of nodes and labels are equal."
                nx.relabel_nodes(G, dict(list(zip(G.nodes(), labels))), copy=False)
            return G



    def check_conformer_hbonds_pairwise(self, conformer_a, conformer_b):
        """
        Parameters
        ----------
        conformer_a : TYPE
            Description
        conformer_b : TYPE
            Description
        """
        pass

    def all_hbonds_by_conformer(conformer):
        pass

        
class HydrogenBondedPath(object):
    """
    """
    def __init__(self, start_residue, end_residue):
        """
        Parameters
        ----------
        start_residue : TYPE
            Description
        end_residue : TYPE
            Description
        """
        self.start_residue = start_residue
        self.end_residue = end_residue
        self.states_analyzed = 0
        self.path_list = []
        self.path_lengths = []
        self.record_indices = []

    def search_for_path(self, hb_network, record_index, start, end):
        """
        Parameters
        ----------
        hb_network : TYPE
            Description
        start : TYPE
            Description
        end : TYPE
            Description
        """
        try:
            path = nx.shortest_path(hb_network, source=start, target=end)
            self.path_list.append(path)
            self.path_lengths.append(len(path))
            self.record_indices.append(record_index)

        except nx.NetworkXNoPath as nopath:
            self.path_list.append(None)
            self.path_lengths.append(np.nan)
            self.record_indices.append(np.nan)

        finally:
            self.states_analyzed += 1

    def filter_paths(self, conformer_data):
        """
        Parameters
        ----------
        conformer_data : TYPE
            Description
        
        """
        path_lengths_array = np.asarray(self.path_lengths)
        path_list_array = np.asarray(self.path_list)
        record_indices_array = np.asarray(self.record_indices)
        filtered_paths = []
        filtered_path_records = []

        min_length = np.nanmin(path_lengths_array)
        shortest_paths = path_list_array[np.where(path_lengths_array == min_length)]
        shortest_path_records = record_indices_array[np.where(path_lengths_array == min_length)]
        for index, path in enumerate(shortest_paths):
            non_grotthuss = [conformer_data[conf][0] for conf in path[1:-1] if not is_grotthuss(conformer_data[conf][0][0:3])]
            #print non_grotthuss
            if len(non_grotthuss) == 0:
                filtered_paths.append(list(path))
                filtered_path_records.append(int(shortest_path_records[index]))
        return filtered_paths, filtered_path_records
