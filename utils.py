"""Some utility functions for microstate hydrogen bond analysis.
"""
from __future__ import print_function

# Globals
import sys
import time
from functools import wraps

GROTTHUSS = ["GLU", "ASP", "SER", "THR", "HIS", "TYR", "ASN", "GLN", "HOH", "PAA", "PDD"]

def is_grotthuss(residue):
    return residue in GROTTHUSS

# Functions
def function_timer(function):
    @wraps(function)
    def function_timer(*args, **kwargs):
        t0 = time.time()
        result = function(*args, **kwargs)
        t1 = time.time()
        print("Total time running %s: %2.2f seconds" %
              (function.__name__, t1 - t0))
        return result
    return function_timer

def print_progress_bar (count, total):
    """
    Create and update progress bar during a loop.
    
    Parameters
    ----------
    iteration : int
        The number of current iteration, used to calculate current progress. 
    total : int
        Total number of iterations
    
    Notes
    -----
        Based on:
        http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = u"\u2588" * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('Progress |%s| %s%s Done\r' % (bar, percents, '%'))
    sys.stdout.flush() 
    if count == total: 
        print()



def generate_cluster_dictionary(sub_state): 
    """
    Generate a dictionary containing key clusters and corresponding residues
    for a given sub_state.  
    Parameters
    ----------
    sub_state : string
        Specifify
    """
    # Lists of residues in each cluster
    clust_dict = {  
    "BNC" : ['CUBK9500', 'HE3K9300'],
    "PLS" : ['PDDK9102', 'PAAK9301', 'PDDK9302', 'TYRA0175',
            'ASPA0407', 'HISA0411', 'ASPB0229', 'ASPA0412', 'THRA0337', 'GLUB0254',
           'SERB0218', 'PAAK9101', 'ARGA0052', 'TYRA0414', 'SERA0497', 'SERA0498'],
    "EXT" : ['HISA0093', 'SERA0156', 'THRA0187', 'SERA0168', 'THRA0100',
            'TYRB0262', 'SERA0186', 'GLUA0182', 'ASPA0188', 'ASPA0485', 'GLUA0488'],
    "GLU" : ["GLUA0286"],
    "EXT_EXP" : ['HISA0093', 'GLUA0182', 'ASPA0188', 'ASPA0485', 'GLUA0488'],
    "EXT_INT" : ['SERA0168', 'THRA0100']
    }

    if sub_state not in ["f1", "f2", "f4"]:
        sys.exit("Only f1, f2 and f4 substates are currently supported.")

    if sub_state == "f1" or sub_state == "f4":
        return clust_dict
    else:
        for k in list(clust_dict.keys()):
            updated_cluster_residues = []
            for res in clust_dict[k]:
                if res[3] == "A":
                    resnum = int(res[4:])
                    modified_resnum = resnum - 11
                    str_resnum = "%04d" % (modified_resnum)
                    modified_res = res[0:4] + str_resnum
                    updated_cluster_residues.append(modified_res)
                else:
                    updated_cluster_residues.append(res)
            clust_dict[k] = updated_cluster_residues
        return clust_dict

def generate_path_pdb(conformers, source_pdb, prefix= "path"):
    """generate a pdb file from path data
    
    Parameters
    ----------
    conformers : TYPE
        Description
    source_pdb : TYPE
        Description
    prefix : str, optional
        Description
    """
    with open(source_pdb, 'r') as pdb:
        pdb_lines = pdb.readlines()
        with open(prefix + ".pdb", "w") as f:
            for conf in conformers:
                for l in pdb_lines:
                    if l.split()[3] == conf[0:3]:
                        #print l.split()[3], l.split()[4]
                        #print conf[0:3], conf[5:]    
                        if l.split()[4] == conf[5:] or l.split()[4] == conf[5:11] + "000":
                            #print l.split()[3], l.split()[4]
                            #print conf[0:3], conf[5:]
                            f.write(l)




