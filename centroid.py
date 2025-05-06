from typing import List
# from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB.Chain import Chain
# from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.internal_coords import *
# from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
# from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test
# from Bio.PDB.SCADIO import write_SCAD
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.PDB.PDBIO import PDBIO
import numpy as np

def calculate_centroid(residues: List[Residue]):
    xs = []
    ys = []
    zs = []

    for res in residues:
        coords = res["CA"].get_coord()
        xs.append(coords[0])
        ys.append(coords[1])
        zs.append(coords[2])

    centroid = [np.average(xs), np.average(ys), np.average(zs)]

    return centroid

def closest_interface_residue(interface_residues, target_res):
    min_distance = np.inf

    # Finding the closest interface residue to the residue of interest
    for i_res in interface_residues:
            min_distance = np.minimum(min_distance, i_res["CA"] - res["CA"])

    return min_distance
