from typing import List

# from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.internal_coords import *
from Bio.PDB.Polypeptide import is_aa
import numpy as np


def calculate_centroid(residues: List[Residue]):
    """
    Calculate the coordinates of the centroid of a list of their residues
    """
    xs = []
    ys = []
    zs = []

    for res in residues:
        if is_aa(res):
            coords = res["CA"].get_coord()
            xs.append(coords[0])
            ys.append(coords[1])
            zs.append(coords[2])

    centroid = [np.average(xs), np.average(ys), np.average(zs)]

    return centroid


def closest_res_distance(interface_residues: List[Residue], target_res: Residue):
    """
    Return the distance in Angstroms from target_res to the closest residue to target_res
    in interface residues
    """
    min_distance = np.inf

    # Finding the closest interface residue to the residue of interest
    for i_res in interface_residues:
        try:
            min_distance = np.minimum(min_distance, i_res["CA"] - target_res["CA"])
        except:
            print("residue didn't have alpha carbon")

    return min_distance


def get_ca_atoms(residues):
    """
    Extracts CA atoms from a list of residues
    """
    return [res["CA"] for res in residues]


def calculate_rmsd(wt_residues: List[Residue], mut_residues: List[Residue]):
    """
    Calculate the RMSD of two sets of residues using their alpha carbons
    """
    # TODO: Verify each atom name matches if
    supermimposer = Superimposer()
    print(wt_residues)
    print(mut_residues)
    wt_atoms = get_ca_atoms(wt_residues)
    mut_atoms = get_ca_atoms(mut_residues)
    if len(wt_atoms) != len(mut_atoms):
        return np.nan

    supermimposer.set_atoms(wt_atoms, mut_atoms)
    supermimposer.apply(mut_atoms)
    return supermimposer.rms
