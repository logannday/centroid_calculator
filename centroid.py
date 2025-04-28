import os, argparse, re
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.internal_coords import *
from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test
from Bio.PDB.SCADIO import write_SCAD
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.PDBIO import PDBIO
import numpy as np

class Target:
    def __init__(self, chain, residues):
        self.chain = chain
        self.residues = residues

def calculate_centroid(structure, targets):
    model = structure[0]

    xs = []
    ys = []
    zs = []
    for target in targets:
        chain = model[target.chain]
        for i in target.residues:
            coords = chain[i]["CA"].get_coord()
            xs.append(coords[0])
            ys.append(coords[1])
            zs.append(coords[2])


    centroid = [np.average(xs), np.average(ys), np.average(zs)]

    return centroid

def closest_interface_residue(structure, targets, res):
    model = structure[0]
    
    alpha_carbon = model[res.chain][res.residues[0]]["CA"] # alpha carbon of residue of interest
    min_distance = np.inf

    # Finding the closest interface residue to the residue of interest
    for target in targets: 
        chain = model[target.chain]
        for i in target.residues:
            for atom in chain[i]:
                min_distance = np.minimum(min_distance, atom - alpha_carbon)

    return min_distance
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb-id', type=str)
    parser.add_argument('pdb-filepath', type=str, help = "filepath to pdb. + ex: \"./1ak4.pdb\"")
    parser.add_argument('--targets', dest='targets', type=str, default="A", required=True, 
                        help = "String defining the interface residues. Example:" +
                        "A1,12,13,15B14,15 indicates chain A residues 1, 12, 13, 14" +
                        "and chain B residues 14 and 15")

    pattern = r'([A-Za-z])((?:\d+,?)*)'
    targets = "A1,12,13,15B14,15"
    matches = re.findall(pattern, targets)

    target_objects = []
    for chain_id, residues_str in matches:
        residues = [int(r) for r in residues_str.split(',') if r]
        target_objects.append(Target(chain_id, residues))

    for obj in target_objects:
        print(obj.chain, obj.residues)

    # load a structure as normal, get first chain
    parser = PDBParser()
    myProtein = parser.get_structure("1ak4", "1ak4.pdb")
    # myChain = myProtein[0]["A"]
    centroid = calculate_centroid(myProtein, target_objects)
    min_distance = closest_interface_residue(myProtein, target_objects, Target("A", [25]))
    print("Centroid: ", centroid)
    print("Min distance:", min_distance)


if __name__ == "__main__":
    main()
