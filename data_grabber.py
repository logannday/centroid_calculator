import os, argparse, re
from Bio.PDB.PDBParser import PDBParser
# from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.internal_coords import *
# from Bio.PDB.PICIO import write_PIC, read_PIC, read_PIC_seq
# from Bio.PDB.ic_rebuild import write_PDB, IC_duplicate, structure_rebuild_test
# from Bio.PDB.SCADIO import write_SCAD
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.PDB.PDBIO import PDBIO
import pandas as pd
# import numpy as np
from centroid import Target, calculate_centroid, closest_interface_residue
from dataclasses import dataclass

@dataclass
class mutant:
    structure: Structure
    insertion: bool
    location: int
    residue: str

def get_structures(path):
    parser = PDBParser()
    myProtein = parser.get_structure("1ak4", "1ak4.pdb")
    path = "./1brs_PDBS/"
    filenames = os.listdir(path)
    mutants = []
    pattern = r'^(\w+)(?:_ins)?_(\d+)_([A-Za-z])\.pdb$'

    # Parse each filename, and create an array of mutant structures
    for filename in filenames:
        match = re.match(pattern, filename)
        if not match:
            raise ValueError(f"Filename {filename} doesn't match the expected format.")

        pdb_name = match.group(1)
        insertion = '_ins' in filename
        residue = int(match.group(2))
        chain = match.group(3)
        struct = parser.get_structure(pdb_name, path + filename)
        mutants.append(mutant(struct, True, residue, chain))
    return mutants

def get_data(mutants):
    target = Target("A", [5, 6, 7, 9, 10])
    targets = [target]
    rows = []
    for mutant in mutants:
        centroid = calculate_centroid(mutant.structure, targets)
        print(mutant)
        rows.append({ 
            "chain": mutant.chain,
        })
    df = pd.DataFrame.from_dict(rows)

def main():
    mutants = get_structures("./1brs_PDBS")
    data = get_data(mutants)

if __name__ == "__main__":
    main()
