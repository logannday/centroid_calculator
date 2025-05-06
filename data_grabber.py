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
from centroid import calculate_centroid, closest_interface_residue
from detect_interface import get_interface_residues
from dataclasses import dataclass

@dataclass
class mutant:
    structure: Structure
    insertion: bool
    location: int
    residue: str

def get_structures(path):
    parser = PDBParser()
    myProtein = parser.get_structure("1brs", "1brs.pdb")
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

# Get the residue defined at the absolute location in the pdb
def get_res_by_absolute_index(struc, location):
    chain_a_len = len(list(struc[0]['A'].get_residues()))
    if chain_a_len > location:
        return struc[0]['A'][location]
    return struc[0]['B'][location - chain_a_len]

def populate_dataframe(wildtype: Structure, mutants):
    rows = []
    for mutant in mutants:
        interface_residues = get_interface_residues(mutant, ["C", "F"])
        centroid = calculate_centroid(interface_residues)
        inserted_residue = get_res_by_absolute_index(mutant.structure, mutant.location)
        cir = closest_interface_residue(interface_residues, inserted_residue)
        print(inserted_residue.get_coord())
        breakpoint()
        rows.append({ 
            "location": mutant.location,
            "residue": mutant.residue,
        })
    df = pd.DataFrame.from_dict(rows)
    return df

def main():
    parser = PDBParser()
    wildtype = parser.get_structure("1brs", "./1brs.pdb")
    for chain in wildtype[0]:
        print(chain.id)
    mutants = get_structures("./1brs_PDBS")
    data = populate_dataframe(wildtype, mutants)
    # print(data.head())

if __name__ == "__main__":
    main()
