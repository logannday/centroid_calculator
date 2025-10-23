from mutant_loading import *
import os

# from auto_normal import plane_fig
from typing import List
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.PDB.internal_coords import *
from dataclasses import dataclass


@dataclass
class Mutant:
    structure: Structure
    insertion: bool
    location: int
    residue: str
    filename: str


def get_mutant(base_path, filename) -> Mutant:
    pattern = r"^(\w+)(?:_ins)?_(\d+)_([A-Za-z])\.pdb$"
    match = re.match(pattern, filename)
    if not match:
        raise ValueError(f"Filename {filename} doesn't match the expected format.")
    insertion = "_ins" in filename
    pdb_id = match.group(1)
    ins_location = int(match.group(2))
    resname = match.group(3)
    parser = PDBParser()
    struct = parser.get_structure(pdb_id, os.path.join(base_path, filename))
    return Mutant(struct, insertion, ins_location, resname, filename)


def get_mutants(structs_names):
    """
    given a list of tuples (struct, filename)
     parse the struct and insertion_or_deletion, mutation location,
     and residue name from the filename into a
     list of Mutant objects

    structs_names: list of tuples (Struture, Filename)
    """
    mutants = []
    pattern = r"^(\w+)(?:_ins)?_(\d+)_([A-Za-z])\.json$"

    # Parse each filename, and create an array of mutant structures
    for struct, filename in structs_names:
        match = re.match(pattern, filename)
        if not match:
            raise ValueError(f"Filename {filename} doesn't match the expected format.")
        insertion = "_ins" in filename
        ins_location = int(match.group(2))
        resname = match.group(3)
        mutants.append(Mutant(struct, insertion, ins_location, resname, filename))
    return mutants


def get_mutants_from_pdbs(base_path) -> List[Mutant]:
    """
    given a list of tuples (struct, filename)
     parse the struct and insertion_or_deletion, mutation location,
     and residue name from the filename into a
     list of Mutant objects
    """
    parser = PDBParser()
    mutants = []
    pattern = r"^(\w+)(?:_ins)?_(\d+)_([A-Za-z])\.pdb$"
    for filename in os.listdir(base_path):
        print(filename)
        match = re.match(pattern, filename)

        if not match:
            raise ValueError(f"Filename {filename} doesn't match the expected format.")

        insertion = "_ins" in filename
        pdb_id = match.group(1)
        ins_location = int(match.group(2))
        resname = match.group(3)
        struct = parser.get_structure(pdb_id, os.path.join(base_path, filename))
        mutants.append(Mutant(struct, insertion, ins_location, resname, filename))

    return mutants
