import os, re
from typing import List

from plotting import plot_centroid_vs_rmsd
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure
from Bio.Align import PairwiseAligner
from Bio.PDB.StructureAlignment import StructureAlignment
from Bio.PDB.Residue import Residue
from Bio.PDB.internal_coords import *
import pandas as pd
from centroid import calculate_centroid, calculate_rmsd, closest_interface_residue
from detect_interface import get_interface_residues
from dataclasses import dataclass

aa_names = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U",  # Selenocysteine
    "PYL": "O",  # Pyrrolysine
    "ASX": "B",  # Aspartic acid or Asparagine
    "GLX": "Z",  # Glutamic acid or Glutamine
    "XAA": "X"   # Unknown or other
}


@dataclass
class Mutant:
    structure: Structure
    insertion: bool
    location: int
    residue: str


def get_mutants(path):
    parser = PDBParser()
    filenames = os.listdir(path)
    mutants = []
    pattern = r"^(\w+)(?:_ins)?_(\d+)_([A-Za-z])\.pdb$"

    # Parse each filename, and create an array of mutant structures
    for filename in filenames:
        match = re.match(pattern, filename)
        if not match:
            raise ValueError(f"Filename {filename} doesn't match the expected format.")

        pdb_name = match.group(1)
        insertion = "_ins" in filename
        ins_location = int(match.group(2))
        resname = match.group(3)
        struct = parser.get_structure(pdb_name, path + filename)
        mutants.append(Mutant(struct, insertion, ins_location, resname))
    return mutants

def get_res_by_absolute_index(struc: Structure, location: int, chains = ("A", "B")) -> Residue:
    """
        get residue at the absolute index given in the protein
        struc = Biopython Structure object
        location = Index of amino acid to return, the absolute index starts at 0 and increments
                   over all the chains
        chains = The chains which are included in the complex of interest
    """
    chain1 = struc[0][chains[0]]
    chain2 = struc[0][chains[1]]

    # Filter out non- amino acid residues
    chain1_res = [res for res in chain1 if is_aa(res)]
    chain2_res = [res for res in chain2 if is_aa(res)]
    chain1_len = len(chain1)

    if chain1_len > location:
        return chain1_res[location]
    return struc[0][chains[1]][location - chain_a_len]
    chain2_res[location]

def get_fasta(struc: Structure, chain_ids = ("A", "B")):
    """
    Return the fasta sequence for the specified chains of a given structure
    """
    model = struc[0]
    sequence = []
    for chain in chain_ids:
        # get single letter code for each residue in the chain
        residues = [aa_names.get(res.resname) for res in model[chain] if res.resname != "HOH"]
        sequence = sequence + residues
    
    return "".join(sequence)


def map_residues(ref_struc: Structure, target_struc: Structure, wt_residues: List[Residue], wt_chains = ("A", "B")):
    """
        map a set of residues in one structure to the same residues in the mutant. 
        ref_struc: The structure of the wildtype
        target_struc: The structure of the mutant
        wt_chains: determines which chains are equivalent to chaines A and B in the mutant. This is necessary because
        the pipeline renames the chains when it produces the mutants
    """
    wt_seq = get_fasta(ref_struc, wt_chains)
    mut_seq = get_fasta(target_struc)

    aligner = PairwiseAligner()
    alignments = aligner.align(wt_seq, mut_seq)
    print(alignments[0])
    alignment = alignments[0]
    struture_alignment = StructureAlignment(alignment, ref_struc, target_struc)
    mut_interface_residues = []

    (res_map, _) = struture_alignment.get_maps()
    mut_interface_residues = [res_map.get(res) for res in wt_residues]

    return mut_interface_residues

def populate_dataframe(wildtype: Structure, mutants: List[Mutant], wt_chains = ("A", "B")):
    wt_interface_residues = get_interface_residues(wildtype, chain_ids=wt_chains)
    rows = []

    for mutant in mutants:
        mut_interface_residues = map_residues(wildtype, mutant.structure, wt_interface_residues, wt_chains)
        centroid = calculate_centroid(wt_interface_residues)
        inserted_residue = get_res_by_absolute_index(wildtype, mutant.location, wt_chains)
        cir = closest_interface_residue(wt_interface_residues, inserted_residue)
        centroid_distance = np.linalg.norm(
            inserted_residue["CA"].get_coord() - centroid
        )
        interface_rmsd = calculate_rmsd(wt_interface_residues, mut_interface_residues)
        rows.append(
            {
                "location": mutant.location,
                "residue": mutant.residue,
                "centroid_distance": centroid_distance,
                "closest_interface_res_distance": cir,
                "interface_rmsd": interface_rmsd,
            }
        )
    df = pd.DataFrame.from_dict(rows)
    return df


def main():
    wt_chains = ("C", "F")
    pdb_id = "1brs"
    parser = PDBParser()
    wildtype = parser.get_structure(pdb_id, f"./{pdb_id}/in.pdb")
    mutants = get_mutants("./1brs_PDBS/")
    data = populate_dataframe(wildtype, mutants, wt_chains)
    plot_centroid_vs_rmsd(data, "test", "./rmsd_data/test.png")


if __name__ == "__main__":
    main()
