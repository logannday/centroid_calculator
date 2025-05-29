from typing import List
import sys, math

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
    points = np.array([res["CA"].get_coord() for res in residues]).T
    centroid = np.mean(points, axis=1, keepdims=True)
    return centroid


def plane_normal(residues: List[Residue]):
    # generate random points on the plane and add random displacement
    points = np.array([res["CA"].get_coord() for res in residues]).T

    # subtract out the centroid and take the SVD
    svd = np.linalg.svd(points - np.mean(points, axis=1, keepdims=True))

    # Extract the left singular vectors
    left = svd[0]
    # the normal vector to the plane
    normal = left[:, -1]
    return normal


def vec_angle(normal1: np.ndarray, normal2: np.ndarray) -> float:
    """
    Compute the angle (in degrees) between two planes defined by their normal vectors.
    """
    # Ensure both normals are unit vectors
    n1 = normal1 / np.linalg.norm(normal1)
    n2 = normal2 / np.linalg.norm(normal2)

    # Compute the dot product and clamp for numerical stability
    dot = np.clip(np.dot(n1, n2), -1.0, 1.0)

    # Return angle in degrees
    angle_rad = np.arccos(abs(dot))
    angle_deg = np.degrees(angle_rad)
    return angle_deg


def face_angle(r1: List[Residue], r2: List[Residue]) -> float:
    """
    compute the angle between the two faces of an interface
    r1: the residues in one participating chain
    r1: the residues in the other participating chain
    """
    norm1 = plane_normal(r1)
    norm2 = plane_normal(r2)
    return vec_angle(norm1, norm2)


def centroid_radius(residues, centroid: List[float]) -> float:
    """
    calculate the radius of the centroid by finding the furthest interface residue from the centroid
    TODO: Check for outliers
    """
    cx, cy, cz = centroid[0], centroid[1], centroid[2]
    radius = 0
    for res in residues:
        if is_aa(res):
            coords = res["CA"].get_coord()
            x, y, z = (coords[0], coords[1], coords[2])
            dist = math.sqrt((cx - x) ** 2 + (cy - y) ** 2 + (cz - z) ** 2)
            radius = math.max(radius, dist)
        else:
            print("Non amino acid input as interface residue")
    return radius


def closest_res_distance(interface_residues: List[Residue], target_res: Residue):
    """
    Return the distance in Angstroms from target_res to the closest residue to target_res
    in interface residues
    """
    min_distance = np.inf

    # Finding the closest interface residue to the residue of interest
    for i_res in interface_residues:
        try:
            min_distance = min(min_distance, i_res["CA"] - target_res["CA"])
        except:
            print("residue didn't have alpha carbon")

    if min_distance == np.inf:
        breakpoint()
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
    wt_atoms = get_ca_atoms(wt_residues)
    mut_atoms = get_ca_atoms(mut_residues)
    if len(wt_atoms) != len(mut_atoms):
        return np.nan

    supermimposer.set_atoms(wt_atoms, mut_atoms)
    supermimposer.apply(mut_atoms)
    return supermimposer.rms


def main():
    print("NA")


def at_near_away(
    interface_residues: List[Residue], flanking_residue: Residue, cir_distance: float
):
    """
    Categorizes an insertion location as at, near, or away from an interface
    flanking residue: a residue adjacent to the insertion location
    interface_residues: The residues determined to be in the interface by change in SASA
    cir_distance: The distance from the flanking residue to the closest interface residue
    """
    if flanking_residue in interface_residues: 
        return "at"
    elif cir_distance < 5:
        return "near"
    else:
        return "away"


if __name__ == "__main__":
    main()
