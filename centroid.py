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

def face_plane(residues: List[Residue]):
    # generate random points on the plane and add random displacement
    points = np.array([res["CA"].get_coord() for res in residues]).T

    # subtract out the centroid and take the SVD
    svd = np.linalg.svd(points - np.mean(points, axis=1, keepdims=True))

    # Extract the left singular vectors
    left = svd[0]
    # the normal vector to the plane
    normal = left[:, -1]
    breakpoint()
    return normal

def plane_angle(normal1: np.ndarray, normal2: np.ndarray) -> float:
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
            dist = math.sqrt((cx-x) ** 2 + (cy-y) ** 2 + (cz-z) **2)
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

def main():
    face_plane(None)

if __name__ == "__main__":
    main()
