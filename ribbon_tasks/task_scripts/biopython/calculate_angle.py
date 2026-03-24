##########################################
### For use in the Biopython Container ###
##########################################

import argparse
import sys
import math
from Bio import PDB
import numpy as np


def calculate_atom_angle(pdb_file,
                         res1_id, chain1_id, atom1_name,
                         res2_id, chain2_id, atom2_name,
                         res3_id, chain3_id, atom3_name):
    """
    Calculate the angle (in degrees) formed at atom2 by atom1-atom2-atom3.
    """
    # Determine if structure is a mmCIF or PDB file
    if pdb_file.endswith('.cif'):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)

    structure = parser.get_structure("structure", pdb_file)

    # Retrieve chains
    chain1 = structure[0][chain1_id]
    chain2 = structure[0][chain2_id]
    chain3 = structure[0][chain3_id]

    # Retrieve residues (handles hetero-residues similarly to other scripts)
    res1 = get_residue(chain1, res1_id)
    res2 = get_residue(chain2, res2_id)
    res3 = get_residue(chain3, res3_id)

    # Retrieve atoms
    atom1 = res1[atom1_name]
    atom2 = res2[atom2_name]
    atom3 = res3[atom3_name]

    # Coordinates as numpy arrays
    p1 = atom1.get_coord()
    p2 = atom2.get_coord()
    p3 = atom3.get_coord()

    # Vectors pointing away from the central atom (atom2)
    v1 = p1 - p2
    v2 = p3 - p2

    # Compute angle between v1 and v2
    dot = np.dot(v1, v2)
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 == 0.0 or norm2 == 0.0:
        raise ValueError("One of the vectors has zero length; cannot compute angle.")

    cos_theta = dot / (norm1 * norm2)
    # Numerical safety
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    angle_rad = math.acos(cos_theta)
    angle_deg = math.degrees(angle_rad)

    return angle_deg


def get_residue(chain, res_id):
    """
    Attempt to retrieve the residue from the chain.
    First, try a standard residue (hetfield = ' ').
    If that fails, search for a hetero-residue with the same residue number.
    """
    try:
        residue = chain[(' ', int(res_id), ' ')]
    except KeyError:
        hetero_residues = [res for res in chain.get_list() if res.id[0] != ' ']
        for res in hetero_residues:
            if res.id[1] == int(res_id):
                return res
        raise KeyError(f"Residue {res_id} not found in chain {chain.id}.")
    return residue


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the angle formed by three atoms (atom1-atom2-atom3) in a PDB/mmCIF file."
    )
    parser.add_argument("pdb_file", type=str, help="Path to the PDB or mmCIF file.")
    parser.add_argument("atom1", type=str, help="First atom in format 'Chain:Residue:Atom' (e.g., 'A:10:CA').")
    parser.add_argument("atom2", type=str, help="Central atom in format 'Chain:Residue:Atom' (e.g., 'A:11:N').")
    parser.add_argument("atom3", type=str, help="Third atom in format 'Chain:Residue:Atom' (e.g., 'A:12:CB').")
    parser.add_argument("--output_file", type=str, default=None, help="Optional output file to write the angle in degrees.")

    args = parser.parse_args()

    try:
        chain1_id, res1_id, atom1_name = args.atom1.split(":")
        chain2_id, res2_id, atom2_name = args.atom2.split(":")
        chain3_id, res3_id, atom3_name = args.atom3.split(":")
    except ValueError:
        print("Error: Each atom specification must be in the format Chain:Residue:Atom (e.g. 'A:10:CA').")
        sys.exit(1)

    try:
        angle = calculate_atom_angle(
            args.pdb_file,
            res1_id, chain1_id, atom1_name,
            res2_id, chain2_id, atom2_name,
            res3_id, chain3_id, atom3_name,
        )
        print(f"The angle between {atom1_name}-{atom2_name}-{atom3_name} is: {angle:.2f} degrees")
        if args.output_file:
            with open(args.output_file, "w") as f:
                f.write(f"{angle:.2f}\n")
    except KeyError as e:
        print(f"Error: Residue, chain, or atom not found. {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
