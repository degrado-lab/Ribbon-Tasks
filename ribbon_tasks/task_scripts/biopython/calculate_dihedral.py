##########################################
### For use in the Biopython Container ###
##########################################

import argparse
import sys
import math
from Bio import PDB
import numpy as np


def calculate_dihedral(pdb_file,
                       res1_id, chain1_id, atom1_name,
                       res2_id, chain2_id, atom2_name,
                       res3_id, chain3_id, atom3_name,
                       res4_id, chain4_id, atom4_name):
    """
    Calculate the dihedral (torsion) angle in degrees for four atoms: atom1-atom2-atom3-atom4.
    Uses vector math similar to standard molecular dihedral calculation.
    """
    # Choose parser based on file extension
    if pdb_file.endswith('.cif'):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)

    structure = parser.get_structure("structure", pdb_file)

    chain1 = structure[0][chain1_id]
    chain2 = structure[0][chain2_id]
    chain3 = structure[0][chain3_id]
    chain4 = structure[0][chain4_id]

    res1 = get_residue(chain1, res1_id)
    res2 = get_residue(chain2, res2_id)
    res3 = get_residue(chain3, res3_id)
    res4 = get_residue(chain4, res4_id)

    atom1 = res1[atom1_name]
    atom2 = res2[atom2_name]
    atom3 = res3[atom3_name]
    atom4 = res4[atom4_name]

    p0 = atom1.get_coord()
    p1 = atom2.get_coord()
    p2 = atom3.get_coord()
    p3 = atom4.get_coord()

    # Vectors between points
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # Normalize b1 for projection
    b1_norm = b1 / np.linalg.norm(b1)

    # Compute normals
    v = b0 - np.dot(b0, b1_norm) * b1_norm
    w = b2 - np.dot(b2, b1_norm) * b1_norm

    # Compute angle between v and w
    x = np.dot(v, w)
    y = np.dot(np.cross(b1_norm, v), w)

    angle_rad = math.atan2(y, x)
    angle_deg = math.degrees(angle_rad)

    return angle_deg


def get_residue(chain, res_id):
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
        description="Calculate the dihedral torsion angle for four atoms (atom1-atom2-atom3-atom4) in a PDB/mmCIF file."
    )
    parser.add_argument("pdb_file", type=str, help="Path to the PDB or mmCIF file.")
    parser.add_argument("atom1", type=str, help="First atom 'Chain:Residue:Atom' (e.g., 'A:1:N').")
    parser.add_argument("atom2", type=str, help="Second atom 'Chain:Residue:Atom' (e.g., 'A:2:CA').")
    parser.add_argument("atom3", type=str, help="Third atom 'Chain:Residue:Atom' (e.g., 'A:2:C').")
    parser.add_argument("atom4", type=str, help="Fourth atom 'Chain:Residue:Atom' (e.g., 'A:3:N').")
    parser.add_argument("--output_file", type=str, default=None, help="Optional output file to write the angle in degrees.")

    args = parser.parse_args()

    try:
        chain1_id, res1_id, atom1_name = args.atom1.split(":")
        chain2_id, res2_id, atom2_name = args.atom2.split(":")
        chain3_id, res3_id, atom3_name = args.atom3.split(":")
        chain4_id, res4_id, atom4_name = args.atom4.split(":")
    except ValueError:
        print("Error: Each atom specification must be in the format Chain:Residue:Atom (e.g. 'A:10:CA').")
        sys.exit(1)

    try:
        angle = calculate_dihedral(
            args.pdb_file,
            res1_id, chain1_id, atom1_name,
            res2_id, chain2_id, atom2_name,
            res3_id, chain3_id, atom3_name,
            res4_id, chain4_id, atom4_name,
        )
        print(f"The dihedral torsion angle between {atom1_name}-{atom2_name}-{atom3_name}-{atom4_name} is: {angle:.2f} degrees")
        if args.output_file:
            with open(args.output_file, "w") as f:
                f.write(f"{angle:.2f}\n")
    except KeyError as e:
        print(f"Error: Residue, chain, or atom not found. {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
