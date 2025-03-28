##########################################
### For use in the Biopython Container ###
##########################################

import argparse
import sys
from Bio import PDB
import scipy

def calculate_atom_distance(pdb_file, res1_id, chain1_id, atom1_name, 
                            res2_id, chain2_id, atom2_name):
    # Determine if structure is a mmCIF or PDB file
    if pdb_file.endswith('.cif'):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)
    
    # Parse the structure
    structure = parser.get_structure("structure", pdb_file)
    
    # Access chains
    chain1 = structure[0][chain1_id]
    chain2 = structure[0][chain2_id]
    
    # Access residues (including hetero-residues)
    res1 = get_residue(chain1, res1_id)
    res2 = get_residue(chain2, res2_id)
    
    # Access atoms
    atom1 = res1[atom1_name]
    atom2 = res2[atom2_name]
    
    # Calculate and return the distance
    distance = atom1 - atom2
    return distance

def get_residue(chain, res_id):
    """
    Attempt to retrieve the residue from the chain. 
    First, try a standard residue (hetfield = ' ').
    If that fails, search for a hetero-residue.
    """
    try:
        # Attempt a standard residue
        residue = chain[(' ', int(res_id), ' ')]
    except KeyError:
        # If not found, look among hetero-residues
        hetero_residues = [res for res in chain.get_list() if res.id[0] != ' ']
        for res in hetero_residues:
            if res.id[1] == int(res_id):
                return res
        # If still not found, raise an error
        raise KeyError(f"Residue {res_id} not found in chain {chain.id}")
    return residue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the distance between two atoms in a PDB/mmCIF file."
    )
    parser.add_argument("pdb_file", type=str,
                        help="Path to the PDB or mmCIF file.")
    parser.add_argument("atom1", type=str,
                        help="Specification of the first atom in the format 'Chain:Residue:Atom' (e.g., 'B:99:CG').")
    parser.add_argument("atom2", type=str,
                        help="Specification of the second atom in the format 'Chain:Residue:Atom' (e.g., 'C:1:C10_1').")
    parser.add_argument("--output_file", type=str,
                        help="Path to the output file (optional).",
                        default=None)
    
    args = parser.parse_args()
    
    # Split the atom specs into chain, residue, and atom name
    try:
        chain1_id, res1_id, atom1_name = args.atom1.split(":")
        chain2_id, res2_id, atom2_name = args.atom2.split(":")
    except ValueError:
        print("Error: Each atom specification must contain exactly two colons, e.g. 'B:99:CG'.")
        sys.exit(1)
    
    # Calculate the distance
    try:
        distance = calculate_atom_distance(args.pdb_file,
                                           res1_id, chain1_id, atom1_name,
                                           res2_id, chain2_id, atom2_name)
        print(
            f"The distance between {atom1_name} in residue {res1_id} (chain {chain1_id}) "
            f"and {atom2_name} in residue {res2_id} (chain {chain2_id}) is: {distance:.2f} Å"
        )
        
        # Write output to file if requested
        if args.output_file:
            with open(args.output_file, "w") as f:
                f.write(f"{distance:.2f}\n")
    except KeyError as e:
        print(f"Error: Residue, chain, or atom not found. {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
