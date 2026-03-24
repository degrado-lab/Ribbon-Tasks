##########################################
### For use in the Reduce Container ###
##########################################

import argparse
import sys
import subprocess
from pathlib import  Path
import shim
from tempfile import NamedTemporaryFile
import shutil

default_reduce_db = "/opt/conda/envs/app/reduce_wwPDB_het_dict.txt"

def create_custom_reduce_db(molecules, outfile, infile=default_reduce_db):
    """
    Create a custom Reduce database from the provided molecules.
    
    Arguments:
    - molecules: A list of tuples containing (residue_name, sdf_file) for each molecule.
    - outfile: The path to the output database file.
    - infile: The path to the input database file (default is stored in default_reduce_db).
    """
    # Copy infile to outfile to start:
    with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
        f_out.write(f_in.read())

    for resname, sdf_file in molecules:

        std_molecule = shim.StandardMolecule(sdf_file = sdf_file)
        temp_db = NamedTemporaryFile(suffix='.txt', delete=False).name
        shim.convert.write_reduce_db(
            standard_mol=std_molecule,
            residue_name=resname,
            chemical_name=resname,
            output_path=temp_db,
        )

        # Read the current outfile content
        with open(outfile, 'r') as f_out:
            existing_content = f_out.read()
        
        # Write temp_db content first, then existing content
        # This ensures the newest entry (custom) is used first, overwriting existing entries of the same resname.
        with open(temp_db, 'r') as f_temp, open(outfile, 'w') as f_out:
            f_out.write(f_temp.read())
            f_out.write('\n')
            f_out.write(existing_content)   

def standardize_structure(structure_file: Path, output_file: Path, molecules, use_hydrogens=False):
    """
    Standardize the structure by renaming ligand atoms and chains to match the reference ligand.

    Arguments:
    - structure_file: The path to the input structure file (PDB or CIF format).
    - output_file: The path to the output standardized structure file.
    - molecules: A list of tuples containing (residue_name, sdf_file) for each molecule.
    - use_hydrogens: Whether to use hydrogens in the standardization process (default is True).
    - rename_from: A list of residue names to rename (default is None).
    """

    #print(structure_file, output_file)

    infile = structure_file
    for resname, sdf_file in molecules:
        
        print(resname, sdf_file)
        print(infile, output_file)
        print()
        shim.fix.fix_atom_names_in_residue(
            infile=infile,
            outfile=output_file,
            resname=resname,
            sdf_file=sdf_file,
            use_hydrogens=use_hydrogens,
        )
        infile = output_file
 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add Hydrogens to a PDB file. Create a custom DB if required."
    )
    parser.add_argument("input_pdb", type=str,
                        help="Path to the input PDB file.")
    parser.add_argument("output_pdb", type=str,
                        help="Path to the output PDB file.")
    parser.add_argument("--flip", action="store_true",
                        help="Flip the orientation of the added hydrogens.")
    parser.add_argument("--custom_ligand_sdfs", type=str, nargs="*",
                        help="List of custom ligand SDF files.")
    parser.add_argument("--custom_ligand_resnames", type=str, nargs="*",
                        help="List of custom ligand residue names.")
    parser.add_argument("--extra_args", type=str,
                        help="Additional arguments for the hydrogen addition.")

    args = parser.parse_args()
    if args.custom_ligand_sdfs is None:
        args.custom_ligand_sdfs = []
    if args.custom_ligand_resnames is None:
        args.custom_ligand_resnames = []
    # Zip the resnames and sdfs together:
    assert len(args.custom_ligand_sdfs) == len(args.custom_ligand_resnames), \
        "Custom ligand SDFs and residue names must have the same length. How did you even do this?"
    molecules = list(zip(args.custom_ligand_resnames, args.custom_ligand_sdfs))

    # Create intermediate pdb file:
    intermediate_pdb = NamedTemporaryFile(suffix='.pdb', delete=False).name
    
    # Create the Custom DB if needed. We only pass this in if we have molecules to add.
    with NamedTemporaryFile(suffix='.txt') as custom_db:

        # Do we need to create a custom DB for the ligands?
        if len(molecules) > 0:
            # Create DB tempfile:
            create_custom_reduce_db(molecules, custom_db.name)

            # Standardize the input file, so it matches the custom DB
            standardize_structure(args.input_pdb, intermediate_pdb, molecules)
        else:
            shutil.copy(args.input_pdb, intermediate_pdb)

        print("Running hydrogen addition...")
        flip_flag = "-FLIP" if args.flip else ""
        custom_DB = "-DB "+custom_db.name if len(molecules) > 0 else ""
        extra_args = "" if args.extra_args is None else args.extra_args
        # Run the command in the shell:
        command = f"reduce {flip_flag} {custom_DB} {extra_args} {intermediate_pdb} > {args.output_pdb}"
        print(f"Running command: {command}")
        subprocess.run(command, shell=True, check=True)
