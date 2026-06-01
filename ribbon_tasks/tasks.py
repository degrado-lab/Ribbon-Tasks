from ribbon.utils import make_directories, make_directory, list_files
from pathlib import Path
from ribbon.runner import Task
from typing import List, Union, Optional, Any, Tuple
import tempfile
import shutil
import json

RESOURCES_DIR = Path(__file__).parent / "resources"

class LigandMPNN(Task):
    def __init__(self, output_dir: Union[str, Path], structure_list : List[Union[str, Path]], num_designs: int = 1, temperature: float =0.1, device: str = 'cpu', extra_args: str = "", dry_run: bool = False):
        """
        Initialize a LigandMPNN task.

        Args:
            output_dir (str): The directory to save the output files.
            structure_list (list): A list of pdb or cif files to use as input structures.
            num_designs (int): The number of designs to generate per input structure.
            device (str): The device to run the task on. Default is 'cpu'.
            extra_args (str): Additional arguments for the task. Default is an empty string.

        Returns:
            None
            Outputs the following directories under output_dir:
                - backbones/: The generated backbones
                - packed/: The packed structures, including sidechains
                - sequences/: The generated sequences as FASTA files
                - seqs_split/: The sequences split into separate FASTA files, one per design
        """
        # Initialize the Task class
        super().__init__(device=device, extra_args=extra_args)

        # This Task name matches the name in the tasks.json file
        self.task_name = "LigandMPNN"
        
        # Your arguments here:
        self.output_dir = output_dir
        self.structure_list = structure_list
        self.num_designs = num_designs
        self.temperature = temperature
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        seqs_dir = self.output_dir / 'seqs'
        seqs_dir.mkdir(parents=True, exist_ok=True)
        
        for pdb_path in self.structure_list:
            pdb_stem = Path(pdb_path).stem
            fasta_path = seqs_dir / f"{pdb_stem}.fasta"
            with open(fasta_path, 'w') as f:
                f.write(">original\nAAAAA\n")
                for i in range(self.num_designs):
                    f.write(f">design_{i}, score=0.1, temp=0.1\nAAAAA\n")
        
        split_dir = self.output_dir / 'seqs_split'
        split_dir.mkdir(parents=True, exist_ok=True)
        for fasta_file in seqs_dir.iterdir():
            with open(fasta_file) as f:
                lines = f.readlines()[2:]
                for i in range(0, len(lines), 2):
                    if not lines[i].startswith('>'):
                        continue
                    name = lines[i].strip()
                    chains = lines[i + 1].strip().split(':')
                    index = 0
                    output_path = split_dir / f'{fasta_file.stem}_{index}.fasta'
                    while output_path.exists():
                        index += 1
                        output_path = split_dir / f'{fasta_file.stem}_{index}.fasta'
                    with open(output_path, 'w') as g:
                        for chain_index, chain in enumerate(chains):
                            name_with_chain = f"{name.split(',')[0]}_{chain_index}" + ', ' + ','.join(name.split(',')[1:])
                            g.write(f'{name_with_chain}\n{chain}\n')

    def run(self):
        if self.dry_run:
            self._run_dry()
            return
        
        ###### HELPER FUNCTIONS #######
        def split_ligandmpnn_fasta(fasta_file, split_output_dir):
            # Each LigandMPNN input produces a FASTA with multiple outputs as > lines. Each line has 1 or more chains separated by ':'.
            # Here, we separate each output into it's own FASTA file, with chains as separate > lines.
            split_output_dir.mkdir(parents=True, exist_ok=True)
            
            with open(fasta_file) as f:
                # Skip the first two lines (original input)
                lines = f.readlines()[2:]
                
                for i in range(0, len(lines), 2):
                    if not lines[i].startswith('>'):
                        continue
                    
                    name = lines[i].strip()
                    chains = lines[i + 1].strip().split(':')
                    index = 0
                    output_path = split_output_dir / f'{fasta_file.stem}_{index}.fasta'
                    
                    # Increment the index to avoid overwriting existing files
                    while output_path.exists():
                        index += 1
                        output_path = split_output_dir / f'{fasta_file.stem}_{index}.fasta'
                    
                    print(f'Writing to {output_path}')
                    with open(output_path, 'w') as g:
                        for chain_index, chain in enumerate(chains):
                            name_with_chain = f"{name.split(',')[0]}_{chain_index}" +', ' + ','.join(name.split(',')[1:])# Add chain to name
                            g.write(f'{name_with_chain}\n{chain}\n')

        # Make directories:
        self.output_dir = make_directory(self.output_dir)

        # Then, write out the files within pdb_input_dir to a json file:
        #pdb_input_json = self.output_dir / 'pdb_input.json'
        # Make a temp file for the json:
        pdb_input_json = tempfile.NamedTemporaryFile(delete=False).name
        with open(pdb_input_json, 'w') as f:
            json.dump(self.structure_list, f)
        
        # Run the task:
        self._run_task(self.task_name, 
                    pdb_input_json = pdb_input_json, 
                    output_dir = self.output_dir, 
                    num_designs = self.num_designs,
                    temperature = self.temperature,
                    extra_args = self.extra_args,
                    device = self.device)
        
        # Split the FASTA files:
        for file in (self.output_dir / 'seqs').iterdir():
            print(f'Splitting {file}')
            split_ligandmpnn_fasta(file, self.output_dir / 'seqs_split')
        
        return 

class LASErMPNN(Task):
    def __init__(self, output_dir: Union[str, Path], structure_list: List[Union[str, Path]], num_designs: int = 1, temperature: float = 0.000001, device: str = 'cpu', fix_beta: bool = False, extra_args: str = "", dry_run: bool = False):
        """
        Initialize a LASErMPNN task.

        Args:
            output_dir (str): The directory to save the output files.
            structure_list (list): A list of pdb or cif files to use as input structures.
            num_designs (int): The number of designs to generate per input structure.
            device (str): The device to run the task on. Default is 'cpu'.
            extra_args (str): Additional arguments for the task. Default is an empty string.

        TODO: Add feature which sets selected residues to have a b-factor of 1.0, and then passes the flag --fix_beta so we can fix those residues.
                For now, the user manually sets the b-factors in the input PDB files.

        Returns:
            None
            Outputs the following directories under output_dir:
                - [input_pdb_stem]/: Directory containing designs for each input structure
                    - design_0.pdb, design_1.pdb, ...: Generated design structures
        """
        # Initialize the Task class
        super().__init__(device=device, extra_args=extra_args)

        # This Task name matches the name in the tasks.json file
        self.task_name = "LASErMPNN"
        
        # Your arguments here:
        self.output_dir = output_dir
        self.structure_list = structure_list
        self.num_designs = num_designs
        self.temperature = temperature
        self.fix_beta = fix_beta
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        import shutil
        for pdb_file in self.structure_list:
            pdb_stem = Path(pdb_file).stem
            pdb_out_dir = self.output_dir / pdb_stem
            pdb_out_dir.mkdir(parents=True, exist_ok=True)
            for i in range(self.num_designs):
                shutil.copy(pdb_file, pdb_out_dir / f"design_{i}.pdb")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return
        
        ###### HELPER FUNCTIONS #######
        
        # Make directories:
        self.output_dir = make_directory(self.output_dir)

        # LaserMPNN expects a single input folder. We'll make a temporary folder with symlinks to our PDBs:
        pdb_input_dir = tempfile.mkdtemp()
        import os
        for pdb_file in self.structure_list:
            os.symlink(Path(pdb_file).resolve(), os.path.join(pdb_input_dir, os.path.basename(pdb_file)))

        # Run the task:
        self._run_task(self.task_name, 
                    input_dir = str(pdb_input_dir),
                    output_dir = self.output_dir, 
                    designs_per_pdb = self.num_designs,
                    extra_args = self.extra_args,
                    temperature = self.temperature,
                    fix_beta = "--fix_beta" if self.fix_beta else "",
                    device = self.device)
        
        return 

class FastRelax(Task):
    def __init__(self, output_dir: Union[str, Path], pdb_input_file: Optional[Union[str, Path]] = None, pdb_input_dir: Optional[Union[str, Path]] = None, nstructs: int = 1, ligand_params_file: Optional[Union[str, Path, List[Union[str, Path]]]] = None, custom_bonds: Optional[List[str]] = None, custom_angles: Optional[List[str]] = None, custom_torsions: Optional[List[str]] = None, constraints_weight: Optional[float] = None, device: str = 'cpu', extra_args: str = "", dry_run: bool = False):
        """
        Initialize a FastRelax task.

        Args:
            output_dir (str): The directory to save the output files.
            pdb_input_file (str, optional): Path to a single PDB file. Default is None.
            pdb_input_dir (str, optional): Path to a directory containing PDB files. Default is None.
            nstructs (int): The number of structures to generate. Default is 1.
            ligand_params_file (str or list of str, optional): Path(s) to ligand .params files. Default is None.
            custom_bonds (list of str, optional): Custom bond/distance constraints.
                Formats:
                  - Default (HARMONIC): "Atom1,Atom2,Target,SD"
                    E.g. ["A:16:ND1,B:1:H1,2.0,0.5"]
                  - Custom function and parameters: "Atom1,Atom2,FUNC_NAME,Param1,Param2,..."
                    E.g. ["A:45:OD2,A:16:HE2,FLAT_HARMONIC,0.0,0.5,6"]
            custom_angles (list of str, optional): Custom angle constraints.
                Formats:
                  - Default (HARMONIC): "Atom1,Atom2,Atom3,Target,SD"
                    E.g. ["A:16:CE1,A:16:ND1,B:1:H1,2.41,0.5"]
                  - Custom function and parameters: "Atom1,Atom2,Atom3,FUNC_NAME,Param1,Param2,..."
                    E.g. ["A:1:CA,A:2:CA,A:3:CA,BOUNDED,1.0,2.0,0.1,TAG"]
            custom_torsions (list of str, optional): Custom dihedral/torsion constraints.
                Formats:
                  - Default (CIRCULARHARMONIC): "Atom1,Atom2,Atom3,Atom4,Target,SD"
                    E.g. ["A:45:OD1,A:45:CG,A:45:OD2,A:16:HE2,3.14,0.35"]
                  - Custom function and parameters: "Atom1,Atom2,Atom3,Atom4,FUNC_NAME,Param1,Param2,..."
            constraints_weight (float, optional): The weight for full-atom constraints. Default is None.
            device (str): The device to run the task on. Default is 'cpu'.
            extra_args (str): Additional arguments for the task. Default is an empty string.

        Raises:
            ValueError: If neither pdb_input_file nor pdb_input_dir is specified.

        Returns:
            None
            Outputs the following files under output_dir:
                - [input_pdb_stem]_0001.pdb, etc.: The relaxed structure files
                - score.sc: Rosetta score file containing energy terms
        """
         # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "FastRelax"
        
        # Task-specific variables
        if pdb_input_file is None and pdb_input_dir is None:
            raise ValueError('Must specify either pdb_input_file or pdb_input_dir')
        self.output_dir = make_directory(output_dir)
        self.pdb_input_file = pdb_input_file
        self.pdb_input_dir = pdb_input_dir
        self.nstructs = nstructs
        self.ligand_params_file = ligand_params_file
        self.custom_bonds = custom_bonds or []
        self.custom_angles = custom_angles or []
        self.custom_torsions = custom_torsions or []
        self.constraints_weight = constraints_weight
        self.device = device
        self.extra_args = extra_args
        self.dry_run = dry_run

    def _parse_atom_selector(self, atom_str: str) -> str:
        """
        Parses atom selector string in format 'Chain:Res:Atom' into Rosetta's PDB constraint format.
        E.g. 'A:16:ND1' -> 'ND1 16A'
        If chain is empty, it formats as pose/continuous numbering (e.g. ':16:ND1' -> 'ND1 16').
        """
        parts = atom_str.split(':')
        if len(parts) != 3:
            raise ValueError(f"Invalid atom selector: '{atom_str}'. Expected format 'Chain:Res:Atom'")
        chain, res, atom = parts
        return f"{atom} {res}{chain}"

    def _parse_constraint(self, constraint_str: str, num_atoms: int, default_func: str) -> str:
        """
        Parses a generic constraint string representing biological constraints.
        Extracts num_atoms first, and then parses the remainder as function name (optional)
        followed by its parameters.
        """
        parts = [p.strip() for p in constraint_str.split(',')]
        if len(parts) < num_atoms + 1:
            raise ValueError(f"Invalid constraint: '{constraint_str}'. Expected at least {num_atoms} atoms and 1 parameter.")
        
        atoms = [self._parse_atom_selector(parts[i]) for i in range(num_atoms)]
        remaining = parts[num_atoms:]
        
        try:
            # Check if first remaining parameter can be parsed as a number
            float(remaining[0])
            func = default_func
            params = remaining
        except ValueError:
            # First parameter is not a number, so it must be the custom function name
            func = remaining[0]
            params = remaining[1:]
            if not params:
                raise ValueError(f"Invalid constraint: '{constraint_str}'. Expected parameters after function '{func}'.")
        
        atoms_str = " ".join(atoms)
        params_str = " ".join(params)
        return f"{atoms_str} {func} {params_str}"

    def _generate_constraints_file(self) -> Optional[Path]:
        """
        Generates and writes a Rosetta constraint file containing custom bonds, angles, and torsions.
        Returns the Path to the written file, or None if no custom constraints were specified.
        """
        cst_lines = []

        # Parse custom bonds
        for bond in self.custom_bonds:
            cst_lines.append(f"AtomPair {self._parse_constraint(bond, 2, 'HARMONIC')}")

        # Parse custom angles
        for angle in self.custom_angles:
            cst_lines.append(f"Angle {self._parse_constraint(angle, 3, 'HARMONIC')}")

        # Parse custom torsions
        for torsion in self.custom_torsions:
            cst_lines.append(f"Dihedral {self._parse_constraint(torsion, 4, 'CIRCULARHARMONIC')}")

        if not cst_lines:
            return None

        cst_file = self.output_dir / "constraints.cst"
        cst_file.write_text("\n".join(cst_lines) + "\n")
        return cst_file

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        if self.pdb_input_file is not None:
            pdb_list = [Path(self.pdb_input_file)]
        else:
            pdb_list = list_files(self.pdb_input_dir, '.pdb')
        
        import shutil
        for pdb_path in pdb_list:
            pdb_stem = Path(pdb_path).stem
            for i in range(1, self.nstructs + 1):
                shutil.copy(pdb_path, self.output_dir / f"{pdb_stem}_{i:04d}.pdb")
        
        with open(self.output_dir / "score.sc", "w") as f:
            f.write("SEQUENCE:\nSCORE:     score description\n")
            for pdb_path in pdb_list:
                pdb_stem = Path(pdb_path).stem
                for i in range(1, self.nstructs + 1):
                    f.write(f"SCORE:     -100.0 {pdb_stem}_{i:04d}\n")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Handle input files
        if self.pdb_input_file is not None:
            temp_dir = tempfile.mkdtemp()
            shutil.copy(self.pdb_input_file, temp_dir)
            pdb_input_dir = temp_dir
        else:
            pdb_input_dir = self.pdb_input_dir

        # Prepare PDB files
        pdb_list = list_files(pdb_input_dir, '.pdb')
        pdb_string = " ".join(map(str, pdb_list)) + " "

        # Dynamic arguments lists
        params_flags = []

        # Generate constraints file
        cst_file = self._generate_constraints_file()
        if cst_file:
            params_flags.append(f"-constraints:cst_fa_file {cst_file.resolve()}")

        # Add constraints weight
        if self.constraints_weight is not None:
            params_flags.append(f"-constraints:cst_fa_weight {self.constraints_weight}")

        # Add ligand params files
        if self.ligand_params_file:
            if isinstance(self.ligand_params_file, (str, Path)):
                params_files = [self.ligand_params_file]
            else:
                params_files = self.ligand_params_file
            
            resolved_files = []
            for pf in params_files:
                pf_path = Path(pf).resolve()
                if not pf_path.exists():
                    raise FileNotFoundError(f"Ligand params file not found: {pf}")
                resolved_files.append(str(pf_path))
            
            if resolved_files:
                params_flags.append(f"-extra_res_fa {' '.join(resolved_files)}")

        # Combine dynamic flags into extra_args
        dynamic_args = " ".join(params_flags)
        combined_extra_args = f"{dynamic_args} {self.extra_args}".strip() if self.extra_args else dynamic_args

        # Run the task
        self._run_task(
            self.task_name,
            pdb_string=pdb_string,
            output_dir=str(self.output_dir),
            nstructs=self.nstructs,
            extra_args=combined_extra_args,
            device=self.device
        )

class Chai1(Task):
    def __init__(self, fasta_file: Union[str, Path], output_dir: Union[str, Path] = '.', smiles_string: Optional[str] = None, num_ligands: int = 1, device: str = 'gpu', dry_run: bool = False):
        """
        Initialize a Chai-1 task.

        Args:
            fasta_file (str): The FASTA file containing the protein sequence (no ligand).
            output_dir (str): The directory to save the output files. Default is '.'.
            smiles_string (str, optional): The SMILES string of the ligand. Default is None.
            num_ligands (int): The number of ligands. Default is 1.
            device (str): The device to run the task on. Default is 'gpu'.

        Returns:
            None
            Outputs the following files under output_dir:
                - [input_fasta_stem]_idx_0.cif, etc.: The predicted 3D structures in CIF format
                - [input_fasta_stem]_idx_0.npz, etc.: Score files containing pTM, ipTM, pLDDT, etc.
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Chai-1"

        # Task-specific variables
        self.fasta_file = fasta_file
        self.smiles_string = smiles_string
        self.output_dir = output_dir
        self.device = device
        self.num_ligands = num_ligands
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        fasta_stem = Path(self.fasta_file).stem
        
        cif_path = self.output_dir / f"{fasta_stem}_idx_0.cif"
        shutil.copy(RESOURCES_DIR / "dummy.cif", cif_path)
        
        npz_path = self.output_dir / f"{fasta_stem}_idx_0.npz"
        import zipfile
        with zipfile.ZipFile(npz_path, "w") as z:
            z.writestr("dummy.npy", b"")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Make the directory:
        self.output_dir = make_directory(self.output_dir)

        # Run the task
        self._run_task(
            self.task_name,
            fasta_file=self.fasta_file,
            smiles_string=self.smiles_string,
            output_dir=str(self.output_dir),
            num_ligands=self.num_ligands,
            device=self.device
        )

class Boltz2(Task):
    def __init__(self, fasta_file: Union[str, Path], output_dir: Union[str, Path] = '.', smiles_list: List[str] = [], calculate_binding: bool = False, use_msa_server: bool = True, device: str = 'gpu', extra_args: str = "", dry_run: bool = False):
        """
        Initialize a Boltz-2 task.

        Args:
            fasta_file (str): The FASTA file containing the protein sequence (no ligand).
            output_dir (str): The directory to save the output files. Default is '.'.
            smiles_list (list of str, optional): List of the SMILES strings of the ligands. Repeat a SMILES string to use multiple copies of the same ligand. 
                Default is no ligands.
            calculate_binding (bool, optional): Whether to calculate binding properties of the ligand. Only calculates for the first ligand in the list. Default is False.
            device (str): The device to run the task on. Default is 'gpu'.

        Returns:
            None
            Outputs the following directories/files under output_dir:
                - boltz_results_[input_yaml_stem]/:
                    - predictions/[input_yaml_stem]/[input_yaml_stem]_model_0.cif, etc.: Predicted structure files
                    - predictions/[input_yaml_stem]/confidence_[input_yaml_stem]_model_0.json: Model confidences and metrics
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Boltz-2"

        # Task-specific variables
        self.fasta_file = fasta_file
        self.smiles_list = smiles_list
        self.output_dir = output_dir
        self.device = device
        self.calculate_binding = calculate_binding
        self.extra_args = extra_args
        self.use_msa_server = use_msa_server
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        fasta_stem = Path(self.fasta_file).stem
        
        pred_dir = self.output_dir / f"boltz_results_{fasta_stem}" / "predictions" / fasta_stem
        pred_dir.mkdir(parents=True, exist_ok=True)
        
        cif_path = pred_dir / f"{fasta_stem}_model_0.cif"
        shutil.copy(RESOURCES_DIR / "dummy.cif", cif_path)
        
        json_path = pred_dir / f"confidence_{fasta_stem}_model_0.json"
        with open(json_path, "w") as f:
            json.dump({"plddt": [90.0], "ptm": 0.8, "iptm": 0.7}, f)

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Make the directory:
        self.output_dir = make_directory(self.output_dir)

        # Only add smiles flag if necessary. 
        ligand_flags = ""
        if self.smiles_list:
            smiles_string = "--ligand_smiles " + " ".join(f'"{smiles}"' for smiles in self.smiles_list)
            ligand_flags += " " + smiles_string

        if self.calculate_binding:
            ligand_flags += " --calculate_binding"

        if self.device == 'cpu':
            self.extra_args += " --accelerator cpu "

        if self.use_msa_server:
            self.use_msa_server = " --use_msa_server "
        else:
            self.use_msa_server = " "

        # Run the task
        self._run_task(
            self.task_name,
            fasta_file=self.fasta_file,
            output_dir=str(self.output_dir),
            use_msa_server=self.use_msa_server,
            device=self.device,
            ligand_flags=ligand_flags,
            extra_args=self.extra_args
        )

class RaptorXSingle(Task):
    def __init__(self, fasta_file_or_dir: Union[str, Path], output_dir: Union[str, Path] = '.', param: str = 'RaptorX-Single-ESM1b.pt', device: str = 'gpu', extra_args: str = "", dry_run: bool = False):
        """
        Initialize a RaptorXSingle task.

        Args:
            fasta_file_or_dir (str): The FASTA file or directory containing multiple FASTA files.
            output_dir (str): The directory to save the output files. Default is '.'.
            param (str): The checkpoint to use. Default is 'RaptorX-Single-ESM1b.pt'.
            device (str): The device to run the task on. Default is 'gpu'.
            extra_args (str): Additional arguments for the task. Default is an empty string.

        Raises:
            ValueError: If an invalid param is specified.

        Returns:
            None
            Outputs the following files under output_dir:
                - [input_fasta_stem].pdb: Predicted structure files
        """
        
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "RaptorXSingle"

        # Task-specific variables
        self.fasta_file = fasta_file_or_dir
        self.output_dir = output_dir
        self.param = 'RaptorX-Single/params/'+param  #This is the directory where the params are stored
        self.extra_args = extra_args
        self.device = device
        self.dry_run = dry_run

        if device == 'gpu':
            self.device_id = '0'
        elif device == 'cpu':
            self.device_id = '-1'
        else:
            self.device_id = str(device)
            self.device = 'cpu' if device == '-1' else 'gpu'

        # Check inputs:
        valid_param_list = [
            'RaptorX-Single-ESM1b.pt',
            'RaptorX-Single-ESM1v.pt',
            'RaptorX-Single-ProtTrans.pt',
            'RaptorX-Single-ESM1b-ESM1v-ProtTrans.pt',
            'RaptorX-Single-ESM1b-Ab.pt',
            'RaptorX-Single-ESM1v-Ab.pt',
            'RaptorX-Single-ProtTrans-Ab.pt',
            'RaptorX-Single-ESM1b-ESM1v-ProtTrans-Ab.pt'
        ]

        if param not in valid_param_list:
            raise ValueError(f'Invalid param: {param}. Must be one of {valid_param_list}')

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        fasta_stem = Path(self.fasta_file).stem
        
        pdb_path = self.output_dir / f"{fasta_stem}.pdb"
        shutil.copy(RESOURCES_DIR / "dummy.pdb", pdb_path)

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Make the directory:
        self.output_dir = make_directory(self.output_dir)
        
        # "python RaptorX-Single/pred.py {fasta_file} RaptorX-Single/params/{param} --plm_param_dir RaptorX-Single/params/ --out_dir {output_dir} --device {device} {extra_args}",

        # Run the task
        self._run_task(
            self.task_name,
            fasta_file=self.fasta_file,
            param=self.param,
            output_dir=str(self.output_dir),
            device=self.device,
            device_id=self.device_id,
            extra_args=self.extra_args
        )

class CalculateDistance(Task):
    def __init__(self, pdb_file: Union[str, Path], atom1: str,
                 atom2: str, output_file: Union[str, Path], device: str = 'cpu', dry_run: bool = False):
        """
        Initialize a CalculateDistance task.
        This calculates the distance between two atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            atom1_name (str): Name of the first atom in the format 'Chain:Residue:Atom'.
            atom2_name (str): Name of the second atom in the format 'Chain:Residue:Atom'.
            output_file (str): Path to the output file. Suffixed with '.dist'.
            device (str): The device to run the task on. Default is 'cpu'.

        Returns:
            None
            Outputs:
                - output_file: Text file containing the calculated distance (in Angstroms), suffixed with '.dist'
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Calculate Distance"
        
        # Task-specific variables
        self.pdb_file = pdb_file
        self.atom1 = atom1
        self.atom2 = atom2
        self.output_file = output_file
        self.device = device
        self.dry_run = dry_run

    def _run_dry(self):
        make_directory(Path(self.output_file).parent)
        with open(self.output_file, 'w') as f:
            f.write("5.0\n")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Ensure output directory exists
        make_directory(Path(self.output_file).parent)

        # Run the task
        self._run_task(
            self.task_name,
            pdb_file=self.pdb_file,
            atom1=self.atom1,
            atom2=self.atom2,
            output_file=self.output_file,
            device=self.device
        )

class CalculateAngle(Task):
    def __init__(self, pdb_file: Union[str, Path], atom1: str,
                 atom2: str, atom3: str, output_file: Union[str, Path], device: str = 'cpu', dry_run: bool = False):
        """
        Initialize a CalculateAngle task.
        This calculates the angle formed by three atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            atom1_name (str): Name of the first atom in the format 'Chain:Residue:Atom'.
            atom2_name (str): Name of the second atom in the format 'Chain:Residue:Atom'.
            atom3_name (str): Name of the third atom in the format 'Chain:Residue:Atom'.
            output_file (str): Path to the output file. Suffixed with '.dist'.
            device (str): The device to run the task on. Default is 'cpu'.

        Returns:
            None
            Outputs:
                - output_file: Text file containing the calculated angle (in degrees), suffixed with '.angle'
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Calculate Angle"
        
        # Task-specific variables
        self.pdb_file = pdb_file
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.output_file = output_file
        self.device = device
        self.dry_run = dry_run

    def _run_dry(self):
        make_directory(Path(self.output_file).parent)
        with open(self.output_file, 'w') as f:
            f.write("109.5\n")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Ensure output directory exists
        make_directory(Path(self.output_file).parent)

        # Run the task
        self._run_task(
            self.task_name,
            pdb_file=self.pdb_file,
            atom1=self.atom1,
            atom2=self.atom2,
            atom3=self.atom3,
            output_file=self.output_file,
            device=self.device
        )

class CalculateDihedral(Task):
    def __init__(self, pdb_file: Union[str, Path], atom1: str,
                 atom2: str, atom3: str, atom4: str, output_file: Union[str, Path], device: str = 'cpu', dry_run: bool = False):
        """
        Initialize a CalculateDihedral task.
        This calculates the dihedral torsion angle for four atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            atom1_name (str): Name of the first atom in the format 'Chain:Residue:Atom'.
            atom2_name (str): Name of the second atom in the format 'Chain:Residue:Atom'.
            atom3_name (str): Name of the third atom in the format 'Chain:Residue:Atom'.
            atom4_name (str): Name of the fourth atom in the format 'Chain:Residue:Atom'.
            output_file (str): Path to the output file. Suffixed with '.dist'.
            device (str): The device to run the task on. Default is 'cpu'.

        Returns:
            None
            Outputs:
                - output_file: Text file containing the calculated dihedral torsion angle (in degrees), suffixed with '.dihedral'
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Calculate Dihedral"
        
        # Task-specific variables
        self.pdb_file = pdb_file
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.output_file = output_file
        self.device = device
        self.dry_run = dry_run

    def _run_dry(self):
        make_directory(Path(self.output_file).parent)
        with open(self.output_file, 'w') as f:
            f.write("180.0\n")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Ensure output directory exists
        make_directory(Path(self.output_file).parent)
        
        # Run the task
        self._run_task(
            self.task_name,
            pdb_file=self.pdb_file,
            atom1=self.atom1,
            atom2=self.atom2,
            atom3=self.atom3,
            atom4=self.atom4,
            output_file=self.output_file,
            device=self.device
        )

class CalculatePairwiseDistance(Task):
    def __init__(self, pdb_file: Union[str, Path], atom_list_A: List[str],
                 atom_list_B: List[str], output_file: Union[str, Path], average: bool = False, device: str = 'cpu', dry_run: bool = False):
        """
        Initialize a CalculateDistance task.
        This calculates the distance between two atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            atom_list_A (list): List of atoms in the first group in the format [ 'Chain:Residue:Atom', ...].
            atom_list_B (list): Name of the second atom in the format [ 'Chain:Residue:Atom', ...].
            output_file (str): Path to the output file. Must be a CSV (.csv) or JSON (.json).
                CSV has the columns: "A_index", "AtomA", "B_index", "AtomB", "Distance".
                JSON is a list of dicts with keys: "A_index", "A_spec", "B_index", "B_spec", "distance"
            average (bool): Whether to return the average distance. Default is False.
                If True, the output file will contain a single value (no other information).
            device (str): The device to run the task on. Default is 'cpu'.
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Calculate Pairwise Distance"

        # Check output file extension:
        if not average:
            if Path(output_file).suffix not in ['.csv', '.json']:
                raise ValueError('Output file must be a CSV (.csv) or JSON (.json) when outputting multiple distances.\n \
                                 If you want to output a single average distance, set average=True.')
        
        # Task-specific variables
        self.pdb_file = pdb_file
        self.atom_list_A = atom_list_A
        self.atom_list_B = atom_list_B
        self.output_file = output_file
        self.average = average
        self.device = device
        self.dry_run = dry_run

    def _run_dry(self):
        make_directory(Path(self.output_file).parent)
        if self.average:
            with open(self.output_file, 'w') as f:
                f.write("5.0\n")
        elif Path(self.output_file).suffix == '.csv':
            with open(self.output_file, 'w') as f:
                f.write("A_index,AtomA,B_index,AtomB,Distance\n1,A:1:CA,2,A:2:CA,5.0\n")
        else: #json
            with open(self.output_file, 'w') as f:
                json.dump([{"A_index": 1, "A_spec": "A:1:CA", "B_index": 2, "B_spec": "A:2:CA", "distance": 5.0}], f)

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Ensure output directory exists
        make_directory(Path(self.output_file).parent)

        atom_group_A = " ".join(self.atom_list_A)
        atom_group_B = " ".join(self.atom_list_B)

        # Add the proper flag to the output file,
        # Since it depends on the file extension
        if self.average:
            self.output_file = f"--average_output {self.output_file}"
        elif Path(self.output_file).suffix == '.csv':
            self.output_file = f"--csv_output {self.output_file}"
        else: #json
            self.output_file = f"--json_output {self.output_file}"
            

        # Run the task
        self._run_task(
            self.task_name,
            pdb_file=self.pdb_file,
            atom_group_A=atom_group_A,
            atom_group_B=atom_group_B,
            output_file=self.output_file,
            device=self.device
        )

class AddHydrogens(Task):
    def __init__(self, input_file: Union[str, Path], output_file: Union[str, Path], selection: str = 'all', dry_run: bool = False):
        """
        Initialize a CalculateDistance task.
        This calculates the distance between two atoms in a PDB file.

        Args:
            input_file (str): Path to the PDB or CIF file.
            output_file (str): Path to the output file.
            selection (str): PyMol selection string to specify what to add hydrogens to. Default is 'all'.

        Returns:
            None
            Outputs:
                - output_file: PDB or CIF structure file containing the added hydrogens
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Add Hydrogens"
        
        # Task-specific variables
        self.input_file = input_file
        self.selection = selection
        self.output_file = output_file
        self.device = 'cpu'
        self.dry_run = dry_run

    def _run_dry(self):
        make_directory(Path(self.output_file).parent)
        import shutil
        shutil.copy(self.input_file, self.output_file)

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Ensure output directory exists
        make_directory(Path(self.output_file).parent)

        # Run the task
        self._run_task(
            self.task_name,
            input_file=self.input_file,
            output_file=self.output_file,
            selection=self.selection,
            device=self.device
        )

class Reduce(Task):
    def __init__(self, pdb_input_file: Union[str, Path], pdb_output_file: Union[str, Path], flip: bool = False, custom_ligands: List[Tuple[str, Union[str, Path]]] = [], dry_run: bool = False):
        """
        Add Hydrogens to a PDB file.

        Args:
            pdb_input_file (str): Path to the input PDB file.
            pdb_output_file (str): Path to the output PDB file.
            flip (bool): Whether to optionally flip N/Q/H residues. Default False.
            custom_ligands (list of tuples): List of custom ligand resnames and SDF files. (E.g. [('KP1', 'kemp1.sdf'), ...] ) Only necessary if there is a ligand which is not already in the Protein Data Bank.

        Returns:
            None
            Outputs:
                - pdb_output_file: PDB structure file containing the added hydrogens
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Reduce"

        # Verify that custom_ligands is a list of tuples:
        if not isinstance(custom_ligands, list):
            raise ValueError("custom_ligands must be a list of tuples.")
        for item in custom_ligands:
            if not isinstance(item, tuple) or len(item) != 2:
                raise ValueError("Each item in custom_ligands must be a tuple of (resname, sdf_file).")

        # Task-specific variables
        self.pdb_input_file = pdb_input_file
        self.pdb_output_file = pdb_output_file
        self.custom_ligands = custom_ligands
        self.flip = flip
        self.device = 'cpu'
        self.dry_run = dry_run

    def _run_dry(self):
        make_directory(Path(self.pdb_output_file).parent)
        import shutil
        shutil.copy(self.pdb_input_file, self.pdb_output_file)

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Ensure output directory exists
        make_directory(Path(self.pdb_output_file).parent)

        # unzip resnames and sdfs:
        custom_resnames = [name for name, _ in self.custom_ligands]
        custom_sdfs = [sdf for _, sdf in self.custom_ligands]

        # turn into a string
        if len(custom_resnames) > 0:
            custom_resnames = "--custom_ligand_resnames " + " ".join(custom_resnames)
        if len(custom_sdfs) > 0:
            custom_sdfs = "--custom_ligand_sdfs " + " ".join([str(custom_sdf) for custom_sdf in custom_sdfs])

        # Run the task
        self._run_task(
            self.task_name,
            input_pdb=self.pdb_input_file,
            output_pdb=self.pdb_output_file,
            custom_ligand_sdfs=     "" if not custom_sdfs else custom_sdfs,
            custom_ligand_resnames= "" if not custom_resnames else custom_resnames,
            flip="--flip" if self.flip else "",
            device=self.device
        )

class CalculateSASA(Task):
    def __init__(self, pdb_file: Union[str, Path], output_file: Union[str, Path], atom_1: str, device: str = 'cpu', dry_run: bool = False):
        """
        Initialize a CalculateSASA task.
        This calculates the Solvent Accessible Surface Area for a set of atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            output_file (str): Path to the output file. Suffixed with '.angle'.
            atom_1 (str): Atom specification in format chain_id:res_id:atom_name.
            device (str): The device to run the task on. Default is 'cpu'.

        Returns:
            None
            Outputs:
                - output_file: Text file containing calculated solvent accessible surface area (SASA) values

        TODO:
            - Implement the task script in ribbon/ribbon_tasks/task_scripts/calculate_sasa.py
        """

        self.dry_run = dry_run
        raise NotImplementedError('This task is not yet implemented .')
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Calculate SASA"
        
        # Task-specific variables
        self.pdb_file = pdb_file
        self.output_file = output_file
        self.atom_1 = atom_1
        self.device = device

    def run(self):
        # Ensure output directory exists
        make_directory(Path(self.output_file).parent)

        # Run the task
        self._run_task(
            self.task_name,
            pdb_file=self.pdb_file,
            output_file=self.output_file,
            atom_1=self.atom_1,
            device=self.device
        )

class RFDiffusionAA(Task):
    def __init__(self, input_structure: Union[str, Path], output_dir: Union[str, Path], contig_map: str, num_designs: int = 1, total_length: str = 'null', ligand: str = 'null',  diffuser_steps: int = 200, deterministic: bool = False, design_startnum: int = 0, force: bool = False, device: str = 'gpu', extra_args: str = "", dry_run: bool = False):
        """
        Initialize a RFDiffusionAA task.

        Args:
            input_structure (str): The input structure file.
            output_dir (str): The directory to save the output files.
            contig_map (str): The contig map. This is formatted as a hydra list. E.g. "[150-150]".
            num_designs (int): The number of designs to generate per input structure.
            total_length (str): The total length of the protein (e.g. "150-150"). Default is 'null'.
            ligand (str): The residue ID of the ligand to dock. (e.g. "OQO"). Default is 'null'.
            diffuser_steps (int): The number of steps to run the diffuser for. Default is 200.
            deterministic (bool): Whether to use a deterministic seed for the diffuser. Default is False.
            design_startnum (int): The starting number for the design output files. Useful for batching jobs in the same out directory. Default is 0.
            force (bool): Whether to overwrite existing output files. Default is False.
            device (str): The device to run the task on. Default is 'gpu'.
            extra_args (str): Additional arguments for the task. Default is an empty string.

        Returns:
            None
            Outputs the following files under output_dir:
                - output_dir/design_0.pdb, etc.: The generated backbone structures
                - output_dir/design_0.trb, etc.: Metadata files containing diffusion traces and alignment details
        """
        # Initialize the Task class
        super().__init__(device=device, extra_args=extra_args)

        # This Task name matches the name in the tasks.json file
        self.task_name = "RFDiffusionAA"
        
        # Your arguments here:
        self.input_structure = input_structure
        self.output_dir = output_dir
        self.contig_map = contig_map
        self.num_designs = num_designs
        self.total_length = total_length
        self.ligand = ligand
        self.diffuser_steps = diffuser_steps
        self.deterministic = deterministic
        self.design_startnum = design_startnum
        self.force = force
        self.dry_run = dry_run
        # we should have a final length flag for config.length

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        import shutil
        for i in range(self.num_designs):
            shutil.copy(self.input_structure, self.output_dir / f"design_{i}.pdb")
            with open(self.output_dir / f"design_{i}.trb", "w") as f:
                f.write("")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return
        
        # Make directories:
        self.output_dir = make_directory(self.output_dir)

        # Convert the output_dir to a prefix:
        output_prefix = Path(self.output_dir) / 'design'
        
        # Run the task:
        self._run_task(self.task_name, 
                    input_structure = self.input_structure,
                    output_prefix = str(output_prefix),  # we provide a prefix (e.g. ./out/sample) instead of dir (e.g. ./out)
                    contig_map = self.contig_map,
                    num_designs = self.num_designs,
                    total_length = self.total_length,
                    ligand = self.ligand,
                    diffuser_steps = self.diffuser_steps,
                    deterministic = self.deterministic,
                    design_startnum = self.design_startnum,
                    cautious = not self.force,  # If force is True, then we overwrite existing files
                    extra_args = self.extra_args,
                    device = self.device)
        
        return 

class EasyMD(Task):
    #"command": "easymd run {input_file} --output {output_prefix} --duration {duration} --relax-duration {relax_duration} --output-frequency {output_frequency} {ligand_files} {forcefield_files} {water_model} {pH} {hydrogen_variants} {ionic_strength} {box_padding} {custom_bonds} {custom_angles} {custom_torsions} {minimize_only} {extra_args}", 
        
    def __init__(self, input_file: Union[str, Path], output_prefix: Union[str, Path], duration: int, relax_duration: int = 1, output_frequency: int = 1, ligand_files: List[Union[str, Path]] = [], forcefield_files: List[str] = ['amber14-all.xml', 'amber14/tip3p.xml'], water_model: str = 'tip3p', pH: float = 7.0, hydrogen_variants: Optional[List[str]] = None, ionic_strength: float = 0.15, box_padding: float = 1.0, custom_bonds: List[str] = [], custom_angles: List[str] = [], custom_torsions: List[str] = [], minimize_only: bool = False, device: str = 'gpu', extra_args: str = "", dry_run: bool = False):
        """
        Initialize a RFDiffusionAA task.

        Args:
            input_file (str): The input structure file.
            output_prefix (str): The prefix for the output files. (Output files will be named [prefix].pdb and [prefix].dcd, etc).
            duration (int): The duration of the simulation in nanoseconds.
            relax_duration (int): The duration of the relaxation phase in nanoseconds. Default is 1 ns.
            output_frequency (int): The frequency of output frames in nanoseconds. Default is 1 ns.
            ligand_files (list): A list of ligand SDF files for assigning bond orders and protons to ligands in the input structure. Default is no ligands added.
            forcefield_files (list): A list of forcefield files to use. Default is ['amber14-all.xml', 'amber14/tip3p.xml'].
            water_model (str): The water model to use. Default is 'tip3p'.
            pH (float): The pH of the system for protonation. Default is 7.0.
            hydrogen_variants (str): List of hydrogen variants to use. Specify the chain and residue number, then the variant. ["A:13=HIE", "B:98=ASH", ...].
            ionic_strength (float): The ionic strength of the system in M. Default is 0.15 M.
            box_padding (float): The padding to add to the box size in nm. Default is 1.0 nm.
            custom_bonds (list): List of custom harmonic bonds to add. Specify the two atoms in the form           
                                    'chain_id:residue_id:atom_name', then the k-value (units kcal/(mol*A^2)), then    
                                    the target distance (units A).  ["A:99:CG,B:99:CG,10.0,2.0", ... ]
            custom_angles (list): List of custom harmonic angles to add. Specify the three atoms in the form       
                                    'chain_id:residue_id:atom_name', then the k-value (units kcal/mol), the           
                                    periodicity, and the phase value (units DEGREES).  E.g.                           
                                    ["A:99:CG,B:99:CG,C:99:CG,10.0,180", ... ]
            custom_torsions (list): List of custom harmonic angles to add. Specify the three atoms in the form       
                                    'chain_id:residue_id:atom_name', then the k-value (units kcal/mol), the          
                                    periodicity, and the phase value (units DEGREES).  E.g.                          
                                    ["A:99:CG,B:99:CG,C:99:CG,D:1:H6,10.0,1,180", ... ]
            minimize_only (bool): Whether to ONLY run a minimization, no simulation. Default is False.
            device (str): The device to run the task on. Default is 'gpu'.
            extra_args (str): Additional arguments for the task.

        Returns:
            None
            Outputs the following files using output_prefix:
                - [output_prefix].pdb: Topology and final coordinate structure file
                - [output_prefix].dcd: Trajectory coordinate file containing simulation frames
        """
        # Initialize the Task class
        super().__init__(device=device, extra_args=extra_args)

        # This Task name matches the name in the tasks.json file
        self.task_name = "EasyMD"
        
        # Your arguments here:
        self.input_file = input_file
        self.output_prefix = output_prefix
        self.duration = duration
        self.relax_duration = relax_duration
        self.output_frequency = output_frequency
        self.ligand_files = "".join([f" -l {ligand}" for ligand in ligand_files])
        self.forcefield_files = " ".join([f" -f {forcefield}" for forcefield in forcefield_files])
        self.water_model = water_model
        self.pH = pH
        self.hydrogen_variants = " ".join([f" -hv {variant}" for variant in hydrogen_variants]) if hydrogen_variants else ""
        self.ionic_strength = ionic_strength
        self.box_padding = box_padding
        self.custom_bonds = " ".join([f" -cb {bond}" for bond in custom_bonds])
        self.custom_angles = " ".join([f" -ca {angle}" for angle in custom_angles])
        self.custom_torsions = " ".join([f" -ct {torsion}" for torsion in custom_torsions])
        self.minimize_only = "--minimize-only" if minimize_only else ""
        self.device = device
        self.dry_run = dry_run

    def _run_dry(self):
        Path(self.output_prefix).parent.mkdir(parents=True, exist_ok=True)
        import shutil
        shutil.copy(self.input_file, f"{self.output_prefix}.pdb")
        with open(f"{self.output_prefix}.dcd", "w") as f:
            f.write("")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return
        
        # Make directories:
        Path(self.output_prefix).parent.mkdir(parents=True, exist_ok=True)
        
        # Run the task:
        self._run_task(self.task_name, 
                    input_file = self.input_file,
                    output_prefix = str(self.output_prefix), 
                    duration = self.duration,
                    relax_duration = self.relax_duration,
                    output_frequency = self.output_frequency,
                    ligand_files = self.ligand_files,
                    forcefield_files = self.forcefield_files,
                    water_model = self.water_model,
                    pH = self.pH,
                    hydrogen_variants = self.hydrogen_variants,
                    ionic_strength = self.ionic_strength,
                    box_padding = self.box_padding,
                    custom_bonds = self.custom_bonds,
                    custom_angles = self.custom_angles,
                    custom_torsions = self.custom_torsions,
                    minimize_only = self.minimize_only,
                    extra_args = self.extra_args,
                    device = self.device)
        
        return 
 
class RosettaLigandPrepare(Task):
    def __init__(self, input_file: Union[str, Path], output_dir: Union[str, Path] = '.', 
                 name: Optional[str] = None, prefix: Optional[str] = None, centroid: bool = False,
                 chain: Optional[str] = None, center: Optional[str] = None, max_confs: Optional[int] = None,
                 root_atom: Optional[int] = None, nbr_atom: Optional[int] = None, kinemage: Optional[str] = None,
                 amino_acid: Optional[str] = None, clobber: bool = False, no_param: bool = False, no_pdb: bool = False,
                 extra_torsion_output: bool = False, keep_names: bool = False, long_names: bool = False,
                 recharge: Optional[int] = None, m_ctrl: Optional[str] = None, mm_as_virt: bool = False,
                 skip_bad_conformers: bool = False, conformers_in_one_file: bool = False,
                 device: str = 'cpu', extra_args: str = "", dry_run: bool = False):
        """
        Initialize a RosettaLigandPrepare task.

        Args:
            input_file (str or Path): Must specify input .mol, .sdf, or .mol2 file!
            output_dir (str or Path): Directory where parameter and coordinates files will be saved. Default is '.'.
            name (str, optional): Name ligand residues NM1,NM2,... instead of LG1,LG2,...
            prefix (str, optional): Prefix for PDB file names.
            centroid (bool): Write files for Rosetta centroid mode too.
            chain (str, optional): The chain letter to use for the output PDB ligand.
            center (str, optional): Translate output PDB coords to have given heavy-atom centroid (format 'X,Y,Z').
            max_confs (int, optional): Don't expand proton chis if above this many total confs.
            root_atom (int, optional): Which atom in the molfile is the root? (indexed from 1).
            nbr_atom (int, optional): Which atom in the molfile is the nbr atom? (indexed from 1).
            kinemage (str, optional): Write ligand topology to FILE.
            amino_acid (str, optional): Set up params file for modified amino acid; .mol2 only; edit chis afterward. Implies --keep-names.
            clobber (bool): Overwrite existing files.
            no_param (bool): Skip writing .params files (for debugging).
            no_pdb (bool): Skip writing .pdb files (for debugging).
            extra_torsion_output (bool): Writing additional torsion files.
            keep_names (bool): Leaves atom names untouched except for duplications.
            long_names (bool): If specified name is longer than 3 letters, keep entire name in param NAME field.
            recharge (int, optional): Ignore existing partial charges, setting total charge to CHG.
            m_ctrl (str, optional): Read additional M control lines from FILE.
            mm_as_virt (bool): Assign mm atom types as VIRT, rather than X.
            skip_bad_conformers (bool): If a conformer has atoms in the wrong order, skip it and continue rather than dying.
            conformers_in_one_file (bool): Output 1st conformer to NAME.pdb and all others to NAME_conformers.pdb.
            device (str): The device to run the task on. Default is 'cpu'.
            extra_args (str): Additional arguments for the task.

        Returns:
            None
            Outputs the following files under output_dir:
                - [name].params: Rosetta parameter topology file for the ligand
                - [name]_0001.pdb: Conformer structure file in PDB format
        """
        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "RosettaLigandPrepare"
        
        self.input_file = input_file
        self.output_dir = output_dir
        self.name = name
        self.prefix = prefix
        self.centroid = centroid
        self.chain = chain
        self.center = center
        self.max_confs = max_confs
        self.root_atom = root_atom
        self.nbr_atom = nbr_atom
        self.kinemage = kinemage
        self.amino_acid = amino_acid
        self.clobber = clobber
        self.no_param = no_param
        self.no_pdb = no_pdb
        self.extra_torsion_output = extra_torsion_output
        self.keep_names = keep_names
        self.long_names = long_names
        self.recharge = recharge
        self.m_ctrl = m_ctrl
        self.mm_as_virt = mm_as_virt
        self.skip_bad_conformers = skip_bad_conformers
        self.conformers_in_one_file = conformers_in_one_file
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        name = self.name or "LG1"
        
        with open(self.output_dir / f"{name}.params", "w") as f:
            f.write(f"NAME {name}\nIO_STRING {name} L\nTYPE LIGAND\n")
            
        shutil.copy(RESOURCES_DIR / "dummy.pdb", self.output_dir / f"{name}_0001.pdb")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        import os
        # Make directories:
        self.output_dir = make_directory(self.output_dir)

        # Resolve paths to absolute paths
        input_file_abs = Path(self.input_file).resolve()
        if not input_file_abs.exists():
            raise FileNotFoundError(f"Input molecular file does not exist: {self.input_file}")

        m_ctrl_abs = Path(self.m_ctrl).resolve() if self.m_ctrl else None
        if m_ctrl_abs and not m_ctrl_abs.exists():
            raise FileNotFoundError(f"M control file does not exist: {self.m_ctrl}")

        # Construct ligand flags
        flags = []
        if self.name:
            flags.append(f"-n {self.name}")
        if self.prefix:
            flags.append(f"-p {self.prefix}")
        if self.centroid:
            flags.append("-c")
        if self.chain:
            flags.append(f"--chain={self.chain}")
        if self.center:
            flags.append(f"--center={self.center}")
        if self.max_confs is not None:
            flags.append(f"-m {self.max_confs}")
        if self.root_atom is not None:
            flags.append(f"--root_atom={self.root_atom}")
        if self.nbr_atom is not None:
            flags.append(f"--nbr_atom={self.nbr_atom}")
        if self.kinemage:
            flags.append(f"-k {self.kinemage}")
        if self.amino_acid:
            flags.append(f"-a {self.amino_acid}")
        if self.clobber:
            flags.append("--clobber")
        if self.no_param:
            flags.append("--no-param")
        if self.no_pdb:
            flags.append("--no-pdb")
        if self.extra_torsion_output:
            flags.append("--extra_torsion_output")
        if self.keep_names:
            flags.append("--keep-names")
        if self.long_names:
            flags.append("--long-names")
        if self.recharge is not None:
            flags.append(f"--recharge={self.recharge}")
        if m_ctrl_abs:
            flags.append(f"--m-ctrl={m_ctrl_abs}")
        if self.mm_as_virt:
            flags.append("--mm-as-virt")
        if self.skip_bad_conformers:
            flags.append("--skip-bad-conformers")
        if self.conformers_in_one_file:
            flags.append("--conformers-in-one-file")

        ligand_flags = " ".join(flags)

        # Run the task by changing working directories to output_dir
        orig_cwd = os.getcwd()
        os.chdir(self.output_dir)
        try:
            self._run_task(
                self.task_name,
                input_file=str(input_file_abs),
                ligand_flags=ligand_flags,
                extra_args=self.extra_args,
                device=self.device
            )
        finally:
            os.chdir(orig_cwd)

        return
 
class Custom(Task):
    def __init__(self, command: str, container: str = 'Ribbon', device: str = 'cpu', dry_run: bool = False):
        """
        Initialize a Custom task.
        This allows the user to run a custom command in a specified container.
        N.B. This allows the user to run arbitrary code; use with caution.

        Args:
            command (str): The command to run.
            container (str): The container to run the command in. Default is 'Ribbon'.
            device (str): The device to run the task on. Default is 'cpu'.
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Custom"
        
        # Task-specific variables
        self.command = command
        self.container = container
        self.device = device
        self.dry_run = dry_run

    def _run_dry(self):
        print(f"[DRY RUN] Custom command: {self.command}")
        import re
        match = re.search(r">\s*(.+)$", self.command)
        if match:
            out_file = Path(match.group(1).strip("'\" "))
            out_file.parent.mkdir(parents=True, exist_ok=True)
            with open(out_file, "w") as f:
                f.write("done\n")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        # Run the task
        self._run_task(
            self.task_name,
            command=self.command,
            container_override=self.container,
            device=self.device
        )