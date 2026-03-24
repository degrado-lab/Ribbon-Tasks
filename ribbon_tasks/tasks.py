from ribbon.utils import make_directories, make_directory, list_files
from pathlib import Path
from ribbon.runner import Task
from typing import List, Union, Optional, Any, Tuple
import tempfile
import shutil
import json

class LigandMPNN(Task):
    def __init__(self, output_dir: Union[str, Path], structure_list : List[Union[str, Path]], num_designs: int = 1, temperature: float =0.1, device: str = 'cpu', extra_args: str = ""):
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
            Outputs the following directories:
                output_dir:
                    - backbones: The generated backbones
                    - packed: The packed structures, including sidechains
                    - sequences: The generated sequences as FASTA files. Multiple chains are separated by ':'. Each line is a different design.
                    - seqs_split: The sequences split into separate FASTA files, one per design. Each line is a different chain.
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

    def run(self):
        
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
    def __init__(self, output_dir: Union[str, Path], structure_list: List[Union[str, Path]], num_designs: int = 1, temperature: float = 0.000001, device: str = 'cpu', fix_beta: bool = False, extra_args: str = ""):
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
            Outputs the following directories:
                output_dir:
                    - backbones: The generated backbones
                    - packed: The packed structures, including sidechains
                    - sequences: The generated sequences as FASTA files. Multiple chains are separated by ':'. Each line is a different design.
                    - seqs_split: The sequences split into separate FASTA files, one per design. Each line is a different chain.
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

    def run(self):
        
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
    def __init__(self, output_dir: Union[str, Path], pdb_input_file: Optional[Union[str, Path]] = None, pdb_input_dir: Optional[Union[str, Path]] = None, nstructs: int = 1, device: str = 'cpu', extra_args: str = ""):
        """
        Initialize a FastRelax task.

        Args:
            output_dir (str): The directory to save the output files.
            pdb_input_file (str, optional): Path to a single PDB file. Default is None.
            pdb_input_dir (str, optional): Path to a directory containing PDB files. Default is None.
            device (str): The device to run the task on. Default is 'cpu'.

        Raises:
            ValueError: If neither pdb_input_file nor pdb_input_dir is specified.
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
        self.device = device
        self.extra_args = extra_args

    def run(self):
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

        # Run the task
        self._run_task(
            self.task_name,
            pdb_string=pdb_string,
            output_dir=str(self.output_dir),
            nstructs=self.nstructs,
            extra_args=self.extra_args,
            device=self.device
        )

class Chai1(Task):
    def __init__(self, fasta_file: Union[str, Path], output_dir: Union[str, Path] = '.', smiles_string: Optional[str] = None, num_ligands: int = 1, device: str = 'gpu'):
        """
        Initialize a Chai-1 task.

        Args:
            fasta_file (str): The FASTA file containing the protein sequence (no ligand).
            output_dir (str): The directory to save the output files. Default is '.'.
            smiles_string (str, optional): The SMILES string of the ligand. Default is None.
            num_ligands (int): The number of ligands. Default is 1.
            device (str): The device to run the task on. Default is 'gpu'.
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

    def run(self):
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
    def __init__(self, fasta_file: Union[str, Path], output_dir: Union[str, Path] = '.', smiles_list: List[str] = [], calculate_binding: bool = False, use_msa_server: bool = True, device: str = 'gpu', extra_args: str = ""):
        """
        Initialize a Boltz-2 task.

        Args:
            fasta_file (str): The FASTA file containing the protein sequence (no ligand).
            output_dir (str): The directory to save the output files. Default is '.'.
            smiles_list (list of str, optional): List of the SMILES strings of the ligands. Repeat a SMILES string to use multiple copies of the same ligand. 
                Default is no ligands.
            calculate_binding (bool, optional): Whether to calculate binding properties of the ligand. Only calculates for the first ligand in the list. Default is False.
            device (str): The device to run the task on. Default is 'gpu'.
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

    def run(self):
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
    def __init__(self, fasta_file_or_dir: Union[str, Path], output_dir: Union[str, Path] = '.', param: str = 'RaptorX-Single-ESM1b.pt', device: str = 'gpu', extra_args: str = ""):
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

    def run(self):
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
                 atom2: str, output_file: Union[str, Path], device: str = 'cpu'):
        """
        Initialize a CalculateDistance task.
        This calculates the distance between two atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            atom1_name (str): Name of the first atom in the format 'Chain:Residue:Atom'.
            atom2_name (str): Name of the second atom in the format 'Chain:Residue:Atom'.
            output_file (str): Path to the output file. Suffixed with '.dist'.
            device (str): The device to run the task on. Default is 'cpu'.
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

    def run(self):
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
                 atom2: str, atom3: str, output_file: Union[str, Path], device: str = 'cpu'):
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

    def run(self):
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
                 atom2: str, atom3: str, atom4: str, output_file: Union[str, Path], device: str = 'cpu'):
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

    def run(self):
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
                 atom_list_B: List[str], output_file: Union[str, Path], average: bool = False, device: str = 'cpu'):
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

    def run(self):
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
    def __init__(self, input_file: Union[str, Path], output_file: Union[str, Path], selection: str = 'all'):
        """
        Initialize a CalculateDistance task.
        This calculates the distance between two atoms in a PDB file.

        Args:
            input_file (str): Path to the PDB or CIF file.
            output_file (str): Path to the output file.
            selection (str): PyMol selection string to specify what to add hydrogens to. Default is 'all'.
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

    def run(self):
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
    def __init__(self, pdb_input_file: Union[str, Path], pdb_output_file: Union[str, Path], flip: bool = False, custom_ligands: List[Tuple[str, Union[str, Path]]] = []):
        """
        Add Hydrogens to a PDB file.

        Args:
            pdb_input_file (str): Path to the input PDB file.
            pdb_output_file (str): Path to the output PDB file.
            flip (bool): Whether to optionally flip N/Q/H residues. Default False.
            custom_ligands (list of tuples): List of custom ligand resnames and SDF files. (E.g. [('KP1', 'kemp1.sdf'), ...] ) Only necessary if there is a ligand which is not already in the Protein Data Bank.
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

    def run(self):
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
    def __init__(self, pdb_file: Union[str, Path], output_file: Union[str, Path], atom_1: str, device: str = 'cpu'):
        """
        Initialize a CalculateSASA task.
        This calculates the Solvent Accessible Surface Area for a set of atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            output_file (str): Path to the output file. Suffixed with '.angle'.
            atom_1 (str): Atom specification in format chain_id:res_id:atom_name.
            device (str): The device to run the task on. Default is 'cpu'.

        TODO:
            - Implement the task script in ribbon/ribbon_tasks/task_scripts/calculate_sasa.py
        """

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
    def __init__(self, input_structure: Union[str, Path], output_dir: Union[str, Path], contig_map: str, num_designs: int = 1, total_length: str = 'null', ligand: str = 'null',  diffuser_steps: int = 200, deterministic: bool = False, design_startnum: int = 0, force: bool = False, device: str = 'gpu', extra_args: str = ""):
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
            Outputs designs in the output_dir.
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
        # we should have a final length flag for config.length

    def run(self):
        
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
        
    def __init__(self, input_file: Union[str, Path], output_prefix: Union[str, Path], duration: int, relax_duration: int = 1, output_frequency: int = 1, ligand_files: List[Union[str, Path]] = [], forcefield_files: List[str] = ['amber14-all.xml', 'amber14/tip3p.xml'], water_model: str = 'tip3p', pH: float = 7.0, hydrogen_variants: Optional[List[str]] = None, ionic_strength: float = 0.15, box_padding: float = 1.0, custom_bonds: List[str] = [], custom_angles: List[str] = [], custom_torsions: List[str] = [], minimize_only: bool = False, device: str = 'gpu', extra_args: str = ""):
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
            Outputs PDB, DCD files using the given output_prefix.
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

    def run(self):
        
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
       
class Custom(Task):
    def __init__(self, command: str, container: str = 'Ribbon', device: str = 'cpu'):
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

    def run(self):

        # Run the task
        self._run_task(
            self.task_name,
            command=self.command,
            container_override=self.container,
            device=self.device
        )