from ribbon.utils import make_directories, make_directory, list_files
from pathlib import Path
from ribbon.runner import Task
import tempfile
import shutil
import json

class LigandMPNN(Task):
    def __init__(self, output_dir, structure_list, num_designs=1, device='cpu', extra_args=""):
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
                    extra_args = self.extra_args,
                    device = self.device)
        
        # Split the FASTA files:
        for file in (self.output_dir / 'seqs').iterdir():
            print(f'Splitting {file}')
            split_ligandmpnn_fasta(file, self.output_dir / 'seqs_split')
        
        return 

class FastRelax(Task):
    def __init__(self, output_dir, pdb_input_file=None, pdb_input_dir=None, device='cpu'):
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
        self.device = device

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
            device=self.device
        )

class Chai1(Task):
    def __init__(self, fasta_file, output_dir='.', smiles_string=None, num_ligands=1, device='gpu'):
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

class RaptorXSingle(Task):
    def __init__(self, fasta_file_or_dir, output_dir='.', param='RaptorX-Single-ESM1b.pt', device='gpu', extra_args=""):
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
    def __init__(self, pdb_file, chain1_id, res1_id, atom1_name,
                 chain2_id, res2_id, atom2_name, output_file, device='cpu'):
        """
        Initialize a CalculateDistance task.
        This calculates the distance between two atoms in a PDB file.

        Args:
            pdb_file (str): Path to the PDB file.
            chain1_id (str): Chain ID of the first atom.
            res1_id (str): Residue ID of the first atom.
            atom1_name (str): Name of the first atom.
            chain2_id (str): Chain ID of the second atom.
            res2_id (str): Residue ID of the second atom.
            atom2_name (str): Name of the second atom.
            output_file (str): Path to the output file. Suffixed with '.dist'.
            device (str): The device to run the task on. Default is 'cpu'.
        """
        # Initialize the Task class
        super().__init__()

        # This Task name matches the name in the tasks.json file
        self.task_name = "Calculate Distance"
        
        # Task-specific variables
        self.pdb_file = pdb_file
        self.chain1_id = chain1_id
        self.res1_id = res1_id
        self.atom1_name = atom1_name
        self.chain2_id = chain2_id
        self.res2_id = res2_id
        self.atom2_name = atom2_name
        self.output_file = output_file
        self.device = device

    def run(self):
        # Ensure output directory exists
        make_directory(Path(self.output_file).parent)

        # Run the task
        self._run_task(
            self.task_name,
            pdb_file=self.pdb_file,
            chain1_id=self.chain1_id,
            res1_id=self.res1_id,
            atom1_name=self.atom1_name,
            chain2_id=self.chain2_id,
            res2_id=self.res2_id,
            atom2_name=self.atom2_name,
            output_file=self.output_file,
            device=self.device
        )

class CalculateSASA(Task):
    def __init__(self, pdb_file, output_file, atom_1, device='cpu'):
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


class Custom(Task):
    def __init__(self, command, container='Ribbon', device='cpu'):
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