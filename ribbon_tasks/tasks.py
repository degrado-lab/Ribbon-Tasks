from ribbon.utils import make_directories, make_directory, list_files
from pathlib import Path
from ribbon.runner import Task
from typing import List, Union, Optional, Any, Tuple, Dict
import tempfile
import shutil
import json
import os

RESOURCES_DIR = Path(__file__).parent / "resources"

class LigandMPNN(Task):
    def __init__(
        self,
        output_dir: Union[str, Path],
        structure_list: List[Union[str, Path]],
        num_designs: int = 1,
        temperature: float = 0.1,
        device: str = "cpu",
        model_type: Optional[str] = None,
        checkpoint_protein_mpnn: Optional[Union[str, Path]] = None,
        checkpoint_ligand_mpnn: Optional[Union[str, Path]] = None,
        checkpoint_per_residue_label_membrane_mpnn: Optional[Union[str, Path]] = None,
        checkpoint_global_label_membrane_mpnn: Optional[Union[str, Path]] = None,
        checkpoint_soluble_mpnn: Optional[Union[str, Path]] = None,
        checkpoint_path_sc: Optional[Union[str, Path]] = None,
        fasta_seq_separation: Optional[str] = None,
        verbose: Optional[Union[int, bool]] = None,
        fixed_residues: Optional[Union[str, List[str]]] = None,
        fixed_residues_multi: Optional[Union[str, Path, dict]] = None,
        redesigned_residues: Optional[Union[str, List[str]]] = None,
        redesigned_residues_multi: Optional[Union[str, Path, dict]] = None,
        bias_AA: Optional[Union[str, dict]] = None,
        bias_AA_per_residue: Optional[Union[str, Path, dict]] = None,
        bias_AA_per_residue_multi: Optional[Union[str, Path, dict]] = None,
        omit_AA: Optional[Union[str, List[str]]] = None,
        omit_AA_per_residue: Optional[Union[str, Path, dict]] = None,
        omit_AA_per_residue_multi: Optional[Union[str, Path, dict]] = None,
        symmetry_residues: Optional[Union[str, List[List[str]]]] = None,
        symmetry_weights: Optional[Union[str, List[List[Union[float, int]]]]] = None,
        homo_oligomer: Optional[Union[int, bool]] = None,
        file_ending: Optional[str] = None,
        zero_indexed: Optional[Union[int, bool]] = None,
        seed: Optional[int] = None,
        batch_size: Optional[int] = None,
        number_of_batches: Optional[int] = None,
        save_stats: Optional[Union[int, bool]] = None,
        ligand_mpnn_use_atom_context: Optional[Union[int, bool]] = None,
        ligand_mpnn_cutoff_for_score: Optional[float] = None,
        ligand_mpnn_use_side_chain_context: Optional[Union[int, bool]] = None,
        chains_to_design: Optional[Union[str, List[str]]] = None,
        parse_these_chains_only: Optional[Union[str, List[str]]] = None,
        transmembrane_buried: Optional[Union[str, List[str]]] = None,
        transmembrane_interface: Optional[Union[str, List[str]]] = None,
        global_transmembrane_label: Optional[Union[int, bool]] = None,
        parse_atoms_with_zero_occupancy: Optional[Union[int, bool]] = None,
        pack_side_chains: Optional[Union[int, bool]] = None,
        number_of_packs_per_design: Optional[int] = None,
        sc_num_denoising_steps: Optional[int] = None,
        sc_num_samples: Optional[int] = None,
        repack_everything: Optional[Union[int, bool]] = None,
        force_hetatm: Optional[Union[int, bool]] = None,
        packed_suffix: Optional[str] = None,
        pack_with_ligand_context: Optional[Union[int, bool]] = None,
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        """
        Initialize a LigandMPNN task.

        Args:
            output_dir (str or Path): The directory to save the output files.
            structure_list (list): A list of pdb or cif files to use as input structures.
            num_designs (int): The number of designs to generate per pass (batch size). Default is 1.
            temperature (float): Sampling temperature for sequence generation. Default is 0.1.
            device (str): The device to run the task on ('cpu', 'gpu', 'gpu_wsl'). Default is 'cpu'.
            model_type (str, optional): Choose model: 'protein_mpnn', 'ligand_mpnn',
                'per_residue_label_membrane_mpnn', 'global_label_membrane_mpnn', 'soluble_mpnn'.
            checkpoint_protein_mpnn (str or Path, optional): Path to protein_mpnn model weights.
            checkpoint_ligand_mpnn (str or Path, optional): Path to ligand_mpnn model weights.
            checkpoint_per_residue_label_membrane_mpnn (str or Path, optional): Path to membrane model weights.
            checkpoint_global_label_membrane_mpnn (str or Path, optional): Path to global membrane model weights.
            checkpoint_soluble_mpnn (str or Path, optional): Path to soluble_mpnn model weights.
            checkpoint_path_sc (str or Path, optional): Path to side chain packer model weights.
            fasta_seq_separation (str, optional): Symbol between sequences from different chains. Default is ':'.
            verbose (int or bool, optional): Print verbose output (0 or 1). Default is 1.
            fixed_residues (str or list of str, optional): Residues to fix, e.g. "A12 A13 B2" or ["A12", "A13"].
            fixed_residues_multi (str, Path, or dict, optional): JSON path or dict mapping PDB paths to fixed residues.
            redesigned_residues (str or list of str, optional): Residues to redesign, fixing all others.
            redesigned_residues_multi (str, Path, or dict, optional): JSON path or dict mapping PDB paths to redesigned residues.
            bias_AA (str or dict, optional): Amino acid bias, e.g. "A:-1.024,P:2.34" or {"A": -1.024, "P": 2.34}.
            bias_AA_per_residue (str, Path, or dict, optional): Per-residue bias mapping, e.g. {"A12": {"G": -0.3}}.
            bias_AA_per_residue_multi (str, Path, or dict, optional): Bias mapping per PDB.
            omit_AA (str or list of str, optional): Global amino acids to omit, e.g. "ACG" or ["A", "C", "G"].
            omit_AA_per_residue (str, Path, or dict, optional): Per-residue omitted amino acids.
            omit_AA_per_residue_multi (str, Path, or dict, optional): Per-residue omitted amino acids per PDB.
            symmetry_residues (str or list of lists, optional): Residues to enforce symmetry, e.g. "A12,A13|C2,C3".
            symmetry_weights (str or list of lists, optional): Symmetry weights matching symmetry_residues.
            homo_oligomer (int or bool, optional): Enable homo-oligomer design with equal weighting.
            file_ending (str, optional): String to append to output filenames.
            zero_indexed (int or bool, optional): Start output PDB numbering with 0 instead of 1.
            seed (int, optional): Random seed for reproducibility.
            batch_size (int, optional): Alias for num_designs.
            number_of_batches (int, optional): Number of batches to run. Default is 1.
            save_stats (int or bool, optional): Save output score statistics into stats/ folder.
            ligand_mpnn_use_atom_context (int or bool, optional): 1 - use atom context, 0 - do not. Default is 1.
            ligand_mpnn_cutoff_for_score (float, optional): Cutoff in Angstroms for reporting score. Default is 8.0.
            ligand_mpnn_use_side_chain_context (int or bool, optional): Use side chain atoms as ligand context. Default is 0.
            chains_to_design (str or list of str, optional): Chains to redesign, keeping others fixed, e.g. "A,B" or ["A", "B"].
            parse_these_chains_only (str or list of str, optional): Chains to parse from input structures.
            transmembrane_buried (str or list of str, optional): Buried residues for membrane MPNN model.
            transmembrane_interface (str or list of str, optional): Interface residues for membrane MPNN model.
            global_transmembrane_label (int or bool, optional): 1 - transmembrane, 0 - soluble.
            parse_atoms_with_zero_occupancy (int or bool, optional): Parse zero-occupancy atoms from input PDBs.
            pack_side_chains (int or bool, optional): Run side chain packer (1 or 0).
            number_of_packs_per_design (int, optional): Side chain packing samples per design. Default is 4.
            sc_num_denoising_steps (int, optional): Denoising/recycling steps for side chain packing. Default is 3.
            sc_num_samples (int, optional): Mixture distribution samples for packer. Default is 16.
            repack_everything (int or bool, optional): Repack all side chains including fixed residues.
            force_hetatm (int or bool, optional): Force ligand atoms to be written as HETATM in packed PDBs.
            packed_suffix (str, optional): Suffix for packed PDB files. Default is '_packed'.
            pack_with_ligand_context (int or bool, optional): Pack side chains using ligand context. Default is 1.
            extra_args (str): Additional arguments passed directly to the command line. Default is "".
            dry_run (bool): If True, generate mock outputs without running container. Default is False.

        Returns:
            None
            Outputs the following directories under output_dir:
                - backbones/: The generated backbones
                - packed/: The packed structures, including sidechains (when packing enabled)
                - stats/: Output statistics .pt files (when save_stats enabled)
                - sequences/ (or seqs/): The generated sequences as FASTA files
                - seqs_split/: The sequences split into separate FASTA files, one per design
        """
        # Initialize the Task class
        super().__init__(device=device, extra_args=extra_args)

        # This Task name matches the name in the tasks.json file
        self.task_name = "LigandMPNN"

        # Handle aliases and case variants
        if batch_size is not None:
            num_designs = batch_size
        if bias_AA is None:
            bias_AA = kwargs.pop("bias_aa", None)
        if bias_AA_per_residue is None:
            bias_AA_per_residue = kwargs.pop("bias_aa_per_residue", None)
        if bias_AA_per_residue_multi is None:
            bias_AA_per_residue_multi = kwargs.pop("bias_aa_per_residue_multi", None)
        if omit_AA is None:
            omit_AA = kwargs.pop("omit_aa", None)
        if omit_AA_per_residue is None:
            omit_AA_per_residue = kwargs.pop("omit_aa_per_residue", None)
        if omit_AA_per_residue_multi is None:
            omit_AA_per_residue_multi = kwargs.pop("omit_aa_per_residue_multi", None)

        if kwargs:
            raise TypeError(f"LigandMPNN.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}")

        self.output_dir = output_dir
        self.structure_list = structure_list
        self.num_designs = num_designs
        self.temperature = temperature
        self.model_type = model_type
        self.checkpoint_protein_mpnn = checkpoint_protein_mpnn
        self.checkpoint_ligand_mpnn = checkpoint_ligand_mpnn
        self.checkpoint_per_residue_label_membrane_mpnn = checkpoint_per_residue_label_membrane_mpnn
        self.checkpoint_global_label_membrane_mpnn = checkpoint_global_label_membrane_mpnn
        self.checkpoint_soluble_mpnn = checkpoint_soluble_mpnn
        self.checkpoint_path_sc = checkpoint_path_sc
        self.fasta_seq_separation = fasta_seq_separation
        self.verbose = verbose
        self.fixed_residues = fixed_residues
        self.fixed_residues_multi = fixed_residues_multi
        self.redesigned_residues = redesigned_residues
        self.redesigned_residues_multi = redesigned_residues_multi
        self.bias_AA = bias_AA
        self.bias_AA_per_residue = bias_AA_per_residue
        self.bias_AA_per_residue_multi = bias_AA_per_residue_multi
        self.omit_AA = omit_AA
        self.omit_AA_per_residue = omit_AA_per_residue
        self.omit_AA_per_residue_multi = omit_AA_per_residue_multi
        self.symmetry_residues = symmetry_residues
        self.symmetry_weights = symmetry_weights
        self.homo_oligomer = homo_oligomer
        self.file_ending = file_ending
        self.zero_indexed = zero_indexed
        self.seed = seed
        self.number_of_batches = number_of_batches
        self.save_stats = save_stats
        self.ligand_mpnn_use_atom_context = ligand_mpnn_use_atom_context
        self.ligand_mpnn_cutoff_for_score = ligand_mpnn_cutoff_for_score
        self.ligand_mpnn_use_side_chain_context = ligand_mpnn_use_side_chain_context
        self.chains_to_design = chains_to_design
        self.parse_these_chains_only = parse_these_chains_only
        self.transmembrane_buried = transmembrane_buried
        self.transmembrane_interface = transmembrane_interface
        self.global_transmembrane_label = global_transmembrane_label
        self.parse_atoms_with_zero_occupancy = parse_atoms_with_zero_occupancy
        self.pack_side_chains = pack_side_chains
        self.number_of_packs_per_design = number_of_packs_per_design
        self.sc_num_denoising_steps = sc_num_denoising_steps
        self.sc_num_samples = sc_num_samples
        self.repack_everything = repack_everything
        self.force_hetatm = force_hetatm
        self.packed_suffix = packed_suffix
        self.pack_with_ligand_context = pack_with_ligand_context
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        seqs_dir = self.output_dir / 'seqs'
        seqs_dir.mkdir(parents=True, exist_ok=True)

        file_ending = self.file_ending or ""
        total_batches = self.number_of_batches if self.number_of_batches is not None else 1
        total_seqs = self.num_designs * total_batches

        for pdb_path in self.structure_list:
            pdb_stem = Path(pdb_path).stem
            fasta_path = seqs_dir / f"{pdb_stem}{file_ending}.fa"
            with open(fasta_path, 'w') as f:
                f.write(">original\nAAAAA\n")
                for i in range(total_seqs):
                    f.write(f">design_{i}, score=0.1, temp={self.temperature}\nAAAAA\n")

        split_dir = self.output_dir / 'seqs_split'
        split_dir.mkdir(parents=True, exist_ok=True)
        sep = self.fasta_seq_separation or ":"
        for fasta_file in seqs_dir.iterdir():
            if not (fasta_file.is_file() and fasta_file.suffix in ['.fa', '.fasta']):
                continue
            with open(fasta_file) as f:
                lines = f.readlines()[2:]
                for i in range(0, len(lines), 2):
                    if not lines[i].startswith('>'):
                        continue
                    if i + 1 >= len(lines):
                        break
                    name = lines[i].strip()
                    chains = lines[i + 1].strip().split(sep)
                    index = 0
                    output_path = split_dir / f'{fasta_file.stem}_{index}.fasta'
                    while output_path.exists():
                        index += 1
                        output_path = split_dir / f'{fasta_file.stem}_{index}.fasta'
                    with open(output_path, 'w') as g:
                        for chain_index, chain in enumerate(chains):
                            name_with_chain = f"{name.split(',')[0]}_{chain_index}" + ', ' + ','.join(name.split(',')[1:])
                            g.write(f'{name_with_chain}\n{chain}\n')

        backbones_dir = self.output_dir / 'backbones'
        backbones_dir.mkdir(parents=True, exist_ok=True)
        is_zero_idx = (self.zero_indexed is True) or (self.zero_indexed is not None and str(self.zero_indexed) == "1")
        for pdb_path in self.structure_list:
            pdb_stem = Path(pdb_path).stem
            for i in range(total_seqs):
                idx = i if is_zero_idx else (i + 1)
                (backbones_dir / f"{pdb_stem}_{idx}{file_ending}.pdb").touch()

        is_pack = (self.pack_side_chains is True) or (self.pack_side_chains is not None and int(self.pack_side_chains) == 1)
        if is_pack:
            packed_dir = self.output_dir / 'packed'
            packed_dir.mkdir(parents=True, exist_ok=True)
            suffix = self.packed_suffix if self.packed_suffix is not None else "_packed"
            num_packs = self.number_of_packs_per_design if self.number_of_packs_per_design is not None else 4
            for pdb_path in self.structure_list:
                pdb_stem = Path(pdb_path).stem
                for d in range(total_seqs):
                    idx = d if is_zero_idx else (d + 1)
                    if num_packs == 0:
                        (packed_dir / f"{pdb_stem}{suffix}_{idx}{file_ending}.pdb").touch()
                    else:
                        for p in range(1, num_packs + 1):
                            (packed_dir / f"{pdb_stem}{suffix}_{idx}_{p}{file_ending}.pdb").touch()

        is_stats = (self.save_stats is True) or (self.save_stats is not None and int(self.save_stats) == 1)
        if is_stats:
            stats_dir = self.output_dir / 'stats'
            stats_dir.mkdir(parents=True, exist_ok=True)
            for pdb_path in self.structure_list:
                pdb_stem = Path(pdb_path).stem
                (stats_dir / f"{pdb_stem}.pt").touch()

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        ###### HELPER FUNCTIONS #######
        def split_ligandmpnn_fasta(fasta_file, split_output_dir):
            split_output_dir.mkdir(parents=True, exist_ok=True)
            sep = self.fasta_seq_separation or ":"
            with open(fasta_file) as f:
                lines = f.readlines()[2:]
                for i in range(0, len(lines), 2):
                    if not lines[i].startswith('>'):
                        continue
                    if i + 1 >= len(lines):
                        break
                    name = lines[i].strip()
                    chains = lines[i + 1].strip().split(sep)
                    index = 0
                    output_path = split_output_dir / f'{fasta_file.stem}_{index}.fasta'
                    while output_path.exists():
                        index += 1
                        output_path = split_output_dir / f'{fasta_file.stem}_{index}.fasta'

                    print(f'Writing to {output_path}')
                    with open(output_path, 'w') as g:
                        for chain_index, chain in enumerate(chains):
                            name_with_chain = f"{name.split(',')[0]}_{chain_index}" + ', ' + ','.join(name.split(',')[1:])
                            g.write(f'{name_with_chain}\n{chain}\n')

        def _bool_or_int_flag(flag_name, val):
            if val is None:
                return ""
            if isinstance(val, bool):
                int_val = 1 if val else 0
            else:
                int_val = int(val)
            return f"--{flag_name} {int_val}"

        def _str_flag(flag_name, val):
            if val is None or val == "":
                return ""
            return f'--{flag_name} "{val}"'

        def _val_flag(flag_name, val):
            if val is None or val == "":
                return ""
            return f"--{flag_name} {val}"

        temp_files = []
        try:
            def _handle_json_mapping(flag_name, val, is_residue_list_map=False):
                if val is None:
                    return ""
                if isinstance(val, dict):
                    temp_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json')
                    temp_files.append(temp_file.name)
                    if is_residue_list_map:
                        out_dict = {
                            str(k): (" ".join(str(x) for x in v) if isinstance(v, (list, tuple)) else str(v))
                            for k, v in val.items()
                        }
                    else:
                        out_dict = val
                    json.dump(out_dict, temp_file)
                    temp_file.close()
                    return f'--{flag_name} "{temp_file.name}"'
                else:
                    return f'--{flag_name} "{Path(val).resolve()}"'

            # Make directories:
            self.output_dir = make_directory(self.output_dir)

            pdb_input_json = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.json')
            temp_files.append(pdb_input_json.name)
            json.dump([str(Path(p).resolve()) for p in self.structure_list], pdb_input_json)
            pdb_input_json.close()

            # Format residue lists
            fixed_res_val = " ".join(str(x) for x in self.fixed_residues) if isinstance(self.fixed_residues, (list, tuple)) else self.fixed_residues
            redesigned_res_val = " ".join(str(x) for x in self.redesigned_residues) if isinstance(self.redesigned_residues, (list, tuple)) else self.redesigned_residues
            chains_to_design_val = ",".join(str(x) for x in self.chains_to_design) if isinstance(self.chains_to_design, (list, tuple)) else self.chains_to_design
            parse_chains_val = ",".join(str(x) for x in self.parse_these_chains_only) if isinstance(self.parse_these_chains_only, (list, tuple)) else self.parse_these_chains_only
            tm_buried_val = " ".join(str(x) for x in self.transmembrane_buried) if isinstance(self.transmembrane_buried, (list, tuple)) else self.transmembrane_buried
            tm_interface_val = " ".join(str(x) for x in self.transmembrane_interface) if isinstance(self.transmembrane_interface, (list, tuple)) else self.transmembrane_interface
            omit_aa_val = "".join(str(x) for x in self.omit_AA) if isinstance(self.omit_AA, (list, tuple)) else self.omit_AA

            # Format bias_AA
            if isinstance(self.bias_AA, dict):
                bias_aa_val = ",".join(f"{k}:{v}" for k, v in self.bias_AA.items())
            else:
                bias_aa_val = self.bias_AA

            # Format symmetry
            if isinstance(self.symmetry_residues, (list, tuple)):
                sym_res_parts = [",".join(str(x) for x in g) if isinstance(g, (list, tuple)) else str(g) for g in self.symmetry_residues]
                sym_res_val = "|".join(sym_res_parts)
            else:
                sym_res_val = self.symmetry_residues

            if isinstance(self.symmetry_weights, (list, tuple)):
                sym_w_parts = [",".join(str(x) for x in g) if isinstance(g, (list, tuple)) else str(g) for g in self.symmetry_weights]
                sym_w_val = "|".join(sym_w_parts)
            else:
                sym_w_val = self.symmetry_weights

            # Build flags
            model_type_flag = _val_flag("model_type", self.model_type)
            checkpoint_protein_mpnn_flag = f'--checkpoint_protein_mpnn "{Path(self.checkpoint_protein_mpnn).resolve()}"' if self.checkpoint_protein_mpnn is not None else ""
            checkpoint_ligand_mpnn_flag = f'--checkpoint_ligand_mpnn "{Path(self.checkpoint_ligand_mpnn).resolve()}"' if self.checkpoint_ligand_mpnn is not None else ""
            checkpoint_per_residue_label_membrane_mpnn_flag = f'--checkpoint_per_residue_label_membrane_mpnn "{Path(self.checkpoint_per_residue_label_membrane_mpnn).resolve()}"' if self.checkpoint_per_residue_label_membrane_mpnn is not None else ""
            checkpoint_global_label_membrane_mpnn_flag = f'--checkpoint_global_label_membrane_mpnn "{Path(self.checkpoint_global_label_membrane_mpnn).resolve()}"' if self.checkpoint_global_label_membrane_mpnn is not None else ""
            checkpoint_soluble_mpnn_flag = f'--checkpoint_soluble_mpnn "{Path(self.checkpoint_soluble_mpnn).resolve()}"' if self.checkpoint_soluble_mpnn is not None else ""
            checkpoint_path_sc_flag = f'--checkpoint_path_sc "{Path(self.checkpoint_path_sc).resolve()}"' if self.checkpoint_path_sc is not None else ""
            fasta_seq_separation_flag = _str_flag("fasta_seq_separation", self.fasta_seq_separation)
            verbose_flag = _bool_or_int_flag("verbose", self.verbose)
            fixed_residues_flag = _str_flag("fixed_residues", fixed_res_val)
            fixed_residues_multi_flag = _handle_json_mapping("fixed_residues_multi", self.fixed_residues_multi, is_residue_list_map=True)
            redesigned_residues_flag = _str_flag("redesigned_residues", redesigned_res_val)
            redesigned_residues_multi_flag = _handle_json_mapping("redesigned_residues_multi", self.redesigned_residues_multi, is_residue_list_map=True)
            bias_aa_flag = _str_flag("bias_AA", bias_aa_val)
            bias_aa_per_residue_flag = _handle_json_mapping("bias_AA_per_residue", self.bias_AA_per_residue)
            bias_aa_per_residue_multi_flag = _handle_json_mapping("bias_AA_per_residue_multi", self.bias_AA_per_residue_multi)
            omit_aa_flag = _str_flag("omit_AA", omit_aa_val)
            omit_aa_per_residue_flag = _handle_json_mapping("omit_AA_per_residue", self.omit_AA_per_residue)
            omit_aa_per_residue_multi_flag = _handle_json_mapping("omit_AA_per_residue_multi", self.omit_AA_per_residue_multi)
            symmetry_residues_flag = _str_flag("symmetry_residues", sym_res_val)
            symmetry_weights_flag = _str_flag("symmetry_weights", sym_w_val)
            homo_oligomer_flag = _bool_or_int_flag("homo_oligomer", self.homo_oligomer)
            file_ending_flag = _str_flag("file_ending", self.file_ending)
            zero_indexed_flag = _bool_or_int_flag("zero_indexed", self.zero_indexed)
            seed_flag = _val_flag("seed", self.seed)
            number_of_batches_flag = _val_flag("number_of_batches", self.number_of_batches)
            save_stats_flag = _bool_or_int_flag("save_stats", self.save_stats)
            ligand_mpnn_use_atom_context_flag = _bool_or_int_flag("ligand_mpnn_use_atom_context", self.ligand_mpnn_use_atom_context)
            ligand_mpnn_cutoff_for_score_flag = _val_flag("ligand_mpnn_cutoff_for_score", self.ligand_mpnn_cutoff_for_score)
            ligand_mpnn_use_side_chain_context_flag = _bool_or_int_flag("ligand_mpnn_use_side_chain_context", self.ligand_mpnn_use_side_chain_context)
            chains_to_design_flag = _str_flag("chains_to_design", chains_to_design_val)
            parse_these_chains_only_flag = _str_flag("parse_these_chains_only", parse_chains_val)
            transmembrane_buried_flag = _str_flag("transmembrane_buried", tm_buried_val)
            transmembrane_interface_flag = _str_flag("transmembrane_interface", tm_interface_val)
            global_transmembrane_label_flag = _bool_or_int_flag("global_transmembrane_label", self.global_transmembrane_label)
            parse_atoms_with_zero_occupancy_flag = _bool_or_int_flag("parse_atoms_with_zero_occupancy", self.parse_atoms_with_zero_occupancy)
            pack_side_chains_flag = _bool_or_int_flag("pack_side_chains", self.pack_side_chains)
            number_of_packs_per_design_flag = _val_flag("number_of_packs_per_design", self.number_of_packs_per_design)
            sc_num_denoising_steps_flag = _val_flag("sc_num_denoising_steps", self.sc_num_denoising_steps)
            sc_num_samples_flag = _val_flag("sc_num_samples", self.sc_num_samples)
            repack_everything_flag = _bool_or_int_flag("repack_everything", self.repack_everything)
            force_hetatm_flag = _bool_or_int_flag("force_hetatm", self.force_hetatm)
            packed_suffix_flag = _str_flag("packed_suffix", self.packed_suffix)
            pack_with_ligand_context_flag = _bool_or_int_flag("pack_with_ligand_context", self.pack_with_ligand_context)

            self._run_task(
                self.task_name,
                pdb_input_json=pdb_input_json.name,
                output_dir=self.output_dir,
                num_designs=self.num_designs,
                temperature=self.temperature,
                model_type_flag=model_type_flag,
                checkpoint_protein_mpnn_flag=checkpoint_protein_mpnn_flag,
                checkpoint_ligand_mpnn_flag=checkpoint_ligand_mpnn_flag,
                checkpoint_per_residue_label_membrane_mpnn_flag=checkpoint_per_residue_label_membrane_mpnn_flag,
                checkpoint_global_label_membrane_mpnn_flag=checkpoint_global_label_membrane_mpnn_flag,
                checkpoint_soluble_mpnn_flag=checkpoint_soluble_mpnn_flag,
                checkpoint_path_sc_flag=checkpoint_path_sc_flag,
                fasta_seq_separation_flag=fasta_seq_separation_flag,
                verbose_flag=verbose_flag,
                fixed_residues_flag=fixed_residues_flag,
                fixed_residues_multi_flag=fixed_residues_multi_flag,
                redesigned_residues_flag=redesigned_residues_flag,
                redesigned_residues_multi_flag=redesigned_residues_multi_flag,
                bias_aa_flag=bias_aa_flag,
                bias_aa_per_residue_flag=bias_aa_per_residue_flag,
                bias_aa_per_residue_multi_flag=bias_aa_per_residue_multi_flag,
                omit_aa_flag=omit_aa_flag,
                omit_aa_per_residue_flag=omit_aa_per_residue_flag,
                omit_aa_per_residue_multi_flag=omit_aa_per_residue_multi_flag,
                symmetry_residues_flag=symmetry_residues_flag,
                symmetry_weights_flag=symmetry_weights_flag,
                homo_oligomer_flag=homo_oligomer_flag,
                file_ending_flag=file_ending_flag,
                zero_indexed_flag=zero_indexed_flag,
                seed_flag=seed_flag,
                number_of_batches_flag=number_of_batches_flag,
                save_stats_flag=save_stats_flag,
                ligand_mpnn_use_atom_context_flag=ligand_mpnn_use_atom_context_flag,
                ligand_mpnn_cutoff_for_score_flag=ligand_mpnn_cutoff_for_score_flag,
                ligand_mpnn_use_side_chain_context_flag=ligand_mpnn_use_side_chain_context_flag,
                chains_to_design_flag=chains_to_design_flag,
                parse_these_chains_only_flag=parse_these_chains_only_flag,
                transmembrane_buried_flag=transmembrane_buried_flag,
                transmembrane_interface_flag=transmembrane_interface_flag,
                global_transmembrane_label_flag=global_transmembrane_label_flag,
                parse_atoms_with_zero_occupancy_flag=parse_atoms_with_zero_occupancy_flag,
                pack_side_chains_flag=pack_side_chains_flag,
                number_of_packs_per_design_flag=number_of_packs_per_design_flag,
                sc_num_denoising_steps_flag=sc_num_denoising_steps_flag,
                sc_num_samples_flag=sc_num_samples_flag,
                repack_everything_flag=repack_everything_flag,
                force_hetatm_flag=force_hetatm_flag,
                packed_suffix_flag=packed_suffix_flag,
                pack_with_ligand_context_flag=pack_with_ligand_context_flag,
                extra_args=self.extra_args,
                device=self.device,
            )
        finally:
            for tf in temp_files:
                try:
                    if os.path.exists(tf):
                        os.remove(tf)
                except Exception:
                    pass

        # Split the FASTA files:
        seqs_dir = self.output_dir / 'seqs'
        if seqs_dir.exists():
            for file in seqs_dir.iterdir():
                if file.is_file() and file.suffix in ['.fa', '.fasta']:
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

class RosettaScripts(Task):
    def __init__(self, output_dir: Union[str, Path], pdb_input_file: Optional[Union[str, Path]] = None, pdb_input_dir: Optional[Union[str, Path]] = None, protocol_xml: Optional[Union[str, Path]] = None, ligand_params_file: Optional[Union[str, Path]] = None, device: str = 'cpu', extra_args: str = "", dry_run: bool = False):
        super().__init__()
        self.task_name = "RosettaScripts"
        
        if pdb_input_file is None and pdb_input_dir is None:
            raise ValueError('Must specify either pdb_input_file or pdb_input_dir')
        self.output_dir = make_directory(output_dir)
        self.pdb_input_file = pdb_input_file
        self.pdb_input_dir = pdb_input_dir
        self.protocol_xml = protocol_xml
        self.ligand_params_file = ligand_params_file
        self.device = device
        self.extra_args = extra_args
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        if self.pdb_input_file is not None:
            pdb_list = [Path(self.pdb_input_file)]
        else:
            pdb_list = list_files(self.pdb_input_dir, '.pdb')
        
        import shutil
        for pdb_path in pdb_list:
            pdb_stem = Path(pdb_path).stem
            shutil.copy(pdb_path, self.output_dir / f"{pdb_stem}_0001.pdb")
        
        with open(self.output_dir / "score.sc", "w") as f:
            f.write("SEQUENCE:\nSCORE:     score pstat description\n")
            for pdb_path in pdb_list:
                pdb_stem = Path(pdb_path).stem
                f.write(f"SCORE:     -100.0 0.65 {pdb_stem}_0001\n")

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        if self.pdb_input_file is not None:
            temp_dir = tempfile.mkdtemp()
            shutil.copy(self.pdb_input_file, temp_dir)
            pdb_input_dir = temp_dir
        else:
            pdb_input_dir = self.pdb_input_dir

        pdb_list = list_files(pdb_input_dir, '.pdb')
        pdb_string = " ".join(map(str, pdb_list)) + " "

        params_flags = []
        if self.ligand_params_file:
            params_flags.append(f"-extra_res_fa {Path(self.ligand_params_file).resolve()}")

        dynamic_args = " ".join(params_flags)
        combined_extra_args = f"{dynamic_args} {self.extra_args}".strip() if self.extra_args else dynamic_args

        self._run_task(
            self.task_name,
            pdb_string=pdb_string,
            protocol_xml=str(Path(self.protocol_xml).resolve()),
            output_dir=str(self.output_dir),
            extra_args=combined_extra_args,
            device=self.device
        )

class LDDM(Task):
    """
    Ligand Denoising Diffusion Model (LDDM) for target-aware de novo molecule generation,
    fragment elaboration (growing/linking), and molecular docking.
    """

    def __init__(
        self,
        output_dir: Union[str, Path],
        protein: Optional[Union[str, Path]] = None,
        ref_ligand: Optional[Union[str, Path]] = None,
        mode: str = "design",
        ligand: Optional[Union[str, Path]] = None,
        checkpoint: Optional[Union[str, Path]] = None,
        output_file: Optional[Union[str, Path]] = None,
        n_samples: int = 10,
        batch_size: Optional[int] = None,
        n_steps: int = 100,
        sampler: str = "ForwardEuler",
        sampling_noise: float = 5.0,
        molecule_size: Optional[Union[int, str]] = None,
        atoms_to_dock: Optional[Union[str, List[int]]] = None,
        n_frames: Optional[int] = None,
        return_projected_final: bool = False,
        seed: Optional[int] = None,
        device: str = "gpu",
        verbose: bool = False,
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        """
        Initialize an LDDM task for de novo design, fragment elaboration, or docking.

        Args:
            output_dir (str or Path): Directory where output files will be written.
            protein (str or Path): Path to protein PDB structure.
            ref_ligand (str or Path): Path to reference ligand SDF defining the pocket.
            mode (str): Sampling mode, either 'design' or 'dock'. Default is 'design'.
                - 'design': Generate new molecules de novo, or grow/link fragments if ligand is provided.
                - 'dock': Predict the binding pose of the molecule provided in ligand.
            ligand (str or Path, optional):
                - In 'design' mode: optional SDF with fragment(s) to grow or link.
                - In 'dock' mode: SDF of molecule to dock (required).
            checkpoint (str or Path, optional): Path to model checkpoint .ckpt file.
                Defaults to the built-in '/app/models/checkpoints/lddm.ckpt'.
                Accepts 'CD+BB' or 'lddm_CDBB' to use the MIT-licensed checkpoint.
            output_file (str or Path, optional): Output SDF file path. Default is output_dir / 'samples.sdf'.
            n_samples (int): Number of molecules or poses to generate. Default is 10.
            batch_size (int, optional): Number of samples per forward pass. Default is n_samples.
            n_steps (int): Number of ODE integration steps. Default is 100.
            sampler (str): ODE sampler ('ForwardEuler' or 'HeunSampler'). Default is 'ForwardEuler'.
            sampling_noise (float): Stochastic noise scale during sampling. Default is 5.0.
            molecule_size (int or str, optional): (Design mode only) Target ligand atom count.
                Integer for exact size, 'uniform_<low>_<high>' for range, or None for pocket histogram.
            atoms_to_dock (list of int or str, optional): (Dock mode only) 0-indexed atom indices
                to dock for partial docking. If omitted, docks the whole molecule.
            n_frames (int, optional): If set, save sampling trajectory with this many frames.
            return_projected_final (bool): When saving trajectory, project each frame onto clean prediction.
            seed (int, optional): Random seed for reproducible sampling.
            device (str): Device to run inference on ('gpu', 'cpu', 'gpu_wsl'). Default is 'gpu'.
            verbose (bool): Enable debug logging. Default is False.
            extra_args (str): Additional CLI arguments passed to lddm-sample. Default is "".
            dry_run (bool): If True, create mock output files without invoking container. Default is False.

        Returns:
            None
            Outputs written to output_file (default: output_dir / 'samples.sdf').
        """
        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "LDDM"

        # Support aliases
        if protein is None:
            protein = kwargs.pop("protein_file", None) or kwargs.pop("pdb_file", None)
        if ref_ligand is None:
            ref_ligand = kwargs.pop("ref_ligand_file", None)
        if ligand is None:
            ligand = kwargs.pop("ligand_file", None) or kwargs.pop("fragment_file", None)
        if "num_samples" in kwargs:
            n_samples = kwargs.pop("num_samples")
        if "num_steps" in kwargs:
            n_steps = kwargs.pop("num_steps")

        if kwargs:
            raise TypeError(f"LDDM.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}")

        if protein is None:
            raise ValueError("Argument 'protein' (path to protein PDB) is required for LDDM.")
        if ref_ligand is None:
            raise ValueError("Argument 'ref_ligand' (path to reference ligand SDF) is required for LDDM.")

        if mode not in ["design", "dock"]:
            raise ValueError(f"Invalid mode '{mode}'. LDDM supports 'design' and 'dock'.")

        if mode == "dock" and ligand is None:
            raise ValueError("In 'dock' mode, 'ligand' (SDF of molecule to dock) is required.")

        self.output_dir = Path(output_dir)
        self.output_file = Path(output_file) if output_file is not None else self.output_dir / "samples.sdf"
        self.protein = protein
        self.ref_ligand = ref_ligand
        self.mode = mode
        self.ligand = ligand
        self.checkpoint = checkpoint
        self.n_samples = n_samples
        self.batch_size = batch_size
        self.n_steps = n_steps
        self.sampler = sampler
        self.sampling_noise = sampling_noise
        self.molecule_size = molecule_size
        self.atoms_to_dock = atoms_to_dock
        self.n_frames = n_frames
        self.return_projected_final = return_projected_final
        self.seed = seed
        self.device = device
        self.verbose = verbose
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)
        self.output_file = Path(self.output_file)
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

        # Write mock SDF molecules
        mock_entry = (
            "mock_ligand\n"
            "     RDKit          3D\n\n"
            "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
            "M  END\n"
            "$$$$\n"
        )
        with open(self.output_file, "w") as f:
            for _ in range(max(1, self.n_samples)):
                f.write(mock_entry)

        if self.n_frames is not None:
            pocket_pdb = self.output_file.with_suffix(".pdb")
            pocket_pdb.touch()

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        self.output_dir = make_directory(self.output_dir)
        self.output_file = Path(self.output_file).resolve()
        self.output_file.parent.mkdir(parents=True, exist_ok=True)

        protein_path = Path(self.protein).resolve()
        if not protein_path.exists():
            raise FileNotFoundError(f"Protein file not found: {protein_path}")

        ref_ligand_path = Path(self.ref_ligand).resolve()
        if not ref_ligand_path.exists():
            raise FileNotFoundError(f"Reference ligand file not found: {ref_ligand_path}")

        # Checkpoint resolution
        if self.checkpoint is None or str(self.checkpoint).lower() in ["lddm", "default"]:
            checkpoint_val = "/app/models/checkpoints/lddm.ckpt"
        elif str(self.checkpoint).lower() in ["cd+bb", "lddm_cdbb", "mit"]:
            checkpoint_val = "/app/models/checkpoints/lddm_CDBB.ckpt"
        else:
            chk = Path(self.checkpoint).resolve()
            if not chk.exists():
                raise FileNotFoundError(f"Checkpoint file not found: {chk}")
            checkpoint_val = str(chk)

        # Optional ligand
        ligand_flag = ""
        if self.ligand is not None:
            lig_path = Path(self.ligand).resolve()
            if not lig_path.exists():
                raise FileNotFoundError(f"Ligand file not found: {lig_path}")
            ligand_flag = f'--ligand "{lig_path}"'

        # Mode specific flags
        molecule_size_flag = ""
        if self.mode == "design" and self.molecule_size is not None:
            molecule_size_flag = f"--molecule_size {self.molecule_size}"

        atoms_to_dock_flag = ""
        if self.mode == "dock" and self.atoms_to_dock is not None:
            if isinstance(self.atoms_to_dock, (list, tuple)):
                indices_str = " ".join(str(x) for x in self.atoms_to_dock)
            else:
                indices_str = str(self.atoms_to_dock)
            atoms_to_dock_flag = f"--atoms_to_dock {indices_str}"

        # Sampling and trajectory flags
        batch_size_flag = f"--batch_size {self.batch_size}" if self.batch_size is not None else ""
        n_frames_flag = f"--n_frames {self.n_frames}" if self.n_frames is not None else ""
        return_projected_final_flag = "--return_projected_final" if self.return_projected_final else ""
        seed_flag = f"--seed {self.seed}" if self.seed is not None else ""
        verbose_flag = "--verbose" if self.verbose else ""

        # Device mapping: Ribbon runner accepts 'gpu', 'gpu_wsl', 'cpu'
        if self.device in ["gpu", "gpu_wsl", "cpu"]:
            ribbon_device = self.device
        elif "cuda" in str(self.device).lower():
            ribbon_device = "gpu"
        else:
            ribbon_device = "cpu"

        if self.device in ["gpu", "gpu_wsl"]:
            lddm_device = "cuda:0"
        elif self.device == "cpu":
            lddm_device = "cpu"
        else:
            lddm_device = self.device
        device_flag = f"--device {lddm_device}"

        self._run_task(
            self.task_name,
            mode=self.mode,
            protein=str(protein_path),
            ref_ligand=str(ref_ligand_path),
            checkpoint=checkpoint_val,
            output_file=str(self.output_file),
            n_samples=self.n_samples,
            n_steps=self.n_steps,
            sampler=self.sampler,
            sampling_noise=self.sampling_noise,
            batch_size_flag=batch_size_flag,
            molecule_size_flag=molecule_size_flag,
            ligand_flag=ligand_flag,
            atoms_to_dock_flag=atoms_to_dock_flag,
            n_frames_flag=n_frames_flag,
            return_projected_final_flag=return_projected_final_flag,
            seed_flag=seed_flag,
            device_flag=device_flag,
            verbose_flag=verbose_flag,
            extra_args=self.extra_args,
            device=ribbon_device,
        )

class BindCraft2(Task):
    """
    BindCraft2 (BC2): Protein binder design workflow combining AlphaFold 2 sequence optimisation,
    ProteinMPNN redesign, structural filtering, and multi-state evaluation for de novo miniproteins,
    peptides, cyclic peptides, and scaffolded modalities (VHH, ARP, Fab, scFv).
    """

    def __init__(
        self,
        output_dir: Union[str, Path],
        settings_file: Optional[Union[str, Path]] = None,
        target_pdb: Optional[Union[str, Path]] = None,
        target_name: Optional[str] = None,
        chains: Optional[Union[str, List[str]]] = None,
        hotspots: Optional[Union[str, List[Union[int, str]]]] = None,
        target: Optional[str] = None,
        targets: Optional[List[Dict[str, Any]]] = None,
        modality: str = "binder",
        binder_lengths: Optional[Union[int, List[int], Tuple[int, int], str]] = None,
        binder_chain: Optional[str] = None,
        copies: Optional[int] = None,
        number_of_final_designs: int = 1,
        max_trajectories: Optional[int] = None,
        campaign_seed: Optional[int] = None,
        campaign_name: Optional[str] = None,
        core: Optional[str] = None,
        bigbang: bool = False,
        disulfide_staple: bool = False,
        forced_targeting: bool = False,
        humanize: bool = False,
        initial_guess: bool = False,
        mixed_topology: bool = False,
        protease_stable: bool = False,
        termini_accessible: bool = False,
        termini_together: bool = False,
        screen_steps: Optional[int] = None,
        refine_steps: Optional[int] = None,
        anneal_steps: Optional[int] = None,
        harden_steps: Optional[int] = None,
        mutate_steps: Optional[int] = None,
        design_recycles: Optional[int] = None,
        validation_recycles: Optional[int] = None,
        design_models: Optional[Union[str, List[str]]] = None,
        validation_models: Optional[Union[str, List[str]]] = None,
        mpnn_variant: Optional[str] = None,
        resume: bool = False,
        metadata_file: Optional[Union[str, Path]] = None,
        settings_overrides: Optional[Dict[str, Any]] = None,
        device: str = "gpu",
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        """
        Initialize a BindCraft2 binder design task.

        Args:
            output_dir (str or Path): Campaign output directory (project_folder).
            settings_file (str or Path, optional): Path to existing campaign settings JSON file.
            target_pdb (str or Path, optional): Path to target PDB structure.
            target_name (str, optional): Target name (defaults to PDB stem).
            chains (str or list of str, optional): Target chain(s) to design against (e.g. 'A' or ['A', 'B']).
            hotspots (str or list, optional): Hotspot residues on target (e.g. '54,56,66,115' or [54, 56, 66]).
            target (str, optional): Shipped target name if not using a custom PDB (e.g. 'hPDL1', 'hIL2R', etc.).
            targets (list of dict, optional): List of target definitions for multi-target campaigns.
            modality (str): Binder format / preset. Default is 'binder'.
                Options: 'binder', 'peptide', 'cyclic_peptide', 'ARP', 'Fab', 'VHH', 'scFv',
                'homo_oligomer', 'multidomain', 'induced_fit', 'fold_switch', 'large_binder'.
            binder_lengths (int, list, tuple, or str, optional): Target binder length or range, e.g. [60, 100] or 80.
            binder_chain (str, optional): Chain ID for the binder. Default is 'B'.
            copies (int, optional): Number of binder copies (for homo-oligomers).
            number_of_final_designs (int): Number of accepted designs to generate. Default is 1.
            max_trajectories (int, optional): Maximum design trajectories to attempt before giving up.
            campaign_seed (int, optional): Random seed for reproducible trajectories.
            campaign_name (str, optional): Campaign name. Defaults to output_dir name.
            core (str, optional): Core profile under presets (e.g. 'benchmark').
            bigbang (bool): Start gradient stages from coordinates (--bigbang).
            disulfide_staple (bool): Require geometrically detected disulfide (--disulfide-staple).
            forced_targeting (bool): Concentrate paratope on target hotspots (--forced-targeting).
            humanize (bool): Human-like sequence preferences and MHC-II anchor score (--humanize).
            initial_guess (bool): Re-predict accepted candidate from folded trajectory (--initial-guess).
            mixed_topology (bool): Enforce secondary structure mix (<=50% helix, >=20% sheet) (--mixed-topology).
            protease_stable (bool): Penalize protease cleavage motifs and exposed loops (--protease-stable).
            termini_accessible (bool): Point chain ends away from target (--termini-accessible).
            termini_together (bool): Bring N and C termini within 10A (--termini-together).
            screen_steps (int, optional): Steps in screening stage.
            refine_steps (int, optional): Steps in refinement stage.
            anneal_steps (int, optional): Steps in annealing stage.
            harden_steps (int, optional): Steps in hardening stage.
            mutate_steps (int, optional): Steps in mutation stage.
            design_recycles (int, optional): AF2 recycles during design.
            validation_recycles (int, optional): AF2 recycles during validation.
            design_models (str or list of str, optional): AF2 design models.
            validation_models (str or list of str, optional): AF2 validation models.
            mpnn_variant (str, optional): ProteinMPNN variant (e.g. 'v_48_020').
            resume (bool): Resume into an existing project folder.
            metadata_file (str or Path, optional): Metadata JSON file (--metadata <file>).
            settings_overrides (dict, optional): Additional settings overrides passed via --set KEY=VALUE.
            device (str): Device to run on ('gpu', 'cpu', 'gpu_wsl'). Default is 'gpu'.
            extra_args (str): Additional raw CLI flags. Default is "".
            dry_run (bool): If True, create mock outputs without running container. Default is False.
        """
        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "BindCraft2"

        # Support aliases
        if target_pdb is None:
            target_pdb = (
                kwargs.pop("pdb_file", None)
                or kwargs.pop("protein_file", None)
                or kwargs.pop("protein", None)
                or kwargs.pop("target_path", None)
            )
        if chains is None:
            chains = kwargs.pop("target_chains", None)
        if hotspots is None:
            hotspots = kwargs.pop("target_hotspots", None)
        if binder_lengths is None:
            binder_lengths = (
                kwargs.pop("binder_length", None)
                or kwargs.pop("lengths", None)
                or kwargs.pop("length", None)
            )
        if "num_designs" in kwargs:
            number_of_final_designs = kwargs.pop("num_designs")
        elif "n_designs" in kwargs:
            number_of_final_designs = kwargs.pop("n_designs")
        elif "final_designs" in kwargs:
            number_of_final_designs = kwargs.pop("final_designs")
        if max_trajectories is None:
            max_trajectories = kwargs.pop("max_attempts", None)
        if campaign_seed is None:
            campaign_seed = kwargs.pop("seed", None)
        if settings_overrides is None:
            settings_overrides = kwargs.pop("set_settings", None)

        if kwargs:
            raise TypeError(f"BindCraft2.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}")

        if (
            settings_file is None
            and target_pdb is None
            and target is None
            and targets is None
        ):
            raise ValueError(
                "Must specify either 'settings_file', 'target_pdb' (path to PDB), "
                "'target' (shipped target name e.g. 'hPDL1'), or 'targets' (list of target dicts)."
            )

        valid_modalities = [
            "binder",
            "peptide",
            "cyclic_peptide",
            "ARP",
            "Fab",
            "VHH",
            "scFv",
            "homo_oligomer",
            "multidomain",
            "induced_fit",
            "fold_switch",
            "large_binder",
        ]
        if modality not in valid_modalities:
            raise ValueError(f"Invalid modality '{modality}'. Expected one of: {valid_modalities}")

        self.output_dir = Path(output_dir)
        self.settings_file = Path(settings_file) if settings_file is not None else None
        self.target_pdb = target_pdb
        self.target_name = target_name
        self.chains = ",".join(chains) if isinstance(chains, (list, tuple)) else chains
        if isinstance(hotspots, (list, tuple)):
            self.hotspots = ",".join(str(h) for h in hotspots)
        else:
            self.hotspots = str(hotspots) if hotspots is not None else None
        self.target = target
        self.targets = targets
        self.modality = modality

        # Format binder_lengths
        if isinstance(binder_lengths, int):
            self.binder_lengths = [binder_lengths, binder_lengths]
        elif isinstance(binder_lengths, (list, tuple)):
            self.binder_lengths = list(binder_lengths)
        elif isinstance(binder_lengths, str):
            self.binder_lengths = [int(x.strip()) for x in binder_lengths.split(",") if x.strip()]
        else:
            self.binder_lengths = None

        self.binder_chain = binder_chain
        self.copies = copies
        self.number_of_final_designs = number_of_final_designs
        self.max_trajectories = max_trajectories
        self.campaign_seed = campaign_seed
        self.campaign_name = campaign_name
        self.core = core

        # Properties
        self.bigbang = bigbang
        self.disulfide_staple = disulfide_staple
        self.forced_targeting = forced_targeting
        self.humanize = humanize
        self.initial_guess = initial_guess
        self.mixed_topology = mixed_topology
        self.protease_stable = protease_stable
        self.termini_accessible = termini_accessible
        self.termini_together = termini_together

        # Optimization knobs
        self.screen_steps = screen_steps
        self.refine_steps = refine_steps
        self.anneal_steps = anneal_steps
        self.harden_steps = harden_steps
        self.mutate_steps = mutate_steps
        self.design_recycles = design_recycles
        self.validation_recycles = validation_recycles
        self.design_models = design_models
        self.validation_models = validation_models
        self.mpnn_variant = mpnn_variant
        self.resume = resume

        self.metadata_file = metadata_file
        self.settings_overrides = settings_overrides or {}
        self.device = device
        self.extra_args = extra_args
        self.dry_run = dry_run

    def _run_dry(self):
        self.output_dir = make_directory(self.output_dir)

        # 3_Ranked
        ranked_dir = self.output_dir / "3_Ranked"
        ranked_dir.mkdir(parents=True, exist_ok=True)
        ranked_csv = ranked_dir / "!_Ranked.csv"
        with open(ranked_csv, "w") as f:
            f.write("rank,design,hash,trajectory,binder_length,i_pDAE,i_pTM,pLDDT,outcome\n")
            for i in range(1, max(1, self.number_of_final_designs) + 1):
                f.write(f"{i},mock_design_{i:04d},hash_{i:04d},{i},80,0.85,0.82,0.88,passed\n")
                (ranked_dir / f"mock_design_{i:04d}_seq1.cif").touch()

        # 2_Refolded
        refolded_dir = self.output_dir / "2_Refolded"
        refolded_dir.mkdir(parents=True, exist_ok=True)
        (refolded_dir / "!_Refolded.csv").touch()

        # 1_Trajectories
        traj_dir = self.output_dir / "1_Trajectories"
        traj_dir.mkdir(parents=True, exist_ok=True)
        (traj_dir / "!_Trajectories.csv").touch()

        # summary.csv
        summary_csv = self.output_dir / "summary.csv"
        with open(summary_csv, "w") as f:
            f.write("campaign,scope,metric,samples,mean,std,min,max\n")
            f.write(f"{self.output_dir.name},campaign,accepted_designs,1,{self.number_of_final_designs},0,{self.number_of_final_designs},{self.number_of_final_designs}\n")

        # campaign_metadata.json
        meta_json = self.output_dir / "campaign_metadata.json"
        with open(meta_json, "w") as f:
            json.dump({"software": "BindCraft2", "designs": self.number_of_final_designs}, f, indent=2)

    def run(self):
        if self.dry_run:
            self._run_dry()
            return

        self.output_dir = make_directory(self.output_dir)

        if self.settings_file is not None:
            settings_path = Path(self.settings_file).resolve()
            if not settings_path.exists():
                raise FileNotFoundError(f"Settings file not found: {settings_path}")
        else:
            settings_path = self.output_dir / "campaign_settings.json"
            cfg: Dict[str, Any] = {
                "project_folder": str(self.output_dir.resolve()),
                "modality": self.modality,
                "number_of_final_designs": self.number_of_final_designs,
            }
            if self.campaign_name:
                cfg["campaign_name"] = self.campaign_name
            else:
                cfg["campaign_name"] = self.output_dir.name or "bindcraft_campaign"

            if self.binder_lengths is not None:
                cfg["binder_lengths"] = self.binder_lengths
            if self.max_trajectories is not None:
                cfg["max_trajectories"] = self.max_trajectories
            if self.campaign_seed is not None:
                cfg["campaign_seed"] = self.campaign_seed
            if self.binder_chain is not None:
                cfg["binder_chain"] = self.binder_chain
            if self.copies is not None:
                cfg["copies"] = self.copies
            if self.screen_steps is not None:
                cfg["screen_steps"] = self.screen_steps
            if self.refine_steps is not None:
                cfg["refine_steps"] = self.refine_steps
            if self.anneal_steps is not None:
                cfg["anneal_steps"] = self.anneal_steps
            if self.harden_steps is not None:
                cfg["harden_steps"] = self.harden_steps
            if self.mutate_steps is not None:
                cfg["mutate_steps"] = self.mutate_steps
            if self.design_recycles is not None:
                cfg["design_recycles"] = self.design_recycles
            if self.validation_recycles is not None:
                cfg["validation_recycles"] = self.validation_recycles
            if self.design_models is not None:
                cfg["design_models"] = self.design_models
            if self.validation_models is not None:
                cfg["validation_models"] = self.validation_models
            if self.mpnn_variant is not None:
                cfg["mpnn_variant"] = self.mpnn_variant
            if self.resume:
                cfg["resume"] = True

            if self.targets is not None:
                cfg["targets"] = self.targets
            elif self.target is not None:
                cfg["target"] = self.target
            elif self.target_pdb is not None:
                t_path = Path(self.target_pdb).resolve()
                if not t_path.exists():
                    raise FileNotFoundError(f"Target PDB file not found: {t_path}")
                t_dict = {
                    "name": self.target_name or t_path.stem,
                    "target_path": str(t_path),
                }
                if self.chains is not None:
                    t_dict["chains"] = self.chains
                if self.hotspots is not None:
                    t_dict["hotspots"] = self.hotspots
                cfg["targets"] = [t_dict]

            if self.settings_overrides:
                cfg.update(self.settings_overrides)

            with open(settings_path, "w") as f:
                json.dump(cfg, f, indent=2)

        # Set flags (passed to --set KEY=VALUE)
        set_list = []
        if self.settings_file is not None:
            set_list.append(f"--set 'project_folder={self.output_dir.resolve()}'")
            if self.number_of_final_designs is not None:
                set_list.append(f"--set 'number_of_final_designs={self.number_of_final_designs}'")
            if self.binder_lengths is not None:
                set_list.append(f"--set 'binder_lengths={json.dumps(self.binder_lengths)}'")
            if self.max_trajectories is not None:
                set_list.append(f"--set 'max_trajectories={self.max_trajectories}'")
            if self.campaign_seed is not None:
                set_list.append(f"--set 'campaign_seed={self.campaign_seed}'")
            if self.settings_overrides:
                for k, v in self.settings_overrides.items():
                    val_str = json.dumps(v) if isinstance(v, (dict, list, bool)) else str(v)
                    set_list.append(f"--set '{k}={val_str}'")
        set_flags = " ".join(set_list)

        core_flag = f"--core {self.core}" if self.core else ""
        modality_flag = f"--modality {self.modality}" if self.modality else ""

        props = []
        if self.bigbang: props.append("--bigbang")
        if self.disulfide_staple: props.append("--disulfide-staple")
        if self.forced_targeting: props.append("--forced-targeting")
        if self.humanize: props.append("--humanize")
        if self.initial_guess: props.append("--initial-guess")
        if self.mixed_topology: props.append("--mixed-topology")
        if self.protease_stable: props.append("--protease-stable")
        if self.termini_accessible: props.append("--termini-accessible")
        if self.termini_together: props.append("--termini-together")
        property_flags = " ".join(props)

        metadata_flag = f'--metadata "{Path(self.metadata_file).resolve()}"' if self.metadata_file else ""

        ribbon_device = self.device if self.device in ["gpu", "gpu_wsl", "cpu"] else "gpu"

        self._run_task(
            self.task_name,
            settings_file=str(settings_path),
            core_flag=core_flag,
            modality_flag=modality_flag,
            property_flags=property_flags,
            set_flags=set_flags,
            metadata_flag=metadata_flag,
            extra_args=self.extra_args,
            device=ribbon_device,
        )


class PyMolScript(Task):
    """
    Run a Python script inside the PyMol container.

    Enables running arbitrary Python scripts that interact with PyMOL (e.g. using
    the `pymol` and `pymol.cmd` modules) within the containerized PyMOL environment.

    Args:
        script_file (Union[str, Path]): Path to the Python (.py) script to execute.
        script_args (Optional[Union[str, List[str], Tuple[str, ...]]]): Command-line arguments
            to pass to the Python script. Can be a formatted string or a sequence of arguments.
        output_dir (Optional[Union[str, Path]]): Optional output directory to ensure exists.
        output_file (Optional[Union[str, Path]]): Optional output file whose parent directory
            will be created if specified.
        device (str): Device to run the task on ('cpu', 'gpu', or 'gpu_wsl'). Default is 'cpu'.
        extra_args (str): Additional CLI flags to append to the command. Default is ''.
        dry_run (bool): If True, runs mock simulation without invoking the container.
    """

    def __init__(
        self,
        script_file: Optional[Union[str, Path]] = None,
        script_args: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
        output_dir: Optional[Union[str, Path]] = None,
        output_file: Optional[Union[str, Path]] = None,
        device: str = "cpu",
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "PyMolScript"

        # 1. Alias Handling
        if script_file is None:
            script_file = (
                kwargs.pop("script", None)
                or kwargs.pop("python_script", None)
                or kwargs.pop("script_path", None)
                or kwargs.pop("input_file", None)
                or kwargs.pop("file", None)
            )
        if script_args is None:
            script_args = (
                kwargs.pop("args", None)
                or kwargs.pop("arguments", None)
                or kwargs.pop("params", None)
            )
        if output_dir is None:
            output_dir = kwargs.pop("out_dir", None) or kwargs.pop("output_folder", None)
        if output_file is None:
            output_file = kwargs.pop("out_file", None)
        if "device_type" in kwargs:
            device = kwargs.pop("device_type")

        # 2. Reject Unexpected Arguments
        if kwargs:
            raise TypeError(
                f"{self.__class__.__name__}.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}"
            )

        # 3. Argument Validation
        if script_file is None:
            raise ValueError(f"Argument 'script_file' is required for {self.task_name}.")
        if device not in ["cpu", "gpu", "gpu_wsl"]:
            raise ValueError(f"Invalid device '{device}'. Supported: 'cpu', 'gpu', 'gpu_wsl'.")

        self.script_file = Path(script_file)
        self.script_args = script_args
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self.output_file = Path(output_file) if output_file is not None else None
        self.device = device
        self.extra_args = extra_args
        self.dry_run = dry_run

    def _compile_script_args(self) -> str:
        if self.script_args is None:
            return ""
        if isinstance(self.script_args, (list, tuple)):
            return " ".join(str(a) for a in self.script_args)
        return str(self.script_args).strip()

    def _run_dry(self):
        """Generate realistic mock outputs for dry-run validation."""
        print(f"[DRY RUN] {self.task_name} running script: {self.script_file}")
        if self.output_dir is not None:
            self.output_dir = make_directory(self.output_dir)
            mock_out = self.output_dir / "output.txt"
            if not mock_out.exists():
                with open(mock_out, "w") as f:
                    f.write(f"Mock output from dry run of {self.task_name}\n")
        if self.output_file is not None:
            self.output_file = Path(self.output_file)
            make_directory(self.output_file.parent)
            with open(self.output_file, "w") as f:
                f.write(f"Mock output from dry run of {self.task_name}\n")

    def run(self):
        """Execute task or simulate dry run."""
        if self.dry_run:
            self._run_dry()
            return

        # 1. Verify and resolve file paths
        script_path = self.script_file.resolve()
        if not script_path.exists():
            raise FileNotFoundError(f"Script file not found: {script_path}")

        if self.output_dir is not None:
            self.output_dir = make_directory(Path(self.output_dir).resolve())

        if self.output_file is not None:
            self.output_file = Path(self.output_file).resolve()
            make_directory(self.output_file.parent)

        # 2. Compile CLI arguments
        script_args_str = self._compile_script_args()

        # 3. Device translation
        ribbon_device = self.device if self.device in ["gpu", "gpu_wsl", "cpu"] else "cpu"

        # 4. Invoke Task runner
        self._run_task(
            self.task_name,
            script_file=str(script_path),
            script_args=script_args_str,
            extra_args=self.extra_args,
            device=ribbon_device,
        )


PyMOLScript = PyMolScript
RunPyMolScript = PyMolScript


class BioPythonScript(Task):
    """
    Run a Python script inside the BioPython container.

    Enables running arbitrary Python scripts that interact with BioPython (e.g. using
    the `Bio` package and related bioinformatic tools) within the containerized environment.

    Args:
        script_file (Union[str, Path]): Path to the Python (.py) script to execute.
        script_args (Optional[Union[str, List[str], Tuple[str, ...]]]): Command-line arguments
            to pass to the Python script. Can be a formatted string or a sequence of arguments.
        output_dir (Optional[Union[str, Path]]): Optional output directory to ensure exists.
        output_file (Optional[Union[str, Path]]): Optional output file whose parent directory
            will be created if specified.
        device (str): Device to run the task on ('cpu', 'gpu', or 'gpu_wsl'). Default is 'cpu'.
        extra_args (str): Additional CLI flags to append to the command. Default is ''.
        dry_run (bool): If True, runs mock simulation without invoking the container.
    """

    def __init__(
        self,
        script_file: Optional[Union[str, Path]] = None,
        script_args: Optional[Union[str, List[str], Tuple[str, ...]]] = None,
        output_dir: Optional[Union[str, Path]] = None,
        output_file: Optional[Union[str, Path]] = None,
        device: str = "cpu",
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "BioPythonScript"

        # 1. Alias Handling
        if script_file is None:
            script_file = (
                kwargs.pop("script", None)
                or kwargs.pop("python_script", None)
                or kwargs.pop("script_path", None)
                or kwargs.pop("input_file", None)
                or kwargs.pop("file", None)
            )
        if script_args is None:
            script_args = (
                kwargs.pop("args", None)
                or kwargs.pop("arguments", None)
                or kwargs.pop("params", None)
            )
        if output_dir is None:
            output_dir = kwargs.pop("out_dir", None) or kwargs.pop("output_folder", None)
        if output_file is None:
            output_file = kwargs.pop("out_file", None)
        if "device_type" in kwargs:
            device = kwargs.pop("device_type")

        # 2. Reject Unexpected Arguments
        if kwargs:
            raise TypeError(
                f"{self.__class__.__name__}.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}"
            )

        # 3. Argument Validation
        if script_file is None:
            raise ValueError(f"Argument 'script_file' is required for {self.task_name}.")
        if device not in ["cpu", "gpu", "gpu_wsl"]:
            raise ValueError(f"Invalid device '{device}'. Supported: 'cpu', 'gpu', 'gpu_wsl'.")

        self.script_file = Path(script_file)
        self.script_args = script_args
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self.output_file = Path(output_file) if output_file is not None else None
        self.device = device
        self.extra_args = extra_args
        self.dry_run = dry_run

    def _compile_script_args(self) -> str:
        if self.script_args is None:
            return ""
        if isinstance(self.script_args, (list, tuple)):
            return " ".join(str(a) for a in self.script_args)
        return str(self.script_args).strip()

    def _run_dry(self):
        """Generate realistic mock outputs for dry-run validation."""
        print(f"[DRY RUN] {self.task_name} running script: {self.script_file}")
        if self.output_dir is not None:
            self.output_dir = make_directory(self.output_dir)
            mock_out = self.output_dir / "output.txt"
            if not mock_out.exists():
                with open(mock_out, "w") as f:
                    f.write(f"Mock output from dry run of {self.task_name}\n")
        if self.output_file is not None:
            self.output_file = Path(self.output_file)
            make_directory(self.output_file.parent)
            with open(self.output_file, "w") as f:
                f.write(f"Mock output from dry run of {self.task_name}\n")

    def run(self):
        """Execute task or simulate dry run."""
        if self.dry_run:
            self._run_dry()
            return

        # 1. Verify and resolve file paths
        script_path = self.script_file.resolve()
        if not script_path.exists():
            raise FileNotFoundError(f"Script file not found: {script_path}")

        if self.output_dir is not None:
            self.output_dir = make_directory(Path(self.output_dir).resolve())

        if self.output_file is not None:
            self.output_file = Path(self.output_file).resolve()
            make_directory(self.output_file.parent)

        # 2. Compile CLI arguments
        script_args_str = self._compile_script_args()

        # 3. Device translation
        ribbon_device = self.device if self.device in ["gpu", "gpu_wsl", "cpu"] else "cpu"

        # 4. Invoke Task runner
        self._run_task(
            self.task_name,
            script_file=str(script_path),
            script_args=script_args_str,
            extra_args=self.extra_args,
            device=ribbon_device,
        )


BiopythonScript = BioPythonScript
RunBioPythonScript = BioPythonScript


class PLACER(Task):
    """
    Protein-Ligand Atomistic Conformational Ensemble Resolver (PLACER).

    PLACER is an all-atom graph neural network that denoises atomic coordinates to predict
    ligand binding poses, dock ligands into protein pockets, predict protein sidechain
    conformations, and model conformational heterogeneity.
    """

    def __init__(
        self,
        output_dir: Union[str, Path] = ".",
        input_file: Optional[Union[str, Path]] = None,
        input_dir: Optional[Union[str, Path]] = None,
        input_files: Optional[List[Union[str, Path]]] = None,
        num_samples: int = 10,
        predict_ligand: Optional[Union[str, List[str]]] = None,
        fixed_ligand: Optional[Union[str, List[str]]] = None,
        predict_multi: bool = False,
        target_res: Optional[str] = None,
        ligand_file: Optional[Union[str, List[str], Dict[str, Union[str, Path]]]] = None,
        rerank: Optional[str] = None,
        bonds: Optional[Union[str, List[str]]] = None,
        mutate: Optional[Union[str, List[str], Dict[str, str]]] = None,
        crop_centers: Optional[Union[str, List[str]]] = None,
        corruption_centers: Optional[Union[str, List[str]]] = None,
        residue_json: Optional[Union[str, Path]] = None,
        fixed_ligand_noise: Optional[float] = None,
        exclude_common_ligands: bool = False,
        ignore_ligand_hydrogens: bool = False,
        use_sm: bool = True,
        suffix: Optional[str] = None,
        cautious: bool = False,
        weights: Optional[Union[str, Path]] = None,
        device: str = "gpu",
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        """
        Initialize a PLACER task.

        Args:
            output_dir (str or Path): Output folder to save models and CSV files. Default is '.'.
            input_file (str or Path, optional): Path to single input PDB or mmCIF file (or list file).
            input_dir (str or Path, optional): Directory containing input PDB files.
            input_files (list of str or Path, optional): List of input PDB/mmCIF files.
            num_samples (int): Number of samples/models to generate per structure. Default is 10.
            predict_ligand (str or list of str, optional): Ligand identifier(s) to predict (e.g. 'D-LDP-501').
            fixed_ligand (str or list of str, optional): Ligand identifier(s) to keep fixed.
            predict_multi (bool): Predict and score all allowed ligands in input simultaneously. Default False.
            target_res (str, optional): Residue used as crop center (e.g. 'A-149' or 'A-149-TRP').
            ligand_file (str, list of str, or dict, optional): Reference ligand definition(s)
                (e.g. 'HEM:ligand.mol2', {'HEM': 'ligand.mol2'}, or 'HEM:CCD').
            rerank (str, optional): Metric to rank output models and CSV ('prmsd', 'plddt', 'plddt_pde').
            bonds (str or list of str, optional): Space-separated list or list of bonds
                (e.g. 'A-42-ALA-CB:B-173-JRP-CL:1.5').
            mutate (str, list of str, or dict, optional): Positions to mutate
                (e.g. '128A:75I', ['128A:75I', '5A:TRP'], or {'128A': '75I'}).
            crop_centers (str or list of str, optional): Atom names used as crop centers.
            corruption_centers (str or list of str, optional): Atom names used as corruption centers.
            residue_json (str or Path, optional): JSON file specifying custom residues to add to CCD library.
            fixed_ligand_noise (float, optional): Noise added to fixed ligand coordinates.
            exclude_common_ligands (bool): Exclude common solvents and crystallographic additives. Default False.
            ignore_ligand_hydrogens (bool): Ignore hydrogen mismatch in ligand definitions. Default False.
            use_sm (bool): Whether to make predictions with small molecule (True for holo, False for apo). Default True.
            suffix (str, optional): Suffix added to output PDB and CSV files.
            cautious (bool): Skip prediction if output CSV already exists. Default False.
            weights (str or Path, optional): Custom model weights .pt file.
            device (str): Execution device ('gpu', 'cpu', 'gpu_wsl'). Default is 'gpu'.
            extra_args (str): Additional commandline arguments passed to run_PLACER.py.
            dry_run (bool): If True, create mock output files without invoking container. Default False.
        """
        # Support aliases
        if input_file is None:
            input_file = (
                kwargs.pop("pdb_file", None)
                or kwargs.pop("pdb_input_file", None)
                or kwargs.pop("input_structure", None)
                or kwargs.pop("cif_file", None)
                or kwargs.pop("ifile", None)
            )
        if input_dir is None:
            input_dir = kwargs.pop("idir", None)
        if "nsamples" in kwargs:
            num_samples = kwargs.pop("nsamples")
        if predict_ligand is None and "predict_ligands" in kwargs:
            predict_ligand = kwargs.pop("predict_ligands")
        if fixed_ligand is None and "fixed_ligands" in kwargs:
            fixed_ligand = kwargs.pop("fixed_ligands")
        if ligand_file is None and "ligand_files" in kwargs:
            ligand_file = kwargs.pop("ligand_files")

        if kwargs:
            raise TypeError(f"PLACER.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}")

        input_count = sum(x is not None for x in [input_file, input_dir, input_files])
        if input_count != 1:
            raise ValueError("Provide exactly one of: input_file, input_dir, input_files")

        if num_samples < 1:
            raise ValueError(f"num_samples must be >= 1, got {num_samples}")

        if rerank is not None and rerank not in ["prmsd", "plddt", "plddt_pde"]:
            raise ValueError(f"Invalid rerank scoreterm '{rerank}'. Must be one of: 'prmsd', 'plddt', 'plddt_pde'")

        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "PLACER"

        self.output_dir = Path(output_dir)
        self.input_file = Path(input_file) if input_file is not None else None
        self.input_dir = Path(input_dir) if input_dir is not None else None
        self.input_files = [Path(f) for f in input_files] if input_files is not None else None
        self.num_samples = num_samples
        self.predict_ligand = predict_ligand
        self.fixed_ligand = fixed_ligand
        self.predict_multi = predict_multi
        self.target_res = target_res
        self.ligand_file = ligand_file
        self.rerank = rerank
        self.bonds = bonds
        self.mutate = mutate
        self.crop_centers = crop_centers
        self.corruption_centers = corruption_centers
        self.residue_json = Path(residue_json) if residue_json is not None else None
        self.fixed_ligand_noise = fixed_ligand_noise
        self.exclude_common_ligands = exclude_common_ligands
        self.ignore_ligand_hydrogens = ignore_ligand_hydrogens
        self.use_sm = use_sm
        self.suffix = suffix
        self.cautious = cautious
        self.weights = Path(weights) if weights is not None else None
        self.device = device
        self.dry_run = dry_run

    def _get_input_labels(self) -> List[str]:
        """Determine base labels for output file naming."""
        labels = []
        if self.input_file is not None:
            name = self.input_file.name
            for ext in [".cif.gz", ".cif", ".pdb", ".ent"]:
                if name.endswith(ext):
                    name = name[:-len(ext)]
                    break
            labels.append(f"{name}_{self.suffix}" if self.suffix else name)
        elif self.input_files is not None:
            for f in self.input_files:
                name = Path(f).name
                for ext in [".cif.gz", ".cif", ".pdb", ".ent"]:
                    if name.endswith(ext):
                        name = name[:-len(ext)]
                        break
                labels.append(f"{name}_{self.suffix}" if self.suffix else name)
        elif self.input_dir is not None:
            if self.input_dir.exists():
                for f in sorted(self.input_dir.glob("*.pdb")):
                    name = f.stem
                    labels.append(f"{name}_{self.suffix}" if self.suffix else name)
            if not labels:
                labels.append(f"model_{self.suffix}" if self.suffix else "model")
        return labels

    @property
    def output_csv(self) -> Optional[Path]:
        """Return primary output CSV score file path (if single input structure)."""
        labels = self._get_input_labels()
        if labels:
            return self.output_dir / f"{labels[0]}.csv"
        return None

    @property
    def output_model_pdb(self) -> Optional[Path]:
        """Return primary output multimodel PDB file path (if single input structure)."""
        labels = self._get_input_labels()
        if labels:
            return self.output_dir / f"{labels[0]}_model.pdb"
        return None

    def _run_dry(self):
        """Generate realistic mock outputs for dry-run validation."""
        self.output_dir = make_directory(self.output_dir)
        labels = self._get_input_labels()
        for label in labels:
            csv_path = self.output_dir / f"{label}.csv"
            pdb_path = self.output_dir / f"{label}_model.pdb"

            with open(csv_path, "w") as f:
                f.write("label,model_idx,fape,lddt,rmsd,kabsch,prmsd,plddt,plddt_pde\n")
                for i in range(1, self.num_samples + 1):
                    f.write(f"{label},{i},2.177,0.852,7.594,1.151,6.062,0.468,0.668\n")

            with open(pdb_path, "w") as f:
                for i in range(1, self.num_samples + 1):
                    f.write(f"MODEL {i}\n")
                    f.write("ATOM      1  N   SER A   1       0.000   0.000   0.000  1.00  0.15           N\n")
                    f.write("ATOM      2  CA  SER A   1       1.458   0.000   0.000  1.00  0.12           C\n")
                    f.write("ATOM      3  C   SER A   1       2.000   1.400   0.000  1.00  0.11           C\n")
                    f.write("ATOM      4  O   SER A   1       1.300   2.400   0.000  1.00  0.23           O\n")
                    f.write("TER\n")
                    f.write("ENDMDL\n")

    def run(self):
        """Execute PLACER inference or generate mock files in dry_run mode."""
        if self.dry_run:
            self._run_dry()
            return

        self.output_dir = make_directory(self.output_dir.resolve())

        # 1. Resolve and validate input
        if self.input_file is not None:
            inp = self.input_file.resolve()
            if not inp.exists():
                raise FileNotFoundError(f"Input file not found: {inp}")
            input_flag = f'-f "{inp}"'
        elif self.input_dir is not None:
            inp_dir = self.input_dir.resolve()
            if not inp_dir.exists():
                raise FileNotFoundError(f"Input directory not found: {inp_dir}")
            input_flag = f'-i "{inp_dir}"'
        elif self.input_files is not None:
            resolved_files = []
            for f in self.input_files:
                p = Path(f).resolve()
                if not p.exists():
                    raise FileNotFoundError(f"Input file not found: {p}")
                resolved_files.append(str(p))
            if len(resolved_files) == 1:
                input_flag = f'-f "{resolved_files[0]}"'
            else:
                list_file = self.output_dir / "placer_input_list.txt"
                with open(list_file, "w") as f:
                    for rf in resolved_files:
                        f.write(f"{rf}\n")
                input_flag = f'-f "{list_file}"'

        # 2. Format optional arguments
        predict_ligand_flag = ""
        if self.predict_ligand is not None:
            if isinstance(self.predict_ligand, (list, tuple)):
                items = " ".join(str(x) for x in self.predict_ligand)
            else:
                items = str(self.predict_ligand)
            predict_ligand_flag = f"--predict_ligand {items}"

        fixed_ligand_flag = ""
        if self.fixed_ligand is not None:
            if isinstance(self.fixed_ligand, (list, tuple)):
                items = " ".join(str(x) for x in self.fixed_ligand)
            else:
                items = str(self.fixed_ligand)
            fixed_ligand_flag = f"--fixed_ligand {items}"

        predict_multi_flag = "--predict_multi" if self.predict_multi else ""
        target_res_flag = f"--target_res {self.target_res}" if self.target_res else ""

        ligand_file_flag = ""
        if self.ligand_file is not None:
            parts = []
            if isinstance(self.ligand_file, dict):
                for name, path in self.ligand_file.items():
                    if str(path).upper() == "CCD":
                        parts.append(f"{name}:CCD")
                    else:
                        p = Path(path).resolve()
                        if not p.exists():
                            raise FileNotFoundError(f"Ligand reference file not found: {p}")
                        parts.append(f"{name}:{p}")
            elif isinstance(self.ligand_file, (list, tuple)):
                for item in self.ligand_file:
                    item_str = str(item)
                    if ":" in item_str:
                        name, path = item_str.split(":", 1)
                        if path.upper() == "CCD":
                            parts.append(f"{name}:CCD")
                        else:
                            p = Path(path).resolve()
                            if not p.exists():
                                raise FileNotFoundError(f"Ligand reference file not found: {p}")
                            parts.append(f"{name}:{p}")
                    else:
                        parts.append(item_str)
            else:
                item_str = str(self.ligand_file)
                if ":" in item_str:
                    name, path = item_str.split(":", 1)
                    if path.upper() == "CCD":
                        parts.append(f"{name}:CCD")
                    else:
                        p = Path(path).resolve()
                        if not p.exists():
                            raise FileNotFoundError(f"Ligand reference file not found: {p}")
                        parts.append(f"{name}:{p}")
                else:
                    parts.append(item_str)
            ligand_file_flag = f"--ligand_file {' '.join(parts)}"

        rerank_flag = f"--rerank {self.rerank}" if self.rerank else ""

        bonds_flag = ""
        if self.bonds is not None:
            if isinstance(self.bonds, (list, tuple)):
                items = " ".join(str(x) for x in self.bonds)
            else:
                items = str(self.bonds)
            bonds_flag = f"--bonds {items}"

        mutate_flag = ""
        if self.mutate is not None:
            if isinstance(self.mutate, dict):
                items = " ".join(f"{pos}:{resn}" for pos, resn in self.mutate.items())
            elif isinstance(self.mutate, (list, tuple)):
                items = " ".join(str(x) for x in self.mutate)
            else:
                items = str(self.mutate)
            mutate_flag = f"--mutate {items}"

        crop_centers_flag = ""
        if self.crop_centers is not None:
            if isinstance(self.crop_centers, (list, tuple)):
                items = " ".join(str(x) for x in self.crop_centers)
            else:
                items = str(self.crop_centers)
            crop_centers_flag = f"--crop_centers {items}"

        corruption_centers_flag = ""
        if self.corruption_centers is not None:
            if isinstance(self.corruption_centers, (list, tuple)):
                items = " ".join(str(x) for x in self.corruption_centers)
            else:
                items = str(self.corruption_centers)
            corruption_centers_flag = f"--corruption_centers {items}"

        residue_json_flag = ""
        if self.residue_json is not None:
            res_p = self.residue_json.resolve()
            if not res_p.exists():
                raise FileNotFoundError(f"Residue JSON file not found: {res_p}")
            residue_json_flag = f'--residue_json "{res_p}"'

        fixed_ligand_noise_flag = f"--fixed_ligand_noise {self.fixed_ligand_noise}" if self.fixed_ligand_noise is not None else ""
        exclude_common_ligands_flag = "--exclude_common_ligands" if self.exclude_common_ligands else ""
        ignore_ligand_hydrogens_flag = "--ignore_ligand_hydrogens" if self.ignore_ligand_hydrogens else ""
        use_sm_flag = "" if self.use_sm else "--no-use_sm"
        suffix_flag = f"--suffix {self.suffix}" if self.suffix else ""
        cautious_flag = "--cautious" if self.cautious else ""

        weights_flag = ""
        if self.weights is not None:
            w_p = self.weights.resolve()
            if not w_p.exists():
                raise FileNotFoundError(f"Weights file not found: {w_p}")
            weights_flag = f'--weights "{w_p}"'

        # 3. Device translation
        ribbon_device = self.device if self.device in ["gpu", "gpu_wsl", "cpu"] else "gpu"

        # 4. Invoke Task runner
        self._run_task(
            self.task_name,
            input_flag=input_flag,
            output_dir=str(self.output_dir),
            num_samples=self.num_samples,
            predict_ligand_flag=predict_ligand_flag,
            fixed_ligand_flag=fixed_ligand_flag,
            predict_multi_flag=predict_multi_flag,
            target_res_flag=target_res_flag,
            ligand_file_flag=ligand_file_flag,
            rerank_flag=rerank_flag,
            bonds_flag=bonds_flag,
            mutate_flag=mutate_flag,
            crop_centers_flag=crop_centers_flag,
            corruption_centers_flag=corruption_centers_flag,
            residue_json_flag=residue_json_flag,
            fixed_ligand_noise_flag=fixed_ligand_noise_flag,
            exclude_common_ligands_flag=exclude_common_ligands_flag,
            ignore_ligand_hydrogens_flag=ignore_ligand_hydrogens_flag,
            use_sm_flag=use_sm_flag,
            suffix_flag=suffix_flag,
            cautious_flag=cautious_flag,
            weights_flag=weights_flag,
            extra_args=self.extra_args,
            device=ribbon_device,
        )


Placer = PLACER


class ESMFold2(Task):
    """
    Biohub ESMFold2 and ESMFold2-Fast.

    All-atom co-folding of protein and DNA/RNA/ligand complexes using the 6B-parameter
    ESMC language model and all-atom diffusion trunk.
    """

    def __init__(
        self,
        output_dir: Union[str, Path] = ".",
        fasta_file: Optional[Union[str, Path]] = None,
        sequence: Optional[str] = None,
        sequences: Optional[Union[Dict[str, str], List[Tuple[str, str]]]] = None,
        config: Optional[Union[str, Path]] = None,
        ligands: Optional[Union[str, List[str], Dict[str, str]]] = None,
        apo: bool = True,
        conditions: Optional[List[str]] = None,
        model: str = "full",
        num_diffusion_samples: int = 5,
        num_loops: Optional[int] = None,
        num_sampling_steps: Optional[int] = None,
        chunk_size: Optional[int] = None,
        seeds: Optional[Union[int, List[int]]] = None,
        overwrite: bool = False,
        no_warmup: bool = False,
        device: str = "gpu",
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        """
        Initialize an ESMFold2 task.

        Args:
            output_dir (str or Path): Output directory to save predicted structures and summaries. Default is '.'.
            fasta_file (str or Path, optional): FASTA file with protein sequence(s). Multi-chain complexes use ':' in sequence.
            sequence (str, optional): Single sequence string, optionally prefixed with 'NAME=' (e.g. '1ubq=MQIF...').
            sequences (dict or list of tuple, optional): Sequence collection as {name: seq} or [(name, seq), ...].
            config (str or Path, optional): JSON matrix configuration file ({sequences, ligands, conditions, seeds, settings}).
            ligands (dict, list, or str, optional): Ligand definitions (e.g. {'BTN': 'ccd:BTN', 'TSA': 'smiles:Cc1...'}
                or ['BTN=ccd:BTN', 'TSA=smiles:...']).
            apo (bool): Include apo (ligand-free) folding condition. Default is True.
            conditions (list of str, optional): Explicit subset of defined conditions to run.
            model (str): Model checkpoint ('full' for ESMFold2, 'fast' for ESMFold2-Fast). Default is 'full'.
            num_diffusion_samples (int): Number of diffusion structure samples per condition. Default is 5.
            num_loops (int, optional): Number of recycle/refinement loops. Default is 10.
            num_sampling_steps (int, optional): Number of diffusion sampling steps. Default is 68.
            chunk_size (int, optional): Sequence chunk size for lower peak VRAM on long targets.
            seeds (int or list of int, optional): Random seed(s) (e.g. 5 for seeds 0..4, or [42, 100]).
            overwrite (bool): Force re-running folds already recorded in manifest.jsonl. Default False.
            no_warmup (bool): Skip initial throwaway warmup fold. Default False.
            device (str): Device to run inference on ('gpu', 'gpu_wsl', 'cpu'). Default is 'gpu'.
            extra_args (str): Additional CLI arguments passed to esmfold2_batch.py.
            dry_run (bool): If True, create mock output files without invoking container. Default False.
        """
        # Aliases
        if fasta_file is None and "input_fasta" in kwargs:
            fasta_file = kwargs.pop("input_fasta")
        if fasta_file is None and "input_file" in kwargs:
            inp = kwargs.pop("input_file")
            if str(inp).endswith(".json"):
                config = inp
            else:
                fasta_file = inp
        if "fast" in kwargs:
            if kwargs.pop("fast"):
                model = "fast"
        if "num_samples" in kwargs:
            num_diffusion_samples = kwargs.pop("num_samples")
        if "samples" in kwargs:
            num_diffusion_samples = kwargs.pop("samples")
        if ligands is None:
            if "smiles" in kwargs:
                sm = kwargs.pop("smiles")
                ligands = {"LIG": f"smiles:{sm}" if not str(sm).startswith("smiles:") else str(sm)}
            elif "smiles_string" in kwargs:
                sm = kwargs.pop("smiles_string")
                ligands = {"LIG": f"smiles:{sm}" if not str(sm).startswith("smiles:") else str(sm)}
            elif "ccd" in kwargs:
                c = kwargs.pop("ccd")
                ligands = {str(c): f"ccd:{c}" if not str(c).startswith("ccd:") else str(c)}

        if kwargs:
            raise TypeError(f"ESMFold2.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}")

        input_count = sum(x is not None for x in [fasta_file, sequence, sequences, config])
        if input_count != 1:
            raise ValueError("Provide exactly one of: fasta_file, sequence, sequences, config")

        if model not in ["full", "fast"]:
            raise ValueError(f"Invalid model '{model}'. Must be 'full' (ESMFold2) or 'fast' (ESMFold2-Fast)")

        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "ESMFold2"

        self.output_dir = Path(output_dir)
        self.fasta_file = Path(fasta_file) if fasta_file is not None else None
        self.sequence = sequence
        self.sequences = sequences
        self.config = Path(config) if config is not None else None
        self.ligands = ligands
        self.apo = apo
        self.conditions = conditions
        self.model = model
        self.num_diffusion_samples = num_diffusion_samples
        self.num_loops = num_loops
        self.num_sampling_steps = num_sampling_steps
        self.chunk_size = chunk_size
        self.seeds = seeds
        self.overwrite = overwrite
        self.no_warmup = no_warmup
        self.device = device
        self.dry_run = dry_run

    def _get_target_names(self) -> List[str]:
        names = []
        if self.fasta_file is not None:
            if self.fasta_file.exists():
                try:
                    for line in self.fasta_file.read_text().splitlines():
                        if line.strip().startswith(">"):
                            names.append(line.strip()[1:].split()[0])
                except Exception:
                    pass
            if not names:
                names = [self.fasta_file.stem]
        elif self.sequence is not None:
            if "=" in self.sequence:
                names = [self.sequence.split("=", 1)[0]]
            else:
                names = ["seq1"]
        elif self.sequences is not None:
            if isinstance(self.sequences, dict):
                names = list(self.sequences.keys())
            else:
                names = [n for n, _ in self.sequences]
        else:
            names = ["target"]
        return names

    @property
    def summary_csv(self) -> Path:
        """Path to output summary CSV table."""
        return self.output_dir / "summary.csv"

    @property
    def output_cif_files(self) -> List[Path]:
        """List of generated mmCIF predicted structure files."""
        return sorted(self.output_dir.glob("*.cif"))

    def _run_dry(self):
        """Generate realistic mock outputs for dry-run validation."""
        self.output_dir = make_directory(self.output_dir)
        names = self._get_target_names()

        for name in names:
            cif_path = self.output_dir / f"{name}__apo__s0__m0.cif"
            cif_path.write_text(
                "data_mock\n"
                "_entry.id mock\n"
                "#\n"
                "_atom_site.group_PDB\n"
                "_atom_site.id\n"
                "_atom_site.type_symbol\n"
                "_atom_site.label_atom_id\n"
                "_atom_site.label_comp_id\n"
                "_atom_site.label_asym_id\n"
                "_atom_site.label_seq_id\n"
                "_atom_site.Cartn_x\n"
                "_atom_site.Cartn_y\n"
                "_atom_site.Cartn_z\n"
                "_atom_site.occupancy\n"
                "_atom_site.B_iso_or_equiv\n"
                "ATOM 1 N N SER A 1 0.000 0.000 0.000 1.00 85.0\n"
            )

        summary_path = self.output_dir / "summary.csv"
        with open(summary_path, "w") as f:
            f.write("ok,name,cond,seed,n_res,n_tokens,checkpoint,seed_samples,plddt_mean,plddt_best,ptm_mean,iptm_mean,fold_s,peak_vram_gb,error,cif\n")
            for name in names:
                f.write(f"True,{name},apo,0,76,76,ESMFold2,1,0.85,0.85,0.80,,3.5,14.0,,{name}__apo__s0__m0.cif\n")

    def run(self):
        """Execute ESMFold2 inference or simulate dry run."""
        if self.dry_run:
            self._run_dry()
            return

        self.output_dir = make_directory(self.output_dir.resolve())

        # Ensure user-owned cache dirs exist
        for d in ["/tmp/xdg_cache_ribbon", "/tmp/triton_cache_ribbon", "/tmp/inductor_cache_ribbon", "/tmp/hf_cache_ribbon"]:
            Path(d).mkdir(parents=True, exist_ok=True)

        # 1. Format input flags
        if self.fasta_file is not None:
            inp = self.fasta_file.resolve()
            if not inp.exists():
                raise FileNotFoundError(f"FASTA file not found: {inp}")
            input_flags = f'--fasta "{inp}"'
        elif self.sequence is not None:
            if "=" in self.sequence:
                seq_str = self.sequence
            else:
                seq_str = f"seq1={self.sequence}"
            input_flags = f'--sequence "{seq_str}"'
        elif self.sequences is not None:
            flags = []
            if isinstance(self.sequences, dict):
                for n, q in self.sequences.items():
                    flags.append(f'--sequence "{n}={q}"')
            else:
                for n, q in self.sequences:
                    flags.append(f'--sequence "{n}={q}"')
            input_flags = " ".join(flags)
        elif self.config is not None:
            cfg = self.config.resolve()
            if not cfg.exists():
                raise FileNotFoundError(f"Config file not found: {cfg}")
            input_flags = f'--config "{cfg}"'

        # 2. Format ligand flags
        lflags = []
        if self.apo:
            lflags.append("--apo")
        if self.ligands is not None:
            if isinstance(self.ligands, dict):
                for name, spec in self.ligands.items():
                    lflags.append(f'--ligand "{name}={spec}"')
            elif isinstance(self.ligands, (list, tuple)):
                for item in self.ligands:
                    lflags.append(f'--ligand "{item}"')
            else:
                lflags.append(f'--ligand "{self.ligands}"')
        if not self.apo and not self.ligands:
            lflags.append("--apo")
        ligand_flags = " ".join(lflags)

        # 3. Model & execution flags
        model_flag = f"--model {self.model}"
        num_diffusion_samples_flag = f"--num-diffusion-samples {self.num_diffusion_samples}"
        num_loops_flag = f"--num-loops {self.num_loops}" if self.num_loops is not None else ""
        num_sampling_steps_flag = f"--num-sampling-steps {self.num_sampling_steps}" if self.num_sampling_steps is not None else ""
        chunk_size_flag = f"--chunk-size {self.chunk_size}" if self.chunk_size is not None else ""

        if self.seeds is not None:
            if isinstance(self.seeds, int):
                seeds_flag = f"--seeds {self.seeds}"
            else:
                seeds_flag = f"--seed-list {','.join(str(s) for s in self.seeds)}"
        else:
            seeds_flag = ""

        overwrite_flag = "--overwrite" if self.overwrite else ""
        no_warmup_flag = "--no-warmup" if self.no_warmup else ""

        # 4. Device translation
        ribbon_device = self.device if self.device in ["gpu", "gpu_wsl", "cpu"] else "gpu"

        # 5. Invoke runner
        self._run_task(
            self.task_name,
            input_flags=input_flags,
            output_dir=str(self.output_dir),
            ligand_flags=ligand_flags,
            model_flag=model_flag,
            num_diffusion_samples_flag=num_diffusion_samples_flag,
            num_loops_flag=num_loops_flag,
            num_sampling_steps_flag=num_sampling_steps_flag,
            chunk_size_flag=chunk_size_flag,
            seeds_flag=seeds_flag,
            overwrite_flag=overwrite_flag,
            no_warmup_flag=no_warmup_flag,
            extra_args=self.extra_args,
            device=ribbon_device,
        )


Esmfold2 = ESMFold2