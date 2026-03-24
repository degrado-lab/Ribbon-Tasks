from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Sequence, Tuple, Union
import os

import argparse
from tempfile import NamedTemporaryFile
import subprocess

# --------------------------
# Minimal YAML emitter (no external deps)
# --------------------------

def _quote_scalar(value: str) -> str:
    # Always single-quote strings; escape single quotes by doubling
    return "'" + value.replace("'", "''") + "'"


def _is_scalar(x: Any) -> bool:
    return isinstance(x, (str, int, float)) or isinstance(x, bool) or x is None


def dump_yaml(data: Any, indent: int = 0) -> str:
    lines: List[str] = []
    ind = " " * indent

    if isinstance(data, dict):
        for i, (k, v) in enumerate(data.items()):
            key = str(k)
            if _is_scalar(v) or (isinstance(v, list) and len(v) == 0):
                lines.append(f"{ind}{key}: {format_scalar(v)}")
            elif isinstance(v, list):
                lines.append(f"{ind}{key}:")
                for item in v:
                    if _is_scalar(item):
                        lines.append(f"{ind}  - {format_scalar(item)}")
                    elif isinstance(item, dict):
                        # emit '- key: value' on one line if single-key dict with scalar value
                        if len(item) == 1:
                            (sk, sv), = item.items()
                            if _is_scalar(sv):
                                lines.append(f"{ind}  - {sk}: {format_scalar(sv)}")
                                continue
                        lines.append(f"{ind}  - {dump_yaml(item, indent + 4).lstrip()}")
                    else:
                        lines.append(f"{ind}  - {dump_yaml(item, indent + 4).lstrip()}")
            else:
                lines.append(f"{ind}{key}:")
                lines.append(dump_yaml(v, indent + 2))
    elif isinstance(data, list):
        for item in data:
            if _is_scalar(item):
                lines.append(f"{ind}- {format_scalar(item)}")
            else:
                lines.append(f"{ind}-")
                lines.append(dump_yaml(item, indent + 2))
    else:
        lines.append(f"{ind}{format_scalar(data)}")

    return "\n".join(lines)


def format_scalar(v: Any) -> str:
    if v is None:
        return "null"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return str(v)
    # strings
    return _quote_scalar(str(v))

# --------------------------
# Write Boltz2 YAML
# --------------------------

def build_boltz_yaml(
    sequences: List[Dict[str, Dict[str, Any]]],
    *,
    constraints: Optional[List[Dict[str, Any]]] = None,
    templates: Optional[List[Dict[str, Any]]] = None,
    properties: Optional[List[Dict[str, Any]]] = None,
) -> str:
    """Assemble, validate, and serialize Boltz input YAML.

    Returns a YAML string.
    """
    # Build top-level mapping in desired order
    data: Dict[str, Any] = {"sequences": sequences}
    if constraints:
        data["constraints"] = constraints
    if templates:
        data["templates"] = templates
    if properties:
        data["properties"] = properties

    #validate_boltz_input(data)
    return dump_yaml(data)


def write_boltz_yaml(
    path: str,
    sequences: List[Dict[str, Dict[str, Any]]],
    *,
    constraints: Optional[List[Dict[str, Any]]] = None,
    templates: Optional[List[Dict[str, Any]]] = None,
    properties: Optional[List[Dict[str, Any]]] = None,
) -> str:
    yaml_text = build_boltz_yaml(
        sequences,
        constraints=constraints,
        templates=templates,
        properties=properties,
    )
    with open(path, "w", encoding="utf-8") as f:
        f.write(yaml_text)
    return os.path.abspath(path)

# --------------------------
# Data
# --------------------------

def protein(
    id: Union[str, List[str]],
    sequence: str,
    msa: Optional[str] = None,
    modifications: Optional[List[Dict[str, Union[int, str]]]] = None,
    cyclic: Optional[bool] = None,
    use_msa_server: bool = False,
) -> Dict[str, Dict[str, Any]]:
    entry: Dict[str, Any] = {"id": id, "sequence": sequence}
    if not use_msa_server:
        if msa is not None:
            entry["msa"] = msa
        else:
            entry["msa"] = "empty"
    if modifications:
        entry["modifications"] = modifications
    if cyclic is not None:
        entry["cyclic"] = bool(cyclic)
    return {"protein": entry}


def ligand_smiles(
    id: Union[str, List[str]],
    smiles: str,
) -> Dict[str, Dict[str, Any]]:
    return {"ligand": {"id": id, "smiles": smiles}}

def affinity(binder_chain_id: str) -> Dict[str, Any]:
    return {"affinity": {"binder": binder_chain_id}}

# ------------------------
# CLI
# ------------------------


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add Hydrogens to a PDB file. Create a custom DB if required."
    )
    parser.add_argument("input_fasta", type=str,
                        help="Path to the input FASTA file.")
    parser.add_argument("output_dir", type=str,
                        help="Path to the output PDB file.")
    parser.add_argument("--ligand_smiles", type=str, nargs="*",
                        help="List of custom ligand SMILES strings.")
    parser.add_argument("--calculate_binding", action="store_true",
                        help="Whether to calculate binding.")
    parser.add_argument("--use_msa_server", action="store_true",
                        help="Whether to use MSA for structure prediction.")
    parser.add_argument("--extra_args", type=str,
                        help="Additional arguments for the hydrogen addition.")

    args = parser.parse_args()
    if args.ligand_smiles is None:
        args.ligand_smiles = []

    # Read the fasta file and store the sequences (names are dropped)
    sequence_list = []
    with open(args.input_fasta, "r") as fasta_file:
        current_sequence = ""
        for line in fasta_file:
            if line.startswith(">"):
                if current_sequence:
                    sequence_list.append(current_sequence)
                    current_sequence = ""
            else:
                current_sequence += line.strip()
        # Add the last sequence if file doesn't end with a header
        if current_sequence:
            sequence_list.append(current_sequence)

    # Now we have a list of sequences (sequence_list) and a list of ligands (args.ligand_smiles).
    # Lets assign chains to each of these (A, B, C...)
    protein_chains = [f"{chr(65 + i)}" for i in range(len(sequence_list))]
    ligand_chains = [f"{chr(65 + len(sequence_list) + i)}" for i in range(len(args.ligand_smiles))]
    # If we have more than 26 it will crap out lol

    # Use this to create the protein objects!
    protein_objects = [
        protein(id=chain, sequence=seq, use_msa_server=args.use_msa_server)
        for chain, seq in zip(protein_chains, sequence_list)
    ]
    # And, create the ligand objects!
    ligand_objects = [
        ligand_smiles(id=chain, smiles=smiles)
        for chain, smiles in zip(ligand_chains, args.ligand_smiles)
    ]
    # If we're calculating properties, make an affinity object for the first ligand!
    if args.calculate_binding:
        affinity_objects = [ affinity(binder_chain_id=ligand_chains[0]) ]
    else:
        affinity_objects = []

    # Write to a temporary YAML:
    input_yaml = NamedTemporaryFile(delete=False, suffix=".yaml")
    write_boltz_yaml(path=input_yaml.name, sequences=protein_objects + ligand_objects, properties=affinity_objects)

    extra_args = "" if args.extra_args is None else args.extra_args
    if args.use_msa_server:
        extra_args += " --use_msa_server "
    # Run the command in the shell:
    command = f"boltz predict {input_yaml.name} --out_dir {args.output_dir} --no_kernels {extra_args}"
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)