from pathlib import Path
from typing import Any, Callable


CaseBuilder = Callable[[Path, Path], dict[str, Any]]

CASES_DIR = Path(__file__).parent / "cases"


def _reduce_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "pdb_input_file": str(case_dir / "inputs" / "minimal.pdb"),
        "pdb_output_file": str(out_dir / "reduced.pdb"),
        "flip": False,
    }


def _distance_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "pdb_file": str(case_dir / "inputs" / "minimal.pdb"),
        "atom1": "A:1:CA",
        "atom2": "A:2:CA",
        "output_file": str(out_dir / "distance.txt"),
        "device": "cpu",
    }


def _add_h_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "input_file": str(case_dir / "inputs" / "minimal.pdb"),
        "output_file": str(out_dir / "with_h.pdb"),
        "selection": "all",
    }


def _rosetta_ligand_prepare_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "input_file": str(case_dir / "inputs" / "minimal.sdf"),
        "output_dir": str(out_dir),
        "name": "LG1",
        "clobber": True,
    }


def _angle_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "pdb_file": str(case_dir / "inputs" / "minimal.pdb"),
        "atom1": "A:1:N",
        "atom2": "A:1:CA",
        "atom3": "A:1:C",
        "output_file": str(out_dir / "angle.txt"),
        "device": "cpu",
    }


def _dihedral_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "pdb_file": str(case_dir / "inputs" / "minimal.pdb"),
        "atom1": "A:1:N",
        "atom2": "A:1:CA",
        "atom3": "A:1:C",
        "atom4": "A:2:N",
        "output_file": str(out_dir / "dihedral.txt"),
        "device": "cpu",
    }


def _pairwise_distance_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "pdb_file": str(case_dir / "inputs" / "minimal.pdb"),
        "atom_list_A": ["A:1:CA"],
        "atom_list_B": ["A:2:CA"],
        "output_file": str(out_dir / "pairwise.json"),
        "average": False,
        "device": "cpu",
    }


def _sasa_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "pdb_file": str(case_dir / "inputs" / "minimal.pdb"),
        "output_file": str(out_dir / "sasa.txt"),
        "atom_1": "A:1:CA",
        "device": "cpu",
    }


def _fast_relax_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "pdb_input_file": str(case_dir / "inputs" / "minimal.pdb"),
        "output_dir": str(out_dir),
        "nstructs": 1,
        "device": "cpu",
    }


def _ligand_mpnn_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "output_dir": str(out_dir),
        "structure_list": [str(case_dir / "inputs" / "minimal.pdb")],
        "num_designs": 1,
        "temperature": 0.1,
        "device": "gpu",
    }


def _laser_mpnn_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "output_dir": str(out_dir),
        "structure_list": [str(case_dir / "inputs" / "minimal.pdb")],
        "num_designs": 1,
        "temperature": 0.000001,
        "device": "gpu",
    }


def _chai1_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "fasta_file": str(case_dir / "inputs" / "minimal.fasta"),
        "output_dir": str(out_dir),
        "device": "gpu",
    }


def _boltz2_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "fasta_file": str(case_dir / "inputs" / "minimal.fasta"),
        "output_dir": str(out_dir),
        "device": "gpu",
    }


def _raptorx_single_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "fasta_file_or_dir": str(case_dir / "inputs" / "minimal.fasta"),
        "output_dir": str(out_dir),
        "device": "gpu",
    }


def _rfdiffusion_aa_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "input_structure": str(case_dir / "inputs" / "minimal.pdb"),
        "output_dir": str(out_dir),
        "contig_map": "[5-5]",
        "num_designs": 1,
        "device": "gpu",
    }


def _easy_md_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "input_file": str(case_dir / "inputs" / "minimal.pdb"),
        "output_prefix": str(out_dir / "traj"),
        "duration": 1,
        "device": "gpu",
    }


def _custom_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
    return {
        "command": f"echo 'done' > {out_dir}/custom.txt",
        "container": "BioPython",
        "device": "cpu",
    }


SMOKE_CASES: list[dict[str, Any]] = [
    {
        "name": "Reduce",
        "class_name": "Reduce",
        "case_dir": CASES_DIR / "reduce",
        "kwargs_builder": _reduce_case,
        "validations": [
            {"kind": "exists", "path": "reduced.pdb"},
            {"kind": "nonempty", "path": "reduced.pdb"},
        ],
        "timeout_s": 180,
        "requires_gpu": False,
    },
    {
        "name": "Calculate Distance",
        "class_name": "CalculateDistance",
        "case_dir": CASES_DIR / "calculate_distance",
        "kwargs_builder": _distance_case,
        "validations": [
            {"kind": "exists", "path": "distance.txt"},
            {"kind": "nonempty", "path": "distance.txt"},
        ],
        "timeout_s": 120,
        "requires_gpu": False,
    },
    {
        "name": "Calculate Angle",
        "class_name": "CalculateAngle",
        "case_dir": CASES_DIR / "calculate_angle",
        "kwargs_builder": _angle_case,
        "validations": [
            {"kind": "exists", "path": "angle.txt"},
            {"kind": "nonempty", "path": "angle.txt"},
        ],
        "timeout_s": 120,
        "requires_gpu": False,
    },
    {
        "name": "Calculate Dihedral",
        "class_name": "CalculateDihedral",
        "case_dir": CASES_DIR / "calculate_dihedral",
        "kwargs_builder": _dihedral_case,
        "validations": [
            {"kind": "exists", "path": "dihedral.txt"},
            {"kind": "nonempty", "path": "dihedral.txt"},
        ],
        "timeout_s": 120,
        "requires_gpu": False,
    },
    {
        "name": "Calculate Pairwise Distance",
        "class_name": "CalculatePairwiseDistance",
        "case_dir": CASES_DIR / "calculate_pairwise_distance",
        "kwargs_builder": _pairwise_distance_case,
        "validations": [
            {"kind": "exists", "path": "pairwise.json"},
            {"kind": "nonempty", "path": "pairwise.json"},
        ],
        "timeout_s": 120,
        "requires_gpu": False,
    },
    {
        "name": "Calculate SASA",
        "class_name": "CalculateSASA",
        "case_dir": CASES_DIR / "calculate_sasa",
        "kwargs_builder": _sasa_case,
        "validations": [
            {"kind": "exists", "path": "sasa.txt"},
        ],
        "timeout_s": 120,
        "requires_gpu": False,
        "expect_not_implemented": True,
    },
    {
        "name": "Add Hydrogens",
        "class_name": "AddHydrogens",
        "case_dir": CASES_DIR / "add_hydrogens",
        "kwargs_builder": _add_h_case,
        "validations": [
            {"kind": "exists", "path": "with_h.pdb"},
            {"kind": "nonempty", "path": "with_h.pdb"},
            {"kind": "contains_text", "path": "with_h.pdb", "text": "ATOM"},
        ],
        "timeout_s": 180,
        "requires_gpu": False,
    },
    {
        "name": "RosettaLigandPrepare",
        "class_name": "RosettaLigandPrepare",
        "case_dir": CASES_DIR / "rosetta_ligand_prepare",
        "kwargs_builder": _rosetta_ligand_prepare_case,
        "validations": [
            {"kind": "exists", "path": "LG1.params"},
            {"kind": "nonempty", "path": "LG1.params"},
            {"kind": "exists", "path": "LG1_0001.pdb"},
            {"kind": "nonempty", "path": "LG1_0001.pdb"},
        ],
        "timeout_s": 180,
        "requires_gpu": False,
    },
    {
        "name": "FastRelax",
        "class_name": "FastRelax",
        "case_dir": CASES_DIR / "fast_relax",
        "kwargs_builder": _fast_relax_case,
        "validations": [
            {"kind": "exists", "path": "score.sc"},
        ],
        "timeout_s": 180,
        "requires_gpu": False,
    },
    {
        "name": "Custom",
        "class_name": "Custom",
        "case_dir": CASES_DIR / "custom",
        "kwargs_builder": _custom_case,
        "validations": [
            {"kind": "exists", "path": "custom.txt"},
            {"kind": "nonempty", "path": "custom.txt"},
        ],
        "timeout_s": 120,
        "requires_gpu": False,
    },
    {
        "name": "LigandMPNN",
        "class_name": "LigandMPNN",
        "case_dir": CASES_DIR / "ligand_mpnn",
        "kwargs_builder": _ligand_mpnn_case,
        "validations": [
            {"kind": "exists", "path": "seqs_split"},
        ],
        "timeout_s": 300,
        "requires_gpu": True,
    },
    {
        "name": "LASErMPNN",
        "class_name": "LASErMPNN",
        "case_dir": CASES_DIR / "laser_mpnn",
        "kwargs_builder": _laser_mpnn_case,
        "validations": [
            {"kind": "exists", "path": "minimal/design_0.pdb"},
        ],
        "timeout_s": 300,
        "requires_gpu": True,
    },
    {
        "name": "Chai-1",
        "class_name": "Chai1",
        "case_dir": CASES_DIR / "chai1",
        "kwargs_builder": _chai1_case,
        "validations": [
            {"kind": "exists", "path": "minimal_idx_0.cif"},
        ],
        "timeout_s": 600,
        "requires_gpu": True,
    },
    {
        "name": "Boltz-2",
        "class_name": "Boltz2",
        "case_dir": CASES_DIR / "boltz2",
        "kwargs_builder": _boltz2_case,
        "validations": [
            {"kind": "glob_exists", "pattern": "boltz_results_*/predictions/*/*_model_0.cif"},
        ],
        "timeout_s": 600,
        "requires_gpu": True,
    },
    {
        "name": "RaptorXSingle",
        "class_name": "RaptorXSingle",
        "case_dir": CASES_DIR / "raptorx_single",
        "kwargs_builder": _raptorx_single_case,
        "validations": [],
        "timeout_s": 300,
        "requires_gpu": True,
    },
    {
        "name": "RFDiffusionAA",
        "class_name": "RFDiffusionAA",
        "case_dir": CASES_DIR / "rfdiffusion_aa",
        "kwargs_builder": _rfdiffusion_aa_case,
        "validations": [],
        "timeout_s": 600,
        "requires_gpu": True,
    },
    {
        "name": "EasyMD",
        "class_name": "EasyMD",
        "case_dir": CASES_DIR / "easy_md",
        "kwargs_builder": _easy_md_case,
        "validations": [],
        "timeout_s": 600,
        "requires_gpu": True,
    },
]
