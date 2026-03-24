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
            {
                "kind": "numeric_tolerance",
                "path": "distance.txt",
                "expected_path": "distance.txt",
                "atol": 0.05,
            },
        ],
        "timeout_s": 120,
        "requires_gpu": False,
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
]

# TODO:
# 1. Add entries for all remaining wrappers in ribbon_tasks/tasks.py.
# 2. Split heavyweight tasks into a separate release lane if needed.
# 3. Add richer validators for additional structured outputs (CSV/JSON/PDB semantics).
