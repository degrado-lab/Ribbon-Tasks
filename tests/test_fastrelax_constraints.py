import pytest
from unittest.mock import MagicMock, patch
from pathlib import Path
import tempfile
import shutil

from ribbon_tasks.tasks import FastRelax

def test_parse_atom_selector():
    # We can create a dummy task instance to test the method
    # FastRelax requires output_dir and either pdb_input_file or pdb_input_dir in __init__
    task = FastRelax(output_dir=".", pdb_input_file="dummy.pdb")
    
    # Standard Chain:Res:Atom
    assert task._parse_atom_selector("A:16:ND1") == "ND1 16A"
    assert task._parse_atom_selector("B:1:H1") == "H1 1B"
    
    # Continuous / no chain numbering
    assert task._parse_atom_selector(":45:OD2") == "OD2 45"
    
    # Invalid formats should raise ValueError
    with pytest.raises(ValueError, match="Invalid atom selector"):
        task._parse_atom_selector("A:16")
    with pytest.raises(ValueError, match="Invalid atom selector"):
        task._parse_atom_selector("A:16:ND1:extra")


def test_generate_constraints_file(tmp_path):
    # Test generation of constraints file
    task = FastRelax(
        output_dir=tmp_path,
        pdb_input_file="dummy.pdb",
        custom_bonds=[
            "A:16:ND1,B:1:H1,2.0,0.5",  # Default HARMONIC
            "A:45:OD2,A:16:HE2,FLAT_HARMONIC,0.0,0.5,6"  # Custom FLAT_HARMONIC with 3 parameters
        ],
        custom_angles=[
            "A:16:CE1,A:16:ND1,B:1:H1,2.41,0.5",  # Default HARMONIC
            "A:1:CA,A:2:CA,A:3:CA,BOUNDED,1.0,2.0,0.1,TAG"  # Custom BOUNDED with 4 parameters
        ],
        custom_torsions=[
            "A:45:OD1,A:45:CG,A:45:OD2,A:16:HE2,3.14,0.35"  # Default CIRCULARHARMONIC
        ]
    )
    
    cst_file = task._generate_constraints_file()
    assert cst_file is not None
    assert cst_file.exists()
    assert cst_file == tmp_path / "constraints.cst"
    
    content = cst_file.read_text().splitlines()
    assert len(content) == 5
    assert content[0] == "AtomPair ND1 16A H1 1B HARMONIC 2.0 0.5"
    assert content[1] == "AtomPair OD2 45A HE2 16A FLAT_HARMONIC 0.0 0.5 6"
    assert content[2] == "Angle CE1 16A ND1 16A H1 1B HARMONIC 2.41 0.5"
    assert content[3] == "Angle CA 1A CA 2A CA 3A BOUNDED 1.0 2.0 0.1 TAG"
    assert content[4] == "Dihedral OD1 45A CG 45A OD2 45A HE2 16A CIRCULARHARMONIC 3.14 0.35"


def test_generate_constraints_file_empty(tmp_path):
    task = FastRelax(output_dir=tmp_path, pdb_input_file="dummy.pdb")
    cst_file = task._generate_constraints_file()
    assert cst_file is None


@patch("ribbon_tasks.tasks.list_files")
@patch("ribbon_tasks.tasks.shutil.copy")
@patch("ribbon_tasks.tasks.tempfile.mkdtemp")
def test_run_compiles_arguments(mock_mkdtemp, mock_copy, mock_list_files, tmp_path):
    mock_mkdtemp.return_value = str(tmp_path / "temp_pdb_dir")
    mock_list_files.return_value = [Path("temp_pdb_dir/input.pdb")]
    
    # Create a dummy ligand params file so we can resolve its path successfully
    dummy_params = tmp_path / "LIG.params"
    dummy_params.write_text("# dummy params")
    
    task = FastRelax(
        output_dir=tmp_path / "out",
        pdb_input_file="dummy.pdb",
        nstructs=2,
        ligand_params_file=dummy_params,
        custom_bonds=["A:16:ND1,B:1:H1,2.0,0.5"],
        constraints_weight=10.0,
        extra_args="-relax:fast"
    )
    
    # Mock the _run_task method inherited from Task
    task._run_task = MagicMock()
    
    task.run()
    
    # Verify that constraints.cst was generated
    expected_cst_path = (tmp_path / "out" / "constraints.cst").resolve()
    assert expected_cst_path.exists()
    
    # Retrieve the args passed to _run_task
    task._run_task.assert_called_once()
    called_kwargs = task._run_task.call_args[1]
    
    assert called_kwargs["pdb_string"] == "temp_pdb_dir/input.pdb "
    assert called_kwargs["output_dir"] == str((tmp_path / "out"))
    assert called_kwargs["nstructs"] == 2
    
    # extra_args should contain constraints file, weight, ligand params, and user extra_args
    extra_args_val = called_kwargs["extra_args"]
    assert f"-constraints:cst_fa_file {expected_cst_path}" in extra_args_val
    assert "-constraints:cst_fa_weight 10.0" in extra_args_val
    assert f"-extra_res_fa {dummy_params.resolve()}" in extra_args_val
    assert "-relax:fast" in extra_args_val


def test_ligand_params_file_not_found(tmp_path):
    dummy_pdb = tmp_path / "dummy.pdb"
    dummy_pdb.touch()
    
    task = FastRelax(
        output_dir=tmp_path,
        pdb_input_file=dummy_pdb,
        ligand_params_file="nonexistent.params"
    )
    
    with pytest.raises(FileNotFoundError, match="Ligand params file not found"):
        task.run()
