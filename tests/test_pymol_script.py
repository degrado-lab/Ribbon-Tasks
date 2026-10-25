from pathlib import Path
from unittest.mock import MagicMock
import pytest
import ribbon
from ribbon_tasks.tasks import PyMolScript, PyMOLScript, RunPyMolScript

TESTS_DIR = Path(__file__).parent
PYMOL_SCRIPT_INPUTS = TESTS_DIR / "smoke" / "cases" / "pymol_script" / "inputs"


def test_pymol_script_init_defaults(tmp_path):
    script = tmp_path / "test.py"
    script.write_text("print('hello')\n")

    task = PyMolScript(script_file=script)
    assert task.task_name == "PyMolScript"
    assert task.script_file == script
    assert task.script_args is None
    assert task.output_dir is None
    assert task.output_file is None
    assert task.device == "cpu"
    assert task.extra_args == ""
    assert task.dry_run is False


def test_pymol_script_aliases(tmp_path):
    script = tmp_path / "test.py"
    script.write_text("print('hello')\n")
    out_dir = tmp_path / "out"
    out_file = out_dir / "result.txt"

    # Alias: script, args, out_dir, out_file, device_type
    task1 = PyMolScript(
        script=script,
        args=["--flag", "val"],
        out_dir=out_dir,
        out_file=out_file,
        device_type="cpu",
    )
    assert task1.script_file == script
    assert task1.script_args == ["--flag", "val"]
    assert task1.output_dir == out_dir
    assert task1.output_file == out_file
    assert task1.device == "cpu"

    # Alias: python_script, arguments, output_folder
    task2 = PyMolScript(
        python_script=script,
        arguments="--flag val",
        output_folder=out_dir,
    )
    assert task2.script_file == script
    assert task2.script_args == "--flag val"
    assert task2.output_dir == out_dir

    # Alias: script_path, params
    task3 = PyMolScript(
        script_path=script,
        params="--quiet",
    )
    assert task3.script_file == script
    assert task3.script_args == "--quiet"

    # Alias: input_file / file
    task4 = PyMolScript(input_file=script)
    assert task4.script_file == script

    task5 = PyMolScript(file=script)
    assert task5.script_file == script


def test_pymol_script_validation_errors(tmp_path):
    # Missing script_file
    with pytest.raises(ValueError, match="Argument 'script_file' is required"):
        PyMolScript()

    script = tmp_path / "test.py"
    script.write_text("print('hello')\n")

    # Invalid device
    with pytest.raises(ValueError, match="Invalid device"):
        PyMolScript(script_file=script, device="tpu")

    # Unexpected keyword argument
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        PyMolScript(script_file=script, unknown_option="foo")


def test_pymol_script_file_not_found(tmp_path):
    missing_script = tmp_path / "does_not_exist.py"
    task = PyMolScript(script_file=missing_script)
    with pytest.raises(FileNotFoundError, match="Script file not found"):
        task.run()


def test_pymol_script_flag_compilation(tmp_path):
    script = tmp_path / "test.py"
    script.write_text("print('hello')\n")
    out_file = tmp_path / "out" / "out.txt"

    # List script_args
    task = PyMolScript(
        script_file=script,
        script_args=["--foo", "bar", "--num", "10"],
        output_file=out_file,
        extra_args="--verbose",
        device="cpu",
    )
    task._run_task = MagicMock()
    task.run()

    task._run_task.assert_called_once_with(
        "PyMolScript",
        script_file=str(script.resolve()),
        script_args="--foo bar --num 10",
        extra_args="--verbose",
        device="cpu",
    )
    assert out_file.parent.exists()


def test_pymol_script_dry_run(tmp_path):
    script = tmp_path / "test.py"
    out_dir = tmp_path / "out_dry"
    out_file = tmp_path / "out_file_dry" / "result.txt"

    task = PyMolScript(
        script_file=script,
        output_dir=out_dir,
        output_file=out_file,
        dry_run=True,
    )
    task.run()

    assert (out_dir / "output.txt").exists()
    assert out_file.exists()
    assert "Mock output from dry run of PyMolScript" in (out_dir / "output.txt").read_text()
    assert "Mock output from dry run of PyMolScript" in out_file.read_text()


def test_pymol_script_root_module_import():
    assert hasattr(ribbon, "PyMolScript")
    assert ribbon.PyMolScript is PyMolScript
    assert hasattr(ribbon, "PyMOLScript")
    assert ribbon.PyMOLScript is PyMolScript
    assert hasattr(ribbon, "RunPyMolScript")
    assert ribbon.RunPyMolScript is PyMolScript


def test_pymol_script_container_smoke(tmp_path):
    container_sif = Path.home() / ".ribbon" / "ribbon_containers" / "pymol_2.3.0.sif"
    if not container_sif.exists():
        pytest.skip("pymol_2.3.0.sif container image not found.")

    script = PYMOL_SCRIPT_INPUTS / "run_pymol.py"
    assert script.exists(), f"Smoke fixture script not found at {script}"

    out_file = tmp_path / "out" / "pymol_result.txt"
    task = PyMolScript(
        script_file=script,
        script_args=[str(out_file)],
        output_file=out_file,
        device="cpu",
    )
    task.run()

    assert out_file.exists(), "Output file missing after PyMolScript run"
    assert "PyMol script executed successfully!" in out_file.read_text()
