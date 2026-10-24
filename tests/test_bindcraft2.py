import json
from pathlib import Path
from unittest.mock import MagicMock

import pytest

import ribbon
from ribbon_tasks.tasks import BindCraft2

TESTS_DIR = Path(__file__).parent
BINDCRAFT2_INPUTS = TESTS_DIR / "smoke" / "cases" / "bindcraft2" / "inputs"


def test_bindcraft2_init_defaults(tmp_path):
    out_dir = tmp_path / "out"
    pdb_file = tmp_path / "target.pdb"
    pdb_file.write_text("ATOM 1\n")

    task = BindCraft2(
        output_dir=out_dir,
        target_pdb=pdb_file,
    )
    assert task.task_name == "BindCraft2"
    assert task.modality == "binder"
    assert task.number_of_final_designs == 1
    assert task.device == "gpu"
    assert task.dry_run is False
    assert task.extra_args == ""
    assert task.core is None
    assert task.bigbang is False
    assert task.disulfide_staple is False
    assert task.forced_targeting is False
    assert task.humanize is False
    assert task.initial_guess is False
    assert task.mixed_topology is False
    assert task.protease_stable is False
    assert task.termini_accessible is False
    assert task.termini_together is False


def test_bindcraft2_aliases(tmp_path):
    out_dir = tmp_path / "out"
    pdb_file = tmp_path / "target.pdb"
    pdb_file.write_text("ATOM 1\n")

    task = BindCraft2(
        output_dir=out_dir,
        pdb_file=pdb_file,
        target_chains=["A", "B"],
        target_hotspots=[54, 56, 66],
        binder_length=70,
        num_designs=5,
        max_attempts=50,
        seed=123,
        set_settings={"autotune": True},
    )
    assert task.target_pdb == pdb_file
    assert task.chains == "A,B"
    assert task.hotspots == "54,56,66"
    assert task.binder_lengths == [70, 70]
    assert task.number_of_final_designs == 5
    assert task.max_trajectories == 50
    assert task.campaign_seed == 123
    assert task.settings_overrides == {"autotune": True}


def test_bindcraft2_validation(tmp_path):
    out_dir = tmp_path / "out"
    pdb_file = tmp_path / "target.pdb"
    pdb_file.write_text("ATOM 1\n")

    # Missing target
    with pytest.raises(ValueError, match="Must specify either"):
        BindCraft2(output_dir=out_dir)

    # Invalid modality
    with pytest.raises(ValueError, match="Invalid modality"):
        BindCraft2(output_dir=out_dir, target_pdb=pdb_file, modality="invalid_modality")

    # Unexpected keyword argument
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        BindCraft2(output_dir=out_dir, target_pdb=pdb_file, unknown_argument=999)


def test_bindcraft2_run_generates_settings_file(tmp_path):
    out_dir = tmp_path / "out"
    pdb_file = tmp_path / "target.pdb"
    pdb_file.write_text("ATOM 1\n")

    task = BindCraft2(
        output_dir=out_dir,
        target_pdb=pdb_file,
        target_name="test_target",
        chains=["A"],
        hotspots=[54, 56, 115],
        modality="VHH",
        binder_lengths=[60, 90],
        number_of_final_designs=3,
        max_trajectories=100,
        campaign_seed=42,
        core="benchmark",
        humanize=True,
        forced_targeting=True,
        screen_steps=2,
        refine_steps=2,
        design_recycles=1,
        settings_overrides={"autotune": False},
        device="gpu",
        extra_args="--test_extra",
    )

    task._run_task = MagicMock()
    task.run()

    settings_file = out_dir / "campaign_settings.json"
    assert settings_file.exists()
    cfg = json.loads(settings_file.read_text())

    assert cfg["project_folder"] == str(out_dir.resolve())
    assert cfg["modality"] == "VHH"
    assert cfg["number_of_final_designs"] == 3
    assert cfg["max_trajectories"] == 100
    assert cfg["campaign_seed"] == 42
    assert cfg["screen_steps"] == 2
    assert cfg["refine_steps"] == 2
    assert cfg["design_recycles"] == 1
    assert cfg["autotune"] is False
    assert len(cfg["targets"]) == 1
    assert cfg["targets"][0]["name"] == "test_target"
    assert cfg["targets"][0]["target_path"] == str(pdb_file.resolve())
    assert cfg["targets"][0]["chains"] == "A"
    assert cfg["targets"][0]["hotspots"] == "54,56,115"

    task._run_task.assert_called_once_with(
        "BindCraft2",
        settings_file=str(settings_file.resolve()),
        core_flag="--core benchmark",
        modality_flag="--modality VHH",
        property_flags="--forced-targeting --humanize",
        set_flags="",
        metadata_flag="",
        extra_args="--test_extra",
        device="gpu",
    )


def test_bindcraft2_run_with_existing_settings_file(tmp_path):
    out_dir = tmp_path / "out"
    settings_file = tmp_path / "existing.json"
    settings_file.write_text('{"target": "hPDL1"}\n')

    task = BindCraft2(
        output_dir=out_dir,
        settings_file=settings_file,
        number_of_final_designs=2,
        binder_lengths=[70, 90],
        max_trajectories=20,
        campaign_seed=99,
        disulfide_staple=True,
        mixed_topology=True,
        settings_overrides={"desperation": False},
        device="gpu",
    )

    task._run_task = MagicMock()
    task.run()

    task._run_task.assert_called_once()
    kwargs = task._run_task.call_args.kwargs
    assert kwargs["settings_file"] == str(settings_file.resolve())
    assert "--disulfide-staple" in kwargs["property_flags"]
    assert "--mixed-topology" in kwargs["property_flags"]
    assert f"--set 'project_folder={out_dir.resolve()}'" in kwargs["set_flags"]
    assert "--set 'number_of_final_designs=2'" in kwargs["set_flags"]
    assert "--set 'binder_lengths=[70, 90]'" in kwargs["set_flags"]
    assert "--set 'max_trajectories=20'" in kwargs["set_flags"]
    assert "--set 'campaign_seed=99'" in kwargs["set_flags"]
    assert "--set 'desperation=false'" in kwargs["set_flags"]


def test_bindcraft2_file_not_found(tmp_path):
    out_dir = tmp_path / "out"

    # Nonexistent settings file
    t1 = BindCraft2(output_dir=out_dir, settings_file=tmp_path / "nonexistent.json")
    with pytest.raises(FileNotFoundError, match="Settings file not found"):
        t1.run()

    # Nonexistent target PDB
    t2 = BindCraft2(output_dir=out_dir, target_pdb=tmp_path / "nonexistent.pdb")
    with pytest.raises(FileNotFoundError, match="Target PDB file not found"):
        t2.run()


def test_bindcraft2_dry_run(tmp_path):
    out_dir = tmp_path / "out"
    pdb_file = tmp_path / "target.pdb"
    pdb_file.write_text("ATOM 1\n")

    task = BindCraft2(
        output_dir=out_dir,
        target_pdb=pdb_file,
        number_of_final_designs=2,
        dry_run=True,
    )
    task.run()

    assert (out_dir / "3_Ranked" / "!_Ranked.csv").exists()
    assert (out_dir / "3_Ranked" / "mock_design_0001_seq1.cif").exists()
    assert (out_dir / "3_Ranked" / "mock_design_0002_seq1.cif").exists()
    assert (out_dir / "2_Refolded" / "!_Refolded.csv").exists()
    assert (out_dir / "1_Trajectories" / "!_Trajectories.csv").exists()
    assert (out_dir / "summary.csv").exists()
    assert (out_dir / "campaign_metadata.json").exists()


def test_bindcraft2_ribbon_module_import():
    assert hasattr(ribbon, "BindCraft2")
    assert ribbon.BindCraft2 is BindCraft2


def test_bindcraft2_gpu_smoke(tmp_path, request):
    if not request.config.getoption("--run-gpu-smoke"):
        pytest.skip("GPU smoke disabled. Use --run-gpu-smoke to enable.")

    hpdl1_pdb = BINDCRAFT2_INPUTS / "hPDL1.pdb"
    assert hpdl1_pdb.exists(), f"hPDL1 PDB missing at {hpdl1_pdb}"

    out_dir = tmp_path / "bindcraft2_gpu_out"

    task = BindCraft2(
        output_dir=out_dir,
        target_pdb=hpdl1_pdb,
        chains="A",
        hotspots="54,56,66,115",
        modality="binder",
        binder_lengths=[35, 35],
        number_of_final_designs=1,
        max_trajectories=1,
        screen_steps=2,
        refine_steps=2,
        anneal_steps=2,
        harden_steps=1,
        mutate_steps=2,
        design_recycles=1,
        validation_recycles=1,
        campaign_seed=42,
        device="gpu",
    )
    task.run()

    assert (out_dir / "summary.csv").exists(), "summary.csv missing after BindCraft2 run"
    assert (out_dir / "summary.csv").stat().st_size > 0, "summary.csv is empty"
    assert (out_dir / "campaign_metadata.json").exists(), "campaign_metadata.json missing"
    assert (out_dir / "1_Trajectories").exists(), "1_Trajectories directory missing"
