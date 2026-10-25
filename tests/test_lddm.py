from pathlib import Path
from unittest.mock import MagicMock

import pytest

import ribbon
from ribbon_tasks.tasks import LDDM

TESTS_DIR = Path(__file__).parent
KRAS_INPUTS = TESTS_DIR / "smoke" / "cases" / "lddm" / "inputs"


def test_lddm_init_defaults(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "protein.pdb"
    ref = tmp_path / "ref.sdf"
    prot.write_text("ATOM 1\n")
    ref.write_text("$$$$\n")

    task = LDDM(
        output_dir=out_dir,
        protein=prot,
        ref_ligand=ref,
    )
    assert task.task_name == "LDDM"
    assert task.mode == "design"
    assert task.ligand is None
    assert task.checkpoint is None
    assert task.output_file == out_dir / "samples.sdf"
    assert task.n_samples == 10
    assert task.batch_size is None
    assert task.n_steps == 100
    assert task.sampler == "ForwardEuler"
    assert task.sampling_noise == 5.0
    assert task.molecule_size is None
    assert task.atoms_to_dock is None
    assert task.n_frames is None
    assert task.return_projected_final is False
    assert task.seed is None
    assert task.device == "gpu"
    assert task.verbose is False
    assert task.extra_args == ""
    assert task.dry_run is False


def test_lddm_aliases(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "kras.pdb"
    ref = tmp_path / "ref.sdf"
    frag = tmp_path / "frag.sdf"
    prot.write_text("ATOM 1\n")
    ref.write_text("$$$$\n")
    frag.write_text("$$$$\n")

    # Test protein_file / ref_ligand_file / ligand_file / num_samples / num_steps
    task1 = LDDM(
        output_dir=out_dir,
        protein_file=prot,
        ref_ligand_file=ref,
        ligand_file=frag,
        num_samples=3,
        num_steps=25,
    )
    assert task1.protein == prot
    assert task1.ref_ligand == ref
    assert task1.ligand == frag
    assert task1.n_samples == 3
    assert task1.n_steps == 25

    # Test pdb_file / fragment_file
    task2 = LDDM(
        output_dir=out_dir,
        pdb_file=prot,
        ref_ligand=ref,
        fragment_file=frag,
    )
    assert task2.protein == prot
    assert task2.ligand == frag


def test_lddm_validation_errors(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "prot.pdb"
    ref = tmp_path / "ref.sdf"

    # Missing protein
    with pytest.raises(ValueError, match="protein"):
        LDDM(output_dir=out_dir, ref_ligand=ref)

    # Missing ref_ligand
    with pytest.raises(ValueError, match="ref_ligand"):
        LDDM(output_dir=out_dir, protein=prot)

    # Invalid mode
    with pytest.raises(ValueError, match="Invalid mode"):
        LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, mode="unsupported")

    # Mode dock without ligand
    with pytest.raises(ValueError, match="dock.*ligand"):
        LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, mode="dock")

    # Unexpected keyword argument
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, invalid_param=123)


def test_lddm_run_compiles_design_mode(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "prot.pdb"
    ref = tmp_path / "ref.sdf"
    frag = tmp_path / "frag.sdf"
    out_sdf = tmp_path / "custom_out.sdf"
    prot.write_text("ATOM 1\n")
    ref.write_text("$$$$\n")
    frag.write_text("$$$$\n")

    task = LDDM(
        output_dir=out_dir,
        protein=prot,
        ref_ligand=ref,
        mode="design",
        ligand=frag,
        output_file=out_sdf,
        n_samples=5,
        batch_size=2,
        n_steps=40,
        sampler="HeunSampler",
        sampling_noise=3.5,
        molecule_size="uniform_15_25",
        n_frames=10,
        return_projected_final=True,
        seed=42,
        device="gpu",
        verbose=True,
        extra_args="--custom_flag 1",
    )

    task._run_task = MagicMock()
    task.run()

    task._run_task.assert_called_once_with(
        "LDDM",
        mode="design",
        protein=str(prot.resolve()),
        ref_ligand=str(ref.resolve()),
        checkpoint="/app/models/checkpoints/lddm.ckpt",
        output_file=str(out_sdf.resolve()),
        n_samples=5,
        n_steps=40,
        sampler="HeunSampler",
        sampling_noise=3.5,
        batch_size_flag="--batch_size 2",
        molecule_size_flag="--molecule_size uniform_15_25",
        ligand_flag=f'--ligand "{frag.resolve()}"',
        atoms_to_dock_flag="",
        n_frames_flag="--n_frames 10",
        return_projected_final_flag="--return_projected_final",
        seed_flag="--seed 42",
        device_flag="--device cuda:0",
        verbose_flag="--verbose",
        extra_args="--custom_flag 1",
        device="gpu",
    )


def test_lddm_run_compiles_dock_mode(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "prot.pdb"
    ref = tmp_path / "ref.sdf"
    lig = tmp_path / "lig.sdf"
    prot.write_text("ATOM 1\n")
    ref.write_text("$$$$\n")
    lig.write_text("$$$$\n")

    task = LDDM(
        output_dir=out_dir,
        protein=prot,
        ref_ligand=ref,
        mode="dock",
        ligand=lig,
        atoms_to_dock=[0, 1, 2, 5],
        device="cpu",
    )

    task._run_task = MagicMock()
    task.run()

    task._run_task.assert_called_once_with(
        "LDDM",
        mode="dock",
        protein=str(prot.resolve()),
        ref_ligand=str(ref.resolve()),
        checkpoint="/app/models/checkpoints/lddm.ckpt",
        output_file=str((out_dir / "samples.sdf").resolve()),
        n_samples=10,
        n_steps=100,
        sampler="ForwardEuler",
        sampling_noise=5.0,
        batch_size_flag="",
        molecule_size_flag="",
        ligand_flag=f'--ligand "{lig.resolve()}"',
        atoms_to_dock_flag="--atoms_to_dock 0 1 2 5",
        n_frames_flag="",
        return_projected_final_flag="",
        seed_flag="",
        device_flag="--device cpu",
        verbose_flag="",
        extra_args="",
        device="cpu",
    )


def test_lddm_checkpoints(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "prot.pdb"
    ref = tmp_path / "ref.sdf"
    custom_ckpt = tmp_path / "custom.ckpt"
    prot.write_text("ATOM 1\n")
    ref.write_text("$$$$\n")
    custom_ckpt.write_text("weights\n")

    # Built-in default
    t_def = LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, checkpoint="default")
    t_def._run_task = MagicMock()
    t_def.run()
    assert t_def._run_task.call_args.kwargs["checkpoint"] == "/app/models/checkpoints/lddm.ckpt"

    # Built-in CD+BB
    t_cdbb = LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, checkpoint="CD+BB")
    t_cdbb._run_task = MagicMock()
    t_cdbb.run()
    assert t_cdbb._run_task.call_args.kwargs["checkpoint"] == "/app/models/checkpoints/lddm_CDBB.ckpt"

    # Custom existing checkpoint
    t_cust = LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, checkpoint=custom_ckpt)
    t_cust._run_task = MagicMock()
    t_cust.run()
    assert t_cust._run_task.call_args.kwargs["checkpoint"] == str(custom_ckpt.resolve())

    # Custom non-existing checkpoint
    t_nonexist = LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, checkpoint=tmp_path / "absent.ckpt")
    with pytest.raises(FileNotFoundError, match="Checkpoint file not found"):
        t_nonexist.run()


def test_lddm_file_not_found_errors(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "prot.pdb"
    ref = tmp_path / "ref.sdf"
    lig = tmp_path / "lig.sdf"

    # Protein not found
    t1 = LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref)
    with pytest.raises(FileNotFoundError, match="Protein file not found"):
        t1.run()

    # Ref ligand not found
    prot.write_text("ATOM 1\n")
    t2 = LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref)
    with pytest.raises(FileNotFoundError, match="Reference ligand file not found"):
        t2.run()

    # Ligand not found
    ref.write_text("$$$$\n")
    t3 = LDDM(output_dir=out_dir, protein=prot, ref_ligand=ref, ligand=lig)
    with pytest.raises(FileNotFoundError, match="Ligand file not found"):
        t3.run()


def test_lddm_dry_run(tmp_path):
    out_dir = tmp_path / "out"
    prot = tmp_path / "prot.pdb"
    ref = tmp_path / "ref.sdf"
    prot.write_text("ATOM 1\n")
    ref.write_text("$$$$\n")

    task = LDDM(
        output_dir=out_dir,
        protein=prot,
        ref_ligand=ref,
        n_samples=3,
        n_frames=5,
        dry_run=True,
    )
    task.run()

    out_file = out_dir / "samples.sdf"
    assert out_file.exists()
    content = out_file.read_text()
    assert content.count("$$$$") == 3
    assert (out_dir / "samples.pdb").exists()


def test_lddm_ribbon_module_import():
    assert hasattr(ribbon, "LDDM")
    assert ribbon.LDDM is LDDM


def test_lddm_gpu_smoke_design(tmp_path, request):
    if not request.config.getoption("--run-gpu-smoke"):
        pytest.skip("GPU smoke disabled. Use --run-gpu-smoke to enable.")

    kras_pdb = KRAS_INPUTS / "kras.pdb"
    ref_ligand_sdf = KRAS_INPUTS / "kras_ref_ligand.sdf"

    assert kras_pdb.exists(), f"KRAS PDB missing at {kras_pdb}"
    assert ref_ligand_sdf.exists(), f"KRAS ref ligand missing at {ref_ligand_sdf}"

    out_dir = tmp_path / "lddm_gpu_design_out"
    out_file = out_dir / "generated.sdf"

    task = LDDM(
        output_dir=out_dir,
        protein=kras_pdb,
        ref_ligand=ref_ligand_sdf,
        mode="design",
        output_file=out_file,
        n_samples=1,
        n_steps=5,
        device="gpu",
    )
    task.run()

    assert out_file.exists(), f"Generated SDF missing: {out_file}"
    assert out_file.stat().st_size > 0, "Generated SDF is empty"
    content = out_file.read_text(encoding="utf-8", errors="ignore")
    assert "$$$$" in content, "Generated SDF missing $$$$ delimiter"


def test_lddm_gpu_smoke_dock(tmp_path, request):
    if not request.config.getoption("--run-gpu-smoke"):
        pytest.skip("GPU smoke disabled. Use --run-gpu-smoke to enable.")

    kras_pdb = KRAS_INPUTS / "kras.pdb"
    ref_ligand_sdf = KRAS_INPUTS / "kras_ref_ligand.sdf"

    assert kras_pdb.exists(), f"KRAS PDB missing at {kras_pdb}"
    assert ref_ligand_sdf.exists(), f"KRAS ref ligand missing at {ref_ligand_sdf}"

    out_dir = tmp_path / "lddm_gpu_dock_out"
    out_file = out_dir / "docked.sdf"

    task = LDDM(
        output_dir=out_dir,
        protein=kras_pdb,
        ref_ligand=ref_ligand_sdf,
        mode="dock",
        ligand=ref_ligand_sdf,
        output_file=out_file,
        n_samples=1,
        n_steps=5,
        device="gpu",
    )
    task.run()

    assert out_file.exists(), f"Docked SDF missing: {out_file}"
    assert out_file.stat().st_size > 0, "Docked SDF is empty"
    content = out_file.read_text(encoding="utf-8", errors="ignore")
    assert "$$$$" in content, "Docked SDF missing $$$$ delimiter"
