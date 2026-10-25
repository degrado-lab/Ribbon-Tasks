import json
import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest

import ribbon
from ribbon_tasks.tasks import LigandMPNN


def test_ligandmpnn_init_defaults(tmp_path):
    task = LigandMPNN(
        output_dir=tmp_path / "out",
        structure_list=["dummy.pdb"],
    )
    assert task.task_name == "LigandMPNN"
    assert task.num_designs == 1
    assert task.temperature == 0.1
    assert task.device == "cpu"
    assert task.model_type is None
    assert task.fixed_residues is None
    assert task.pack_side_chains is None
    assert task.dry_run is False


def test_ligandmpnn_alias_and_casing(tmp_path):
    task = LigandMPNN(
        output_dir=tmp_path / "out",
        structure_list=["dummy.pdb"],
        batch_size=5,
        bias_aa={"A": -1.0},
        omit_aa=["C", "M"],
        bias_aa_per_residue={"A1": {"G": 0.5}},
    )
    assert task.num_designs == 5
    assert task.bias_AA == {"A": -1.0}
    assert task.omit_AA == ["C", "M"]
    assert task.bias_AA_per_residue == {"A1": {"G": 0.5}}

    with pytest.raises(TypeError, match="unexpected keyword argument"):
        LigandMPNN(
            output_dir=tmp_path / "out",
            structure_list=["dummy.pdb"],
            invalid_argument_xyz=123,
        )


def test_ligandmpnn_run_compiles_all_flags(tmp_path):
    out_dir = tmp_path / "out"
    pdb_file = tmp_path / "input.pdb"
    pdb_file.write_text("ATOM      1  N   ALA A   1\nTER\nEND\n")

    task = LigandMPNN(
        output_dir=out_dir,
        structure_list=[pdb_file],
        num_designs=3,
        number_of_batches=2,
        temperature=0.2,
        model_type="ligand_mpnn",
        checkpoint_protein_mpnn=tmp_path / "pmpnn.pt",
        checkpoint_ligand_mpnn=tmp_path / "lmpnn.pt",
        checkpoint_per_residue_label_membrane_mpnn=tmp_path / "mem_res.pt",
        checkpoint_global_label_membrane_mpnn=tmp_path / "mem_glob.pt",
        checkpoint_soluble_mpnn=tmp_path / "sol.pt",
        checkpoint_path_sc=tmp_path / "sc.pt",
        fasta_seq_separation=":",
        verbose=1,
        fixed_residues=["A1", "A2"],
        fixed_residues_multi={str(pdb_file): ["A1", "A2"]},
        redesigned_residues=["A3", "A4"],
        redesigned_residues_multi={str(pdb_file): ["A3", "A4"]},
        bias_AA={"A": -1.5, "P": 2.0},
        bias_AA_per_residue={"A1": {"G": -0.3}},
        bias_AA_per_residue_multi={str(pdb_file): {"A1": {"G": -0.3}}},
        omit_AA=["C", "W"],
        omit_AA_per_residue={"A1": "PG"},
        omit_AA_per_residue_multi={str(pdb_file): {"A1": "PG"}},
        symmetry_residues=[["A1", "A2"], ["A3", "A4"]],
        symmetry_weights=[[1.0, 1.0], [2.0, 2.0]],
        homo_oligomer=True,
        file_ending="_test",
        zero_indexed=True,
        seed=42,
        save_stats=True,
        ligand_mpnn_use_atom_context=True,
        ligand_mpnn_cutoff_for_score=6.5,
        ligand_mpnn_use_side_chain_context=True,
        chains_to_design=["A"],
        parse_these_chains_only=["A", "B"],
        transmembrane_buried=["A1"],
        transmembrane_interface=["A2"],
        global_transmembrane_label=True,
        parse_atoms_with_zero_occupancy=True,
        pack_side_chains=True,
        number_of_packs_per_design=2,
        sc_num_denoising_steps=4,
        sc_num_samples=8,
        repack_everything=True,
        force_hetatm=True,
        packed_suffix="_custom_packed",
        pack_with_ligand_context=True,
        extra_args="--custom_extra_flag",
        device="cpu",
    )

    task._run_task = MagicMock()

    task.run()

    task._run_task.assert_called_once()
    kwargs = task._run_task.call_args[1]

    assert kwargs["output_dir"] == out_dir
    assert kwargs["num_designs"] == 3
    assert kwargs["temperature"] == 0.2
    assert kwargs["model_type_flag"] == "--model_type ligand_mpnn"
    assert kwargs["seed_flag"] == "--seed 42"
    assert kwargs["number_of_batches_flag"] == "--number_of_batches 2"
    assert kwargs["save_stats_flag"] == "--save_stats 1"
    assert kwargs["zero_indexed_flag"] == "--zero_indexed 1"
    assert kwargs["file_ending_flag"] == '--file_ending "_test"'
    assert kwargs["homo_oligomer_flag"] == "--homo_oligomer 1"
    assert kwargs["pack_side_chains_flag"] == "--pack_side_chains 1"
    assert kwargs["repack_everything_flag"] == "--repack_everything 1"
    assert kwargs["force_hetatm_flag"] == "--force_hetatm 1"
    assert kwargs["packed_suffix_flag"] == '--packed_suffix "_custom_packed"'
    assert kwargs["pack_with_ligand_context_flag"] == "--pack_with_ligand_context 1"
    assert kwargs["number_of_packs_per_design_flag"] == "--number_of_packs_per_design 2"
    assert kwargs["sc_num_denoising_steps_flag"] == "--sc_num_denoising_steps 4"
    assert kwargs["sc_num_samples_flag"] == "--sc_num_samples 8"
    assert kwargs["ligand_mpnn_use_atom_context_flag"] == "--ligand_mpnn_use_atom_context 1"
    assert kwargs["ligand_mpnn_cutoff_for_score_flag"] == "--ligand_mpnn_cutoff_for_score 6.5"
    assert kwargs["ligand_mpnn_use_side_chain_context_flag"] == "--ligand_mpnn_use_side_chain_context 1"
    assert kwargs["global_transmembrane_label_flag"] == "--global_transmembrane_label 1"
    assert kwargs["parse_atoms_with_zero_occupancy_flag"] == "--parse_atoms_with_zero_occupancy 1"
    assert kwargs["verbose_flag"] == "--verbose 1"
    assert kwargs["fasta_seq_separation_flag"] == '--fasta_seq_separation ":"'
    assert kwargs["fixed_residues_flag"] == '--fixed_residues "A1 A2"'
    assert kwargs["redesigned_residues_flag"] == '--redesigned_residues "A3 A4"'
    assert kwargs["chains_to_design_flag"] == '--chains_to_design "A"'
    assert kwargs["parse_these_chains_only_flag"] == '--parse_these_chains_only "A,B"'
    assert kwargs["transmembrane_buried_flag"] == '--transmembrane_buried "A1"'
    assert kwargs["transmembrane_interface_flag"] == '--transmembrane_interface "A2"'
    assert kwargs["omit_aa_flag"] == '--omit_AA "CW"'
    assert kwargs["bias_aa_flag"] == '--bias_AA "A:-1.5,P:2.0"'
    assert kwargs["symmetry_residues_flag"] == '--symmetry_residues "A1,A2|A3,A4"'
    assert kwargs["symmetry_weights_flag"] == '--symmetry_weights "1.0,1.0|2.0,2.0"'
    assert "--checkpoint_protein_mpnn" in kwargs["checkpoint_protein_mpnn_flag"]
    assert "--checkpoint_ligand_mpnn" in kwargs["checkpoint_ligand_mpnn_flag"]
    assert "--checkpoint_path_sc" in kwargs["checkpoint_path_sc_flag"]
    assert "--fixed_residues_multi" in kwargs["fixed_residues_multi_flag"]
    assert "--redesigned_residues_multi" in kwargs["redesigned_residues_multi_flag"]
    assert "--bias_AA_per_residue" in kwargs["bias_aa_per_residue_flag"]
    assert "--bias_AA_per_residue_multi" in kwargs["bias_aa_per_residue_multi_flag"]
    assert "--omit_AA_per_residue" in kwargs["omit_aa_per_residue_flag"]
    assert "--omit_AA_per_residue_multi" in kwargs["omit_aa_per_residue_multi_flag"]
    assert kwargs["extra_args"] == "--custom_extra_flag"


def test_ligandmpnn_temp_files_cleaned_up(tmp_path):
    out_dir = tmp_path / "out"
    pdb_file = tmp_path / "input.pdb"
    pdb_file.write_text("ATOM      1  N   ALA A   1\nTER\nEND\n")

    task = LigandMPNN(
        output_dir=out_dir,
        structure_list=[pdb_file],
        fixed_residues_multi={str(pdb_file): ["A1"]},
        bias_AA_per_residue={"A1": {"G": 0.5}},
    )

    created_temp_files = []

    def fake_run_task(task_name, **kwargs):
        # Record temp files while task is executing
        for k, v in kwargs.items():
            if isinstance(v, str) and "/tmp/" in v:
                # Extract filename from flag string
                path_str = v.split('"')[1] if '"' in v else v
                created_temp_files.append(path_str)

    task._run_task = MagicMock(side_effect=fake_run_task)
    task.run()

    assert len(created_temp_files) >= 2
    # Verify that all temp files recorded during execution have been cleaned up afterwards
    for tf in created_temp_files:
        assert not os.path.exists(tf), f"Temp file {tf} was not cleaned up!"


def test_ligandmpnn_dry_run_generates_all_expected_artifacts(tmp_path):
    out_dir = tmp_path / "dry_out"
    pdb_file = tmp_path / "test_prot.pdb"
    pdb_file.write_text("ATOM      1  N   ALA A   1\nTER\nEND\n")

    task = LigandMPNN(
        output_dir=out_dir,
        structure_list=[pdb_file],
        num_designs=2,
        number_of_batches=2,
        pack_side_chains=True,
        number_of_packs_per_design=2,
        save_stats=True,
        dry_run=True,
    )
    task.run()

    # Total designs = 2 * 2 = 4
    seqs_dir = out_dir / "seqs"
    assert (seqs_dir / "test_prot.fa").exists()

    seqs_split_dir = out_dir / "seqs_split"
    assert seqs_split_dir.exists()
    split_files = list(seqs_split_dir.glob("*.fasta"))
    assert len(split_files) == 4

    backbones_dir = out_dir / "backbones"
    assert backbones_dir.exists()
    for i in range(1, 5):
        assert (backbones_dir / f"test_prot_{i}.pdb").exists()

    packed_dir = out_dir / "packed"
    assert packed_dir.exists()
    for d in range(1, 5):
        for p in range(1, 3):
            assert (packed_dir / f"test_prot_packed_{d}_{p}.pdb").exists()

    stats_dir = out_dir / "stats"
    assert stats_dir.exists()
    assert (stats_dir / "test_prot.pt").exists()


@pytest.mark.smoke
def test_ligandmpnn_real_execution_with_new_features(tmp_path, request):
    """
    Integration smoke test running LigandMPNN with multiple new features against real container on GPU.
    """
    if not request.config.getoption("--run-gpu-smoke"):
        pytest.skip("GPU smoke disabled. Use --run-gpu-smoke to enable.")

    minimal_pdb = Path(__file__).parent / "smoke" / "cases" / "ligand_mpnn" / "inputs" / "minimal.pdb"
    out_dir = tmp_path / "real_out"

    task = ribbon.LigandMPNN(
        output_dir=out_dir,
        structure_list=[minimal_pdb],
        num_designs=2,
        seed=111,
        fixed_residues=["A1", "A2"],
        bias_AA={"A": 2.0},
        omit_AA=["C"],
        pack_side_chains=True,
        number_of_packs_per_design=1,
        save_stats=True,
        device="gpu",
    )
    task.run()

    assert (out_dir / "seqs_split").exists()
    assert len(list((out_dir / "seqs_split").glob("*.fasta"))) == 2
    assert (out_dir / "backbones").exists()
    assert (out_dir / "packed").exists()
    assert (out_dir / "stats").exists()
    assert (out_dir / "stats" / "minimal.pt").exists()
