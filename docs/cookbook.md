# Cookbook

!!! warning 
    🚧 This page is under construction! We are still adding recipes for several `Tasks`. 🚧

## Chai-1
Chai-1 is a high accuracy, ligand-aware protein folding tool.
It takes as input a FASTA file, and outputs a directory of folded PDB structures.
Chai-1 scores are output as .npz files.

### Simple
Fold a protein from a sequence
``` python
ribbon.Chai1(
        fasta_file = 'my_sequence.fasta',   # A single input FASTA. If there are multiple sequences, they will be folded in the same structure.
        output_dir = './out'                # Where the outputs will be stored
    ).run()
```

### Advanced
Fold a protein from a sequence. Include two copies of ligand 
``` python
ribbon.Chai1(
        fasta_file = 'my_sequence.fasta',   # A single input FASTA. If there are multiple sequences, they will be folded in the same structure.
        output_dir = './out',               # Where the outputs will be stored
        smiles_string = 'C1=CC=CC=C1',      # SMILES string of our ligand
        num_ligands = 2,                    # How many copies of our ligand?
        device = 'gpu'                      # Run on GPU (necessary for Chai-1)
    ).run()
```

## LigandMPNN
LigandMPNN is a ligand-aware tool for designing sequences for a backbone structure.
It takes as input a list of PDB files, and outputs a sequence (or sequences) that are predicted to fold into that structure.
While LigandMPNN is gpu-accelerated, it seems to run fast on the CPU as well.

### Simple
Design a sequence for a backbone:
``` python
ribbon.LigandMPNN(
        structure_list = ['my_structure.pdb'],  # List of PDB files
        output_dir = './out'                    # Output directory
        num_designs = 5                         # How many sequences should we generate?
    ).run()
```
The output folder will have the following structure:
``` python
output_dir/
├─ backbones/    # Backbone structures with labeled AAs (but no sidechains)
├─ packed/       # (Optional) Backbone structures with packed sidechains
├─ seqs/         # A single FASTA containing the reference sequence and all designed sequences
└─ seqs_split/   # Individual FASTA files for each designed sequence
```

### Advanced
Design a homodimeric sequence, keeping crucial residues fixed.
This example uses the `extra_args` parameter to add extra parameters into your run command.
Note that this can inject arbitrary code into your container - use with caution!
```python
ligandmpnn_task = ribbon.LigandMPNN(
        structure_list = ['my_structure.pdb'],  # List of PDB files
        output_dir = './out'                    # Output directory
        num_designs = 5                         # How many sequences should we generate?
		extra_args= '--fixed_residues \"' + RESIDUES_TO_KEEP + '\" --homo_oligomer 1'	# Make sure to keep my catalytric residues, and make two chains identical.
	)
```

## RaptorX-Single
RaptorX-Single is a fast protein folding tool. It can fold small structures in as little as 5 seconds, after an initial loading period.
It takes as input a FASTA file (or directory containing FASTAs), and outputs a directory of folded PDB structures.

### Simple
Fold a directory of FASTA files
```python
    ribbon.RaptorXSingle(
            fasta_file_or_dir = './my_FASTA_directory/',
            output_dir = './out'
    ).run()
```

### Advanced
Fold a directory of FASTA files using a non-default model checkpoint (param). Run on the CPU (slower; GPU is default).
```python
    ribbon.RaptorXSingle(
            fasta_file_or_dir = './my_FASTA_directory/',
            output_dir = './out',
            param = 'RaptorX-Single-ESM1b-ESM1v-ProtTrans-Ab.pt',
            device='cpu'
    ).run()
```
The available model parameters are:
```python
'RaptorX-Single-ESM1b.pt',
'RaptorX-Single-ESM1v.pt',
'RaptorX-Single-ProtTrans.pt',
'RaptorX-Single-ESM1b-ESM1v-ProtTrans.pt',
'RaptorX-Single-ESM1b-Ab.pt',
'RaptorX-Single-ESM1v-Ab.pt',
'RaptorX-Single-ProtTrans-Ab.pt',
'RaptorX-Single-ESM1b-ESM1v-ProtTrans-Ab.pt'
```

## PLACER
PLACER is an all-atom graph neural network that denoises atomic coordinates to predict ligand binding poses, dock ligands into protein pockets, predict protein sidechain conformations, and generate conformational ensembles.
It takes a PDB/mmCIF structure (with ligand present in the input file or defined via ligand reference files) and outputs a multimodel PDB of predicted conformers (`<label>_model.pdb`) along with confidence and geometric scores in a CSV file (`<label>.csv`).

### Simple
Predict ligand binding conformations in a protein-ligand complex and rerank poses by predicted RMSD (pRMSD):
```python
ribbon.PLACER(
    input_file="protein_ligand.pdb",
    output_dir="./placer_out",
    num_samples=10,
    predict_ligand="D-LDP-501",           # Ligand ID to predict and dock
    rerank="prmsd",                       # Rerank output models by predicted RMSD (prmsd, plddt, or plddt_pde)
    device="gpu"
).run()
```

### Advanced
Predict multiple ligands simultaneously in a pocket with reference topology definitions, mutating specific pocket positions, adding covalent bond constraints, and fixing a co-factor:
```python
ribbon.PLACER(
    input_file="complex.pdb",
    output_dir="./placer_advanced_out",
    num_samples=50,
    predict_ligand=["D-LDP-501", "C-HEM-500"],  # Predict inhibitor and heme together
    predict_multi=True,                         # Dock and score multiple ligands simultaneously
    ligand_file={                               # Provide SDF/MOL2 files for exact hybridization/bonds
        "HEM": "ligands/HEM.mol2",
        "LDP": "ligands/LDP.sdf"
    },
    mutate={"128A": "TRP", "130A": "GLY"},      # Mutate pocket residues
    bonds=["A-42-CYS-SG:D-501-LDP-C1:1.8"],     # Covalent bond constraint (atom1:atom2:distance)
    target_res="A-149",                         # Center crop sphere on a specific residue
    rerank="prmsd",
    device="gpu"
).run()
```

## ESMFold2
ESMFold2 is an all-atom co-folding tool for proteins, nucleic acids, and small molecule ligands developed by the Chan Zuckerberg Biohub. It pairs the 6B-parameter ESMC protein language model with an all-atom diffusion trunk.
It accepts FASTA files, sequence strings, or multimers (`:` separated chains), and outputs predicted mmCIF structures (`<name>__<cond>__s<seed>__m<k>.cif`) alongside a `summary.csv` and `manifest.jsonl` recording pLDDT, pTM, ipTM, and runtime metrics.

### Simple
Fold a protein sequence in apo mode using the distilled ESMFold2-Fast model:
```python
ribbon.ESMFold2(
    sequence="MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    output_dir="./esmfold2_out",
    model="fast",                         # 'fast' (ESMFold2-Fast) or 'full' (ESMFold2 6B trunk)
    num_diffusion_samples=2,              # Number of diffusion samples to generate
    apo=True,
    device="gpu"
).run()
```

### Advanced
Co-fold a multi-chain protein complex with both small-molecule ligands and cofactors across multiple random seeds, using sequence chunking for memory efficiency:
```python
ribbon.ESMFold2(
    fasta_file="heterodimer.fasta",       # FASTA containing sequence(s), with ':' separating chains for complexes
    output_dir="./esmfold2_complex_out",
    model="full",                         # Full 6B-parameter ESMFold2 model
    ligands={                             # Define ligands by CCD code or SMILES string
        "BTN": "ccd:BTN",                 # Biotin from Chemical Component Dictionary
        "INHIBITOR": "smiles:Cc1cc(=O)oc2c1ccc1nn[nH]c12"
    },
    apo=True,                             # Also evaluate the apo state in the same batch run
    num_diffusion_samples=5,              # 5 structural samples per condition
    num_loops=10,                         # Recycling loops
    num_sampling_steps=68,                # Diffusion sampling steps
    chunk_size=128,                       # Chunk size for long targets to fit in GPU memory
    seeds=[42, 100],                      # Multiple random seeds for sampling diversity
    device="gpu"
).run()
```
