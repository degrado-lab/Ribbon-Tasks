# Smoke Case Layout

Each smoke case should have its own folder:

`tests/smoke/cases/<case_name>/`

Recommended structure:

1. `inputs/`: task-specific input files.
2. `expected/`: expected output files or reference values used by validators.

Example:

`tests/smoke/cases/calculate_distance/inputs/minimal.pdb`  
`tests/smoke/cases/calculate_distance/expected/distance.txt`

Reference each case folder in `tests/smoke/manifest.py` via `case_dir`.
