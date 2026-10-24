---
name: script-ribbon-workflows
description: Script and orchestrate structural biology workflows using the Ribbon Python package and its containerized tasks.
---

# Skill: Script Ribbon Workflows

This skill guides AI coding agents in designing, scripting, and orchestrating computational structural biology and protein design workflows using the **Ribbon** python framework.

---

## 1. Context Discovery (Mandatory First Step)

When you are asked to script or modify Ribbon workflows in a user's workspace, **the source code of Ribbon-Tasks may not be visible in your primary workspace explorer**. Ribbon-Tasks is imported dynamically as a subpackage of the standard pip-installed `ribbon` package.

To bootstrap your understanding of the environment and find where all code and configurations reside, you **MUST** run the following commands:

### A. Locate Active Configuration & Paths
Run the Ribbon CLI info command (using the active Python virtual environment or conda environment, such as `~/miniforge3/envs/my_ribbon_env/bin/ribbon`):
```bash
ribbon info
```
If you do not know which environment has ribbon installed, query the user.
From the tabular output, identify and record:
*   **Tasks directory**: The exact absolute path where active task code resides (e.g., `~/.ribbon/ribbon_tasks`).
*   **Active tasks version**: The currently loaded version of task definitions (e.g., `v0.1.5`).

### B. List Registered Tasks
Run the Ribbon CLI list command:
```bash
ribbon list tasks
```
This prints the full set of registered biological tools (e.g., `LigandMPNN`, `RFDiffusionAA`, `FastRelax`, `Chai-1`, `EasyMD`).

---

## 2. Parameter Tracing (Do Not Hallucinate Signatures)

Because Ribbon dynamically loads task classes, you **must not guess or estimate task parameters**.

### PREFERRED METHOD: Inspect pydocs

1.  **Catalogue Classes**:
To list the active Task classes installed in the current version of Ribbon-Tasks, use the pydoc command: `python -m pydoc ribbon.ribbon_tasks.tasks | cat | head -n35` to list the beginning of the API reference for the `tasks.py` file.

2.  **Inspect Class Signatures**:
Determine which Tasks are required for the current programming task. Then, review the `__init__` constructor and docstrings for the specific task you want to invoke using the command `python -m pydoc ribbon.ribbon_tasks.tasks.<TaskName> | cat`.

For example, `python -m pydoc ribbon.ribbon_tasks.tasks.RFDiffusionAA | cat`


### ALTERNATE METHOD: Inspect Source code.

Once you have identified the active **Tasks directory** from `ribbon info` (translating any user `~` symbol to their absolute home path, e.g., `/home/username/`):

1.  **Inspect Class Signatures**:
    Open the Python definitions file at `<Tasks directory>/<version>/Ribbon-Tasks/ribbon_tasks/tasks.py` using file-viewer tools. Review the `__init__` constructor and docstrings for the specific task you want to invoke (e.g., check `class RFDiffusionAA(Task)` to see exactly what keys it expects).
2.  **Inspect Command & Container Layouts**:
    Open `<Tasks directory>/<version>/Ribbon-Tasks/ribbon_tasks/tasks.json` to see how arguments are mapped to terminal strings inside containers.
3.  **Inspect Container registries**:
    Open `<Tasks directory>/<version>/Ribbon-Tasks/ribbon_tasks/containers.json` to review the Apptainer `.sif` mappings and ORAS download urls.

---

## 3. Scripting & Implementation Protocols

When writing or editing Python scripts to automate Ribbon workflows:

### A. Root Module Import Model
Always import `ribbon` directly, and instantiate tasks as attributes of the root `ribbon` namespace. **Do not** import from submodules like `ribbon_tasks.tasks`.
```python
import ribbon

# Preferred syntax:
relax_task = ribbon.FastRelax(pdb_input_file="input.pdb", output_dir="out/")
```

### B. Implement Pipeline-Level `dry_run`
Always add a `dry_run` toggle at the entry-point function of any pipeline script, and propagate it to every task constructor:
```python
def my_pipeline(input_pdb, output_dir, dry_run=False):
    # Pass dry_run flag to each stage
    prep = ribbon.RosettaLigandPrepare(input_file=..., dry_run=dry_run)
    prep.run()
```
*   **Dry-run behavior**: Instantiating with `dry_run=True` generates mock file outputs (PDBs, CIFs, Fastas, score logs) instantly, validating parameter formatting and directory routing without invoking container download processes, GPU runtimes, or HPC scheduling.

### C. Local vs. Cluster Queueing
*   Use `.run()` on a task instance to execute it immediately using the local machine's Apptainer engine.
*   Use `.queue(scheduler="SLURM", time="H:MM:SS", mem="XXG", gpus=X)` to serialize the task and schedule it dynamically on cluster nodes.

### D. Best Practices
* For more detailed instructions on scripting, read [RIBBON WORKFLOW SCRIPTING GUIDE](RIBBON_WORKFLOW_SCRIPTING_GUIDE.md).

---

## 4. Customizing or Extending Ribbon-Tasks
If the user asks you to add support for a new biological tool:
1.  Navigate to the active tasks folder found via `ribbon info`.
2.  Register the new container image path in `containers.json`.
3.  Add the task's command CLI template in `tasks.json`.
4.  Subclass `Task` and implement `__init__`, `_run_dry()`, and `run()` in `tasks.py`.
