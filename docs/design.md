# Task Definitions

## What Are Ribbon Tasks?

Ribbon uses a companion repository called **Ribbon-Tasks**, which contains the `ribbon_tasks` Python module. This module defines custom `Task` objects that Ribbon can run—both locally and when queueing jobs remotely.

Ribbon dynamically imports this module at runtime, allowing you to define and reuse custom tasks flexibly. To point Ribbon to your copy of the Ribbon-Tasks repo, set the `$RIBBON_TASKS_DIR` environment variable. For details, see [Customizing Your Tasks](customizing_tasks.md).

---

## Structure of a Ribbon Task

A Ribbon Task is defined by three components in the `ribbon_tasks` module:

1. **`tasks.py`** — Python wrapper for each Task
2. **`tasks.json`** — Defines the shell command and container to use for each Task
3. **`containers.json`** — Maps container names to their local filenames and remote sources

Let’s walk through each component.

---

## 1. `tasks.py`: Writing the Python Wrapper

Each Task is defined as a class that inherits from `Task`. The class must:

- Define `__init__()` to accept user inputs
- Set the `task_name` to match its entry in `tasks.json`
- Call `super().__init__()` to initialize the base `Task` class
- Store user inputs as instance variables
- Implement a `run()` method that calls `_run_task()` with named arguments

Here's an example for the `Chai-1` Task:

<details>
<summary>Example: Chai-1 Task Class</summary>

```python
class Chai1(Task):
    def __init__(self, fasta_file, output_dir='.', smiles_string=None, num_ligands=1, device='gpu'):
        super().__init__()
        self.task_name = "Chai-1"
        self.fasta_file = fasta_file
        self.output_dir = output_dir
        self.smiles_string = smiles_string
        self.num_ligands = num_ligands
        self.device = device

    def run(self):
        self.output_dir = make_directory(self.output_dir)
        self._run_task(
            self.task_name,
            fasta_file=self.fasta_file,
            smiles_string=self.smiles_string,
            output_dir=str(self.output_dir),
            num_ligands=self.num_ligands,
            device=self.device
        )
```

</details>


### Notes on `_run_task()`

- This method executes the actual job by filling in the command template from `tasks.json` using the provided arguments.
- Only variables defined in `tasks.json` can be passed into `_run_task()`.
- You can optionally pass an `extra_args` string for additional software-specific flags. Be cautious: this can inject arbitrary code if misused.
- Any lightweight pre/post-processing (like file creation or renaming) can also be handled in `run()`. For heavier logic, use a separate script and invoke it from the command string.

---

## 2. `tasks.json`: Defining the Command and Environment

This JSON file defines each Task’s metadata, including:

- The name and description
- The container to use
- Environment variables (if any)
- The shell command to run

<details>
<summary>Example: Chai-1 Task Definition</summary>

```json
"Chai-1": {
  "name": "Chai-1",
  "description": "Run Chai-1 structure prediction on a FASTA file and SMILES.",
  "container": "Chai-1",
  "command": "python $RIBBON_TASKS_MODULE_DIR/task_scripts/chai-1/predict_structure.py {fasta_file} \"{smiles_string}\" {output_dir} --num_ligands {num_ligands} {extra_args}",
  "requirements": "None",
  "environment_variables": {
    "TRANSFORMERS_OFFLINE": "1"
  }
}
```

</details>

### Command String

The `command` field is a string template. Curly-brace `{}` placeholders correspond to argument names passed into `_run_task()`.

For example:
```python
_run_task(fasta_file=self.fasta_file, output_dir=self.output_dir)
```
...would fill in `{fasta_file}` and `{output_dir}` in the command string.

### Using Custom Scripts

Instead of calling software directly, your command can invoke a custom Python script stored under `task_scripts/<container_name>/`. For example:
```bash
python $RIBBON_TASKS_MODULE_DIR/task_scripts/chai-1/predict_structure.py ...
```
- `$RIBBON_TASKS_MODULE_DIR` is automatically set within the container to the location of the `ribbon_tasks` module (which is within the `Ribbon-Tasks` repo).

This approach keeps your Task class clean while allowing complex logic to live in a dedicated script.

---

## 3. `containers.json`: Mapping Container Names

This file maps container names to a list containing:

1. The local `.sif` filename (used when the container is downloaded)
2. The remote location (usually via ORAS / DockerHub)

<details>
<summary>Example: Chai-1 Container</summary>

```json
"Chai-1": [
  "chai-1_0.5.2.sif",
  "oras://docker.io/nicholasfreitas/chai-1:0.5.2"
]
```

</details>

We recommend versioning containers explicitly to avoid compatibility issues—don’t use `latest`. Both the remote tag and local filename should reflect the version number.
