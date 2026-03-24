# Contributing to Ribbon-Tasks

This repository defines task interfaces for Ribbon.
The main goal is to let users switch Ribbon-Task versions without changing the Ribbon runner itself.

## What Lives Where

- `ribbon_tasks/containers.json`: Container registry for tasks. Maps a container name to `[local_sif_name, oras_image_ref]`.
- `ribbon_tasks/tasks.json`: Declarative task specs. Each task defines `container`, `command`, and optional environment variables.
- `ribbon_tasks/tasks.py`: Python wrappers that present user-facing APIs and call `self._run_task(...)`.
- `ribbon_tasks/task_scripts/`: Helper scripts used when a tool needs input/output adaptation beyond a single CLI call.

## Task Architecture

Each task is defined in layers:

1. A container alias in `containers.json`.
2. A command template in `tasks.json`.
3. A Python wrapper class in `tasks.py`.
4. Optional script(s) in `task_scripts/` for custom preprocessing/postprocessing.

The wrapper class sets `self.task_name` to match a key in `tasks.json` and passes formatted arguments to `self._run_task(...)`.

## Add a New Task

1. Add or reuse a container in `ribbon_tasks/containers.json`.
2. Add a task entry in `ribbon_tasks/tasks.json`.
3. Add a wrapper class in `ribbon_tasks/tasks.py`.
4. Add helper script(s) under `ribbon_tasks/task_scripts/<tool>/` if needed.
5. Verify naming consistency across all layers.

## Minimal Wrapper Template

```python
from ribbon.runner import Task
from ribbon.utils import make_directory
from pathlib import Path

class MyTask(Task):
    def __init__(self, input_file, output_dir=".", device="cpu", extra_args=""):
        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "MyTask"  # Must match tasks.json key
        self.input_file = input_file
        self.output_dir = output_dir

    def run(self):
        self.output_dir = make_directory(self.output_dir)
        self._run_task(
            self.task_name,
            input_file=self.input_file,
            output_dir=str(self.output_dir),
            extra_args=self.extra_args,
            device=self.device,
        )
```

## `tasks.json` Template

```json
"MyTask": {
  "name": "MyTask",
  "description": "Short description of what this task does.",
  "container": "MyContainer",
  "command": "python3 tool.py {input_file} --out {output_dir} {extra_args}",
  "requirements": "None"
}
```

## Consistency Checklist

Before opening a PR, confirm all items below:

1. `tasks.py` wrapper `self.task_name` exactly matches a key in `tasks.json`.
2. `tasks.json` `container` value exists in `containers.json`.
3. Every placeholder in `tasks.json` `command` is provided by wrapper `_run_task(...)`.
4. Wrapper and task script argument names are aligned.
5. Output directories/files are created before command execution when needed.
6. Any optional flags are emitted only when set.
7. New tasks are importable from `ribbon_tasks/__init__.py` (currently `from .tasks import *`).

## Common Pitfalls

1. Naming drift between wrapper class and `tasks.json` key.
2. Declaring a task in `tasks.json` without a usable wrapper (or vice versa).
3. Mutable default arguments like `[]` in function signatures.
4. Fragile shell-string command assembly when arguments contain spaces/special characters.
5. Missing docs/examples for expected input formats (FASTA, PDB, atom selectors, etc.).

## Code Style Notes

1. Keep wrappers thin: validation + argument shaping + `_run_task(...)`.
2. Put heavy formatting logic into `task_scripts/`.
3. Prefer explicit argument names over positional coupling.
4. Use `Path` and `make_directory(...)` for filesystem safety.

## I/O Standard (Recommended)

Use these conventions for wrapper arguments unless a tool has a hard requirement that conflicts:

1. Inputs:
- `input_file`: one input file
- `input_dir`: one input directory
- `input_files`: multiple input files
- Use domain-specific names only when they carry important meaning (for example `fasta_file`).
2. Outputs:
- `output_file`: one output artifact
- `output_dir`: directory output
- `output_prefix`: basename/prefix for multiple generated artifacts
3. Optional runtime args:
- `device` and `extra_args` should be named consistently when supported.
4. Validation:
- If a task supports file-or-dir input, validate that exactly one is provided.
- Validate suffixes for output files when format matters.
5. Filesystem behavior:
- Normalize paths in `__init__`.
- Create output directories in `run()`.
- Keep output path values as pure paths (do not embed CLI flags in them).

Suggested `run()` flow:

1. Validate/normalize input mode.
2. Create output target(s).
3. Build optional flags/derived args in separate local variables.
4. Call `_run_task(...)`.

Utilities for these patterns are in `ribbon_tasks/io_utils.py`.

## Manual Test Guidance

1. Run a minimal successful example for the new task.
2. Run one error-path example (invalid input or missing required argument).
3. Verify expected output files are created and named correctly.
4. If you added a helper script, run it directly once to verify CLI parsing.

## Smoke Test Harness

Container smoke tests live in `tests/smoke/` and are intended for pre-release validation.

1. `tests/smoke/manifest.py` defines task cases.
2. `tests/smoke/test_tasks_smoke.py` runs each case end-to-end with real containers.
3. `tests/smoke/cases/<case_name>/inputs/` holds task-specific inputs.
4. `tests/smoke/cases/<case_name>/expected/` holds expected outputs/reference values.
5. `tests/smoke/validators.py` defines reusable validation kinds.

Supported validation kinds:

1. `exists`
2. `nonempty`
3. `contains_text`
4. `exact_text`
5. `json_equals`
6. `csv_columns`
7. `numeric_tolerance`

Useful commands:

1. `pytest -m smoke`
2. `pytest -m smoke --smoke-task Reduce`
3. `pytest -m smoke --run-gpu-smoke`
4. `pytest -m smoke --smoke-timeout-multiplier 2`

## Pull Request Expectations

Include the following in the PR description:

1. What task(s) were added or changed.
2. Which container image tag(s) were introduced or updated.
3. Example invocation used for validation.
4. Any known limitations or follow-up work.
