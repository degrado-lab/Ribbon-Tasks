---
name: develop-ribbon-tasks
description: Guide for creating, enhancing, configuring, and testing containerized Task classes in Ribbon-Tasks.
---

# Skill: Develop Ribbon Tasks

This skill guides AI coding agents in designing, implementing, registering, and testing containerized **Task** classes within the **Ribbon-Tasks** repository for the **Ribbon** python framework.

---

## 1. Architectural Overview

Every task in Ribbon is defined by a cohesive **triad** across three files in `Ribbon-Tasks/ribbon_tasks/`:

```
Ribbon-Tasks/
├── ribbon_tasks/
│   ├── containers.json    # 1. Container image & registry mappings
│   ├── tasks.json         # 2. Command templates & task metadata
│   ├── tasks.py           # 3. Python Task subclasses & logic
│   └── __init__.py        # Re-exports all Task classes for ribbon root import
└── tests/
    ├── smoke/
    │   ├── cases/<task>/inputs/  # Test inputs & fixtures
    │   └── manifest.py           # Automated smoke test registry
    └── test_<task>.py            # Dedicated unit & GPU integration tests
```

### The Three Components

1. **`containers.json`**:
   Maps an identifier (e.g. `"BindCraft2"`) to its Apptainer `.sif` image filename on disk and its remote ORAS registry pull URL:
   ```json
   "BindCraft2": [
       "bindcraft2_1.0.2.sif",
       "oras://docker.io/nicholasfreitas/bindcraft2:1.0.2"
   ]
   ```
   > [!IMPORTANT]
   > The `.sif` filename in `containers.json` must **exactly match** the file in `~/.ribbon/ribbon_containers/`. Watch out for hyphens vs. underscores (e.g. `bindcraft2_1.0.2.sif`).

2. **`tasks.json`**:
   Defines the task name, description, container linkage, and parameterized terminal command template:
   ```json
   "BindCraft2": {
       "name": "BindCraft2",
       "description": "Design protein binders, miniproteins, VHHs, and peptides using BindCraft2.",
       "container": "BindCraft2",
       "command": "bindcraft design {settings_file} {core_flag} {modality_flag} {property_flags} {set_flags} {metadata_flag} {extra_args}",
       "requirements": "None"
   }
   ```
   - Placeholder `{param}` strings are formatted dynamically by `Task._run_task(**kwargs)`.
   - Prefer passing optional CLI flags as compiled flag strings (e.g. `{property_flags}`, `{batch_size_flag}`) rather than hardcoding static flags, allowing Python to omit them when `None`.

3. **`tasks.py`**:
   Implements a subclass of `ribbon.runner.Task`. Each subclass must:
   - Set `self.task_name` matching the key in `tasks.json`.
   - Expose explicit, typed arguments for all upstream tool options (avoid dumping everything into `extra_args`).
   - Accept common aliases and normalize them in `__init__`.
   - Implement `_run_dry(self)` to generate mock outputs instantly without container execution.
   - Implement `run(self)` to validate inputs, build dynamic configs/flags, resolve file paths, and invoke `self._run_task(...)`.

---

## 2. Environment & Workspace Discovery

Ribbon-Tasks are run with the parent `ribbon` python package. Ensure package `ribbon` is installed in your python environment, or query the user if you need the proper environment.

Before writing or editing code:

1. **Find Active Ribbon Environment & Tasks Version**:
   ```bash
   ribbon info
   ```
   Identify:
   - **Active tasks version**: e.g. `v0.1.6`.
   - **Tasks directory**: e.g. `/home/username/.ribbon/ribbon_tasks/v0.1.6/Ribbon-Tasks`.
   - Ensure you are working in the active version's repository.

2. **Inspect the Target Container**:
   Containers reside in `~/.ribbon/ribbon_containers/`. Inspect the container before writing Python code:
   ```bash
   # Check container metadata and labels
   apptainer inspect ~/.ribbon/ribbon_containers/<image>.sif

   # Check container runscript & default environment
   apptainer inspect --runscript ~/.ribbon/ribbon_containers/<image>.sif
   apptainer exec ~/.ribbon/ribbon_containers/<image>.sif env

   # Check entrypoints, executables, baked models, and test fixtures
   apptainer exec ~/.ribbon/ribbon_containers/<image>.sif ls -la /app
   apptainer exec ~/.ribbon/ribbon_containers/<image>.sif which <executable>
   apptainer exec ~/.ribbon/ribbon_containers/<image>.sif <executable> --help
   ```

---

## 3. Step-by-Step Task Development Workflow

### Step 1: Register Container in `containers.json`
Add the container key to `ribbon_tasks/containers.json`:
```json
"MyTool": [
    "mytool_1.0.0.sif",
    "oras://docker.io/organization/mytool:1.0.0"
]
```

### Step 2: Register Command Template in `tasks.json`
Add the command format to `ribbon_tasks/tasks.json`:
- Use clean flag placeholders like `{output_flag}`, `{device_flag}`, `{extra_args}`.
- Always append `{extra_args}` at the end to allow user escape-hatch CLI parameters.

### Step 3: Implement Task Subclass in `tasks.py`

Attempt to port all the functionality of the original tool to the Ribbon-Tasks framework. Look at the original tool's source code and documentation for help. Ask the user for a github link to the tool, for the most thorough documentation.

Follow this canonical template:

```python
class MyTool(Task):
    """
    Class docstring describing the tool and its primary scientific application.
    """

    def __init__(
        self,
        output_dir: Union[str, Path],
        input_pdb: Optional[Union[str, Path]] = None,
        num_designs: int = 1,
        mode: str = "design",
        device: str = "gpu",
        extra_args: str = "",
        dry_run: bool = False,
        **kwargs,
    ):
        super().__init__(device=device, extra_args=extra_args)
        self.task_name = "MyTool"

        # 1. Alias Handling
        if input_pdb is None:
            input_pdb = kwargs.pop("pdb_file", None) or kwargs.pop("protein", None)
        if "n_designs" in kwargs:
            num_designs = kwargs.pop("n_designs")

        # 2. Reject Unexpected Arguments
        if kwargs:
            raise TypeError(f"MyTool.__init__() got unexpected keyword argument(s): {list(kwargs.keys())}")

        # 3. Argument Validation
        if input_pdb is None:
            raise ValueError("Argument 'input_pdb' is required for MyTool.")
        if mode not in ["design", "predict"]:
            raise ValueError(f"Invalid mode '{mode}'. Supported: 'design', 'predict'.")

        self.output_dir = Path(output_dir)
        self.input_pdb = Path(input_pdb)
        self.num_designs = num_designs
        self.mode = mode
        self.device = device
        self.dry_run = dry_run

    def _run_dry(self):
        """Generate realistic mock outputs for dry-run validation."""
        self.output_dir = make_directory(self.output_dir)
        
        # Create expected directory hierarchy and mock files
        out_file = self.output_dir / "result.pdb"
        with open(out_file, "w") as f:
            f.write("HEADER    MOCK DRY RUN STRUCTURE\nATOM      1  CA  ALA A   1       0.000   0.000   0.000\nEND\n")
        
        with open(self.output_dir / "summary.csv", "w") as f:
            f.write("design,score\nmock_001,0.95\n")

    def run(self):
        """Execute task or simulate dry run."""
        if self.dry_run:
            self._run_dry()
            return

        self.output_dir = make_directory(self.output_dir)

        # 1. Verify and resolve file paths
        input_path = self.input_pdb.resolve()
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_path}")

        # 2. Compile conditional CLI flags
        device_flag = "--device cuda:0" if self.device in ["gpu", "gpu_wsl"] else "--device cpu"
        ribbon_device = self.device if self.device in ["gpu", "gpu_wsl", "cpu"] else "gpu"

        # 3. Invoke Task runner
        self._run_task(
            self.task_name,
            input_pdb=str(input_path),
            output_dir=str(self.output_dir.resolve()),
            num_designs=self.num_designs,
            device_flag=device_flag,
            extra_args=self.extra_args,
            device=ribbon_device,
        )
```

### Step 4: Setup Fixtures & Smoke Manifest
1. Copy small, fast test structures (e.g. KRAS, PD-L1) to:
   `tests/smoke/cases/<task_name>/inputs/`
2. Register kwargs builder and validators in `tests/smoke/manifest.py`:
   ```python
   def _mytool_case(case_dir: Path, out_dir: Path) -> dict[str, Any]:
       return {
           "output_dir": str(out_dir),
           "input_pdb": str(case_dir / "inputs" / "minimal.pdb"),
           "num_designs": 1,
           "device": "gpu",
       }

   # In SMOKE_CASES list:
   {
       "name": "MyTool",
       "class_name": "MyTool",
       "case_dir": CASES_DIR / "mytool",
       "kwargs_builder": _mytool_case,
       "validations": [
           {"kind": "exists", "path": "result.pdb"},
           {"kind": "nonempty", "path": "result.pdb"},
       ],
       "timeout_s": 300,
       "requires_gpu": True,
   },
   ```

### Step 5: Write Pytest Test Suite (`tests/test_<task_name>.py`)
Always create a dedicated test file covering:
- **Default arguments**: Verify all default attributes.
- **Aliases & casing**: Test alternative argument names.
- **Input validation**: Missing inputs, invalid modes, unexpected kwargs.
- **Flag compilation**: Use `unittest.mock.MagicMock()` on `_run_task` to assert exact kwargs and flags passed.
- **File not found checks**: Nonexistent inputs raise `FileNotFoundError`.
- **Dry-run generation**: Verify mock output files are created with expected contents.
- **Root module import**: Assert `hasattr(ribbon, "MyTool")`.
- **GPU smoke test**: Guarded with `@pytest.mark.skipif(not request.config.getoption("--run-gpu-smoke"), ...)`.

---

## 4. Key Conventions & Traps to Avoid

### 1. Device Translation
- Ribbon's internal runner expects `device` to be `'gpu'`, `'cpu'`, or `'gpu_wsl'`.
- If `device == 'gpu'`, Ribbon passes `--nv` to Apptainer.
- Inside the container, tools usually expect `--device cuda:0` (or `cuda`) vs `--device cpu`.
- **Trap**: Do not pass `--device gpu` to tools that require PyTorch device strings (`cuda:0` / `cpu`). Map them separately in `run()`.

### 2. Path Resolution & Bind Mounts
- Apptainer automatically mounts the current working directory and `/tmp`.
- Always call `.resolve()` on file paths (`Path(f).resolve()`) so absolute host paths are passed to the container runscript.

### 3. JAX / AlphaFold Cache Pre-warming
- Tools utilizing JAX (e.g. BindCraft2, ColabDesign) compile computation graphs on first run.
- Set `JAX_COMPILATION_CACHE_DIR=/tmp/bindcraft_xla_cache` so repeated test runs do not recompile graphs from scratch.

### 4. Background Processes & Asynchronous Testing
- When testing live container runs in agent sessions, tools may execute as background tasks.
- Wait for task completion notifications before verifying output directories.
- Keep test sample/iteration counts minimal (e.g. `n_samples=1, n_steps=5`) for fast smoke verification.

---

## 5. Verification Checklist

Before finishing any Ribbon Task implementation:
- [ ] `ribbon list tasks` lists the new task.
- [ ] `python -c "import ribbon; print(ribbon.MyTask)"` imports cleanly.
- [ ] Unit tests pass: `pytest tests/test_<task>.py -v`.
- [ ] Smoke dry-run passes: `pytest tests/smoke/test_tasks_smoke.py -k "<TaskName> and dry_run"`.
- [ ] Full smoke dry-run suite has zero regressions: `pytest tests/smoke/test_tasks_smoke.py -k "dry_run"`.
- [ ] GPU live run passes (if GPU is present): `pytest tests/test_<task>.py --run-gpu-smoke -v -s`.
