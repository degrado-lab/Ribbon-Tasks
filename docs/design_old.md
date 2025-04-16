
# Task Definitions

## How does Ribbon use Ribbon-Tasks?
The **Ribbon-Tasks** repo contains the `ribbon_tasks` python module, which is imported at runtime by Ribbon. This dynamic import allows for flexibility - you can use the same user-defined `Tasks` when running manually or when queueing a job, as long as **Ribbon** knows where to find the **Ribbon-Tasks** repo. This is set using the `$RIBBON_TASKS_DIR` environment variable, which can be set manually. (See [Customizing your Tasks](customizing_tasks.md)).

## Anatomy of a Task

The `ribbon_tasks` module defines `Tasks` in 3 main components:

1. `tasks.py` - contains the python wrapper code for each `Task`
2. `tasks.json` - defines what command is run for each `Task`, and in which container.
3. `containers.json` - defines the local name and cloud location for each container.

Let's go through these one-by-one.

## 1. `tasks.py`

When the user creates a `Task` in Ribbon, they are instantiating an object defined by a class in `tasks.py`. 

Let's take a look at the Chai-1 Class:
<details>
    <summary>The Chai-1 Task Class</summary>
    <br>
    ```python
    class Chai1(Task):
        def __init__(self, fasta_file, output_dir='.', smiles_string=None, num_ligands=1, device='gpu'):
            """
            Initialize a Chai-1 task.
            
            Args:
                fasta_file (str): The FASTA file containing the protein sequence (no ligand).
                output_dir (str): The directory to save the output files. Default is '.'.
                smiles_string (str, optional): The SMILES string of the ligand. Default is None.
                num_ligands (int): The number of ligands. Default is 1.
                device (str): The device to run the task on. Default is 'gpu'.
            """
            # Initialize the Task class
            super().__init__()

            # This Task name matches the name in the tasks.json file
            self.task_name = "Chai-1"

            # Task-specific variables
            self.fasta_file = fasta_file
            self.smiles_string = smiles_string
            self.output_dir = output_dir
            self.device = device
            self.num_ligands = num_ligands

        def run(self):
            # Make the directory:
            self.output_dir = make_directory(self.output_dir)

            # Run the task
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

As you can see, each `Task` definition contains two functions: the `__init__()` function, and the `run()` function.

The `__init()__` function *must* take in Task arguments, as well as initialize the parent `Task` class and define the Task name using the following code. 

```python
# Initialize the Task class
super().__init__()

# This Task name matches the name in the tasks.json file
self.task_name = "Chai-1"
```

Finally, it must save the user arguments as object variables, which are then referenced when the Task is run:
```python
# Task-specific variables
self.fasta_file = fasta_file
self.smiles_string = smiles_string
...
```


The `run()` function *must* call the `self._run_task()` function, passing in the user-defined object variables. `_run_task()` is what launches the Task by running the correct command in the container (both defined in `tasks.json`). For this reason, when you pass in variables to `_run_task()`, the arguments names are the command arguments defined in `tasks.json`, which are surrounded by curly braces in the command string.

Additionally, `_run_task()` should contain any simple helper code to run a job, including making the correct file structures, renaming output files, etc. However, if heavy python processing is required, we recommend creating a dedicated script and running that as the job itself.

Finally, it's impractical to add every argument for a piece of software as an input to `__init__()`. For this reason, I add `_run_task()` excepts an `extra_args` input, which is a string which is often appended to the end of the command string. This allows the user to add any other software arguments as a simple string. However, this allows for the injection of arbitrary commands - **so use with caution**!

## 2. `tasks.json`

The `tasks.json` file houses a dictionary in JSON format with information on each Task, including the task name, it's container, any environment variables to pass to the job, and the command string to run inside the container.

Let's take a look at the json definition for Chai-1:
<details>
    <summary>The Chai-1 Task Dict</summary>
    <br>
    ```python
    "Chai-1": {
        "name": "Chai-1",
        "description": "Run Chai-1 structure prediction on a FASTA file and SMILES.",
        "container": "Chai-1",
        "command": "python $RIBBON_TASKS_MODULE_DIR/task_scripts/chai-1/predict_structure.py {fasta_file} \" {smiles_string} \" {output_dir} --num_ligands {num_ligands} {extra_args}",
        "requirements": "None",
        "environment_variables": {
            "TRANSFORMERS_OFFLINE": "1"
        }
    },
    ```
</details>

The `container` name points to a container in `containers.json`. 

The `command` string is what gets run inside the container - to pass in user variables, you surround named arguments with brackets. These arguments are digested by the `Task` parent class, and are used as the argument names within the `_run_task()` function. 

So, if your command looks like `"my_command {input1} {input2}"`, you'd start the task using something like `_run_task( input1=self.user_input1, input2=self.user_input2 )`, where the self variables are defined in the `__init__()` function.

Custom Scripts:
- Often, the command string will run a terminal command which is already installed in the container. Sometimes, we instead write custom scripts to do a Task. For example, the Chai-1 task is run using a custom script which processess the input FASTA files into the correct format, before running the folding software using the python API. This "heavy lifting" should be done in a custom script to keep the Task classes lightweight and clear.
- To store the custom scripts in a consistent location, we use a new directory for each **container** under the `task_scripts` directory. Place your script into the directory, and then you can reference it from the command string using something like `"python $RIBBON_TASKS_MODULE_DIR/task_scripts/chai-1/predict_structure.py ..."`
- The `$RIBBON_TASKS_MODULE_DIR` is automatically set within the container, and refers to the location of the `ribbon_tasks` folder.

## 3. `containers.json`
This JSON file contains a dictionary keyed with the name of each container. Each entry in the dictionary is a list that contains the local name of the container (once it's downloaded), and the Dockerhub link (provided by oras, since we'll be downloading using the `apptainer pull` command under the hood).

For example:
```"Chai-1":     ["chai-1_0.5.2.sif", "oras://docker.io/nicholasfreitas/chai-1:0.5.2"],```

To maintain correct version compatibility, it's important to reference each container on Dockerhub with its specific version tag (rather than simply using `latest`), and to set the local container name to similarly have this information.
