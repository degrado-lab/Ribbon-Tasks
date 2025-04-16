# Customizing Tasks

You may want to add new types of Tasks to your install, or contribute to the Ribbon-Tasks repo for other users. Ribbon is built to make this process easy.

## Where are the Task definitions stored?

By default, the Ribbon-Tasks repo is downloaded to `$USER/.ribbon/ribbon_tasks`. The `.ribbon` directory is also where you'll find the cache of submitted ribbon tasks, as well as downloaded containers.

## I want to install it somwhere else!

If you want to create your own copy of the repo, you can download or clone the repo manually from [GitHub](https://github.com/degrado-lab/Ribbon-Tasks/releases).

Then, you must set the `RIBBON_TASKS_DIR` environment variable to this new location using `export RIBBON_TASKS_DIR=[your new location]`. For this change to be saved across multiple terminal sessions, I recommend adding this line to your `~/.bashrc` file - then, the change will be loaded each time you start a terminal.

## The Tasks were updated in GitHub! How do I update them?

For now, you'll have to download the newest [release](https://github.com/degrado-lab/Ribbon-Tasks/releases), and use the above process to point Ribbon towards the new task Repo. In the future, I would like to add a Task updating feature.

## I want to make my own Tasks!

All Tasks are defined in a simple python file and two JSON files. Please take a look at `Design Notes` for information on how the Tasks are defined`.