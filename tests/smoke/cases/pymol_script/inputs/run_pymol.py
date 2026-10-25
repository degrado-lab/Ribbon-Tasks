import sys
from pathlib import Path
import pymol
from pymol import cmd

pymol.finish_launching(['pymol', '-cq'])

out_file = sys.argv[1] if len(sys.argv) > 1 else "output.txt"
out_path = Path(out_file)
out_path.parent.mkdir(parents=True, exist_ok=True)
with open(out_path, "w") as f:
    f.write("PyMol script executed successfully!\n")

cmd.quit()
