import sys
from pathlib import Path
from Bio.Seq import Seq

out_file = sys.argv[1] if len(sys.argv) > 1 else "output.txt"
out_path = Path(out_file)
out_path.parent.mkdir(parents=True, exist_ok=True)
seq = Seq("ACGT")
with open(out_path, "w") as f:
    f.write(f"BioPython script executed: rev_comp={seq.reverse_complement()}\n")
