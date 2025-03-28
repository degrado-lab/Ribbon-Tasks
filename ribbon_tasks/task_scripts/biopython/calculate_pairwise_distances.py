#!/usr/bin/env python3

import argparse
import json
import csv
import sys
from itertools import repeat

from Bio import PDB

def parse_args():
    parser = argparse.ArgumentParser(
        description="Match atoms from two groups via minimal-distance bipartite matching (Hungarian algorithm)."
    )
    parser.add_argument("structure_file",
                        help="Path to a PDB or mmCIF file.")
    parser.add_argument("--groupA",
                        nargs="+",
                        required=True,
                        help="Atoms for group A (format 'Chain:Residue:Atom'), space-separated.")
    parser.add_argument("--groupB",
                        nargs="+",
                        required=True,
                        help="Atoms for group B (format 'Chain:Residue:Atom'), space-separated.")
    parser.add_argument("--csv_output",
                        default=None,
                        help="Path to CSV output file (optional).")
    parser.add_argument("--json_output",
                        default=None,
                        help="Path to JSON output file (optional).")
    parser.add_argument("--average_output",
                        default=None,
                        help="Path to a text file containing the average distance for matched pairs (and no other information).")
    return parser.parse_args()

def parse_atom_spec(atom_spec):
    """
    Expects 'Chain:Residue:AtomName', e.g. 'B:99:CG'.
    Returns (chain_id, residue_id, atom_name).
    """
    parts = atom_spec.split(":")
    if len(parts) != 3:
        raise ValueError(f"Invalid atom spec: {atom_spec}. Use 'Chain:Residue:AtomName'.")
    chain_id, residue_id, atom_name = parts
    return chain_id, residue_id, atom_name

def get_residue(chain, res_id):
    """
    Retrieve a residue from the chain (standard or hetero).
    """
    try:
        residue = chain[(" ", int(res_id), " ")]
    except KeyError:
        # If not found, look among hetero-residues
        hetero_residues = [r for r in chain.get_list() if r.id[0] != " "]
        for r in hetero_residues:
            if r.id[1] == int(res_id):
                return r
        raise KeyError(f"Residue {res_id} not found in chain {chain.id}.")
    return residue

def get_atom_distance(structure, atomA_spec, atomB_spec):
    """
    Given two atom specs and a structure, return the distance.
    """

    model = structure[0]  # use first model
    chainA, resA, nameA = parse_atom_spec(atomA_spec)
    chainB, resB, nameB = parse_atom_spec(atomB_spec)

    if chainA == '~' or chainB == '~':
        return float("inf")  # We pad with a wildcard ~. Set these at arbitrarily large distance.

    chainA_obj = model[chainA]
    chainB_obj = model[chainB]

    residueA = get_residue(chainA_obj, resA)
    residueB = get_residue(chainB_obj, resB)

    try:
        atomA_obj = residueA[nameA]
        atomB_obj = residueB[nameB]
    except KeyError as e:
        print(f"Atom not found: {atomA_spec} or {atomB_spec}")
        print(f"Atoms in residue {resA}: {[a.name for a in residueA.get_atoms()]}")
        print(f"Atoms in residue {resB}: {[a.name for a in residueB.get_atoms()]}")
        raise e

    return atomA_obj - atomB_obj

def build_distance_matrix(structure, groupA, groupB):
    """
    Returns a list-of-lists (matrix) with distances.
    matrix[i][j] = distance between groupA[i] and groupB[j].
    If an atom cannot be found, we set distance to a large fallback.
    """
    matrix = []
    for A_spec in groupA:
        row = []
        for B_spec in groupB:
            dist = get_atom_distance(structure, A_spec, B_spec)
            row.append(dist)
        matrix.append(row)
    return matrix

######################################################################
# Hungarian (Kuhn–Munkres) algorithm in pure Python (no SciPy needed) #
######################################################################
def hungarian_algorithm(cost_matrix, match_length=None):
    """
    An implementation of the Hungarian (Kuhn–Munkres) algorithm in pure Python.

    cost_matrix: list of lists, shape NxN (square)
    returns: (row_indices, col_indices)
        row_indices[i] is matched with col_indices[i].
    """
    n = len(cost_matrix)
    # For safety, ensure it's square
    for row in cost_matrix:
        if len(row) != n:
            raise ValueError("Cost matrix must be NxN.")

    # Step 1: Subtract row minima
    for i in range(n):
        row_min = min(cost_matrix[i])
        for j in range(n):
            cost_matrix[i][j] -= row_min

    # Step 2: Subtract column minima
    for j in range(n):
        col_min = min(cost_matrix[i][j] for i in range(n))
        for i in range(n):
            cost_matrix[i][j] -= col_min

    # The rest is iterative: we try to cover all zeroes with a minimum number of lines,
    # then adjust the matrix if not enough lines exist. We'll keep track with:
    # - 'match_row' and 'match_col': which row/col is matched with which col/row.
    #   If match_row[i] = j, means row i is matched with col j.
    #   If match_col[j] = i, means col j is matched with row i.

    # Initialize match arrays
    match_row = [-1] * n  # match for each row
    match_col = [-1] * n  # match for each column

    # Markers used for the BFS/DFS steps in the Hungarian algorithm
    def find_match(row, row_seen, col_seen):
        """
        Attempt to find a match for 'row' using DFS or BFS (here DFS).
        row_seen, col_seen keep track of visited rows/columns in this search.
        """
        for col in range(n):
            # If cost is zero and not yet visited in the DFS:
            if cost_matrix[row][col] == 0 and not col_seen[col]:
                col_seen[col] = True
                # If col is not matched OR we can re-match the matched row
                if match_col[col] == -1 or find_match(match_col[col], row_seen, col_seen):
                    match_row[row] = col
                    match_col[col] = row
                    return True
        return False

    # We'll try to find a match for each row
    # (This is effectively a "maximum bipartite matching" on the zero edges.)
    for row in range(n):
        row_seen = [False] * n
        col_seen = [False] * n
        find_match(row, row_seen, col_seen)

    # Count how many matched pairs we have
    matches = sum(1 for x in match_row if x != -1)
    if matches == n:
        # We have a perfect matching
        row_ind = []
        col_ind = []
        for r, c in enumerate(match_row):
            row_ind.append(r)
            col_ind.append(c)
        return row_ind, col_ind

    # If not all matched, we do the Hungarian updates (cover lines, modify matrix) repeatedly.
    # Implementation detail: We'll do repeated expansions until we get a full matching.

    def min_line_cover():
        """
        The standard approach:
        1) Start with all unmatched rows as 'marked'.
        2) For each marked row, mark any col with zero in that row.
        3) For each marked col, mark the row that is matched with that col.
        4) Repeat until no new marks.
        Then lines to cover zeros: unmarked rows + marked cols.
        """
        marked_rows = [False] * n
        marked_cols = [False] * n

        # Step 1: Mark all unmatched rows
        unmatched_rows = [r for r in range(n) if match_row[r] == -1]
        stack = unmatched_rows[:]
        for r in unmatched_rows:
            marked_rows[r] = True

        # BFS
        while stack:
            r = stack.pop()
            # For each zero in row r, mark that column
            for c in range(n):
                if cost_matrix[r][c] == 0 and not marked_cols[c]:
                    marked_cols[c] = True
                    # Also mark the row matched with col c (if it exists)
                    row_match = match_col[c]
                    if row_match != -1 and not marked_rows[row_match]:
                        marked_rows[row_match] = True
                        stack.append(row_match)

        return marked_rows, marked_cols

    if match_length is None:
        match_length = n

    attempts = 0
    while matches < match_length:
        attempts += 1
        if attempts > 100:
            raise ValueError("Failed to find a full matching.")
        marked_rows, marked_cols = min_line_cover()

        # Count how many lines needed:
        # lines = (rows that are NOT marked) + (cols that ARE marked)
        # i.e. we "cover" all unmarked rows with a horizontal line
        # and all marked cols with a vertical line.
        count_lines = sum(not mr for mr in marked_rows) + sum(mc for mc in marked_cols)

        if count_lines == n:
            # We have a minimal cover, so we can build a perfect matching
            break

        # Otherwise, find the smallest un-covered value and subtract it from
        # all uncovered elements, then add it to elements covered by two lines.
        min_uncovered_val = float("inf")
        for r in range(n):
            for c in range(n):
                # if row r is not marked and col c is not marked
                if marked_rows[r] == False and marked_cols[c] == False:
                    if cost_matrix[r][c] < min_uncovered_val:
                        min_uncovered_val = cost_matrix[r][c]

        # Subtract from all uncovered elements, add to covered by two lines
        for r in range(n):
            for c in range(n):
                if marked_rows[r] == False and marked_cols[c] == False:
                    cost_matrix[r][c] -= min_uncovered_val
                elif marked_rows[r] == True and marked_cols[c] == True:
                    cost_matrix[r][c] += min_uncovered_val

        # Now try finding a matching again with the updated matrix
        match_row = [-1] * n
        match_col = [-1] * n
        matches = 0
        for row in range(n):
            row_seen = [False] * n
            col_seen = [False] * n
            if find_match(row, row_seen, col_seen):
                matches += 1

    # By now, we should have a complete matching
    row_ind, col_ind = [], []
    for r, c in enumerate(match_row):
        row_ind.append(r)
        col_ind.append(c)
    return row_ind, col_ind

######################################################################

def main():
    args = parse_args()

    # Parse structure
    if args.structure_file.endswith(".cif"):
        parser = PDB.MMCIFParser(QUIET=True)
    else:
        parser = PDB.PDBParser(QUIET=True)

    structure = parser.get_structure("my_structure", args.structure_file)

    groupA = args.groupA
    groupB = args.groupB

    # For the Hungarian algorithm, groups must be same length if you need perfect 1-to-1 matching.
    # If they aren't, you can do an NxN matrix by padding the smaller group with "dummy" atoms
    # or handle partial matching. Here we'll assume lenA == lenB for a perfect match.
    expected_match_length = min(len(groupA), len(groupB))
    user_group_lengths = (len(groupA), len(groupB))
    if len(groupA) != len(groupB):
        print("Warning: For perfect bipartite matching, groupA and groupB should have the same size.", file=sys.stderr)
        # We'll proceed anyway by forming an NxN matrix with N = max(lenA, lenB).
        N = max(len(groupA), len(groupB))
        # We'll pad whichever is shorter with repeats of a wildcard (that will yield infinite cost in the distance function).
        # That ensures we have an NxN matrix.
        if len(groupA) < N:
            groupA += ['~:~:~'] * (N - len(groupA))
        if len(groupB) < N:
            groupB += ['~:~:~'] * (N - len(groupB))
    else:
        N = len(groupA)

    # 1) Build the actual distance matrix
    original_matrix = build_distance_matrix(structure, groupA, groupB)

    # 2) Deep-copy that matrix for the Hungarian routine
    import copy
    cost_matrix = copy.deepcopy(original_matrix)

    # 3) Hungarian match on the copy
    row_ind, col_ind = hungarian_algorithm(cost_matrix, match_length=expected_match_length)

    results = []
    for i in range(len(row_ind)):
        r = row_ind[i]
        c = col_ind[i]
        # Now fetch the *actual* distance from 'original_matrix'
        dist = original_matrix[r][c]
        results.append({
            "A_index": r,
            "A_spec": groupA[r],
            "B_index": c,
            "B_spec": groupB[c],
            "distance": float(dist) # Convert from float32 to default float
        })

    #If we have fewer A than B (or vice-versa), throw out the extra results:
    results_truncated = []
    for result in results:
        if 0 <= result["A_index"] < user_group_lengths[0] and 0 <= result["B_index"] < user_group_lengths[1]:
            results_truncated.append(result)
    results = results_truncated

    # 4) Print a human-readable table
    print("\nMatched Pairs (Hungarian Minimal-Distance Assignment):\n")
    print(" A_idx |         Atom A        | B_idx |         Atom B        | Distance")
    print("-------+-----------------------+-------+-----------------------+---------")
    for row in results:
        dist_str = f"{row['distance']:.2f}" if row["distance"] is not None else "N/A"
        print(
            f"{row['A_index']:>5}  | "
            f"{row['A_spec']:<21} | "
            f"{row['B_index']:>5}  | "
            f"{row['B_spec']:<21} | "
            f"{dist_str}"
        )

    # 5) Optionally save CSV
    if args.csv_output:
        with open(args.csv_output, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["A_index", "AtomA", "B_index", "AtomB", "Distance"])
            for row in results:
                dist_str = f"{row['distance']:.2f}" if row["distance"] else ""
                writer.writerow([row["A_index"], row["A_spec"],
                                 row["B_index"], row["B_spec"],
                                 dist_str])
        print(f"\nWrote CSV output to {args.csv_output}")

    # 6) Optionally save JSON
    if args.json_output:
        # Distances might be None if something wasn't found
        with open(args.json_output, "w") as jf:
            json.dump(results, jf, indent=2)
        print(f"Wrote JSON output to {args.json_output}")
    
    # 7) If user wants average distance output, compute & write to a text file
    if args.average_output:
        # Filter out None distances (or 1e9 if you prefer to skip those)
        valid_distances = [r["distance"] for r in results if r["distance"] is not None]
        if len(valid_distances) > 0:
            avg_dist = sum(valid_distances) / len(valid_distances)
        else:
            # No valid distances -> you can decide to output 0.00, or skip
            avg_dist = 0.0

        with open(args.average_output, "w") as f:
            # Write ONLY the numeric distance, with no extra text
            f.write(f"{avg_dist:.2f}\n")
        print(f"Wrote average distance to {args.average_output}")
    

if __name__ == "__main__":
    main()
