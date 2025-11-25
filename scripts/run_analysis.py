"""
Parse SDF file to extract docking scores and write to CSV.

snakemake input: docked SDF file
snakemake output: scores CSV file
"""

import os

in_sdf = snakemake.input[0]
out_csv = snakemake.output[0]

os.makedirs(os.path.dirname(out_csv), exist_ok=True)

scores = []
pose_id = 0
capture_score = False

with open(in_sdf) as fh:
    for line in fh:
        line = line.strip()
        if line == "$$$$":  # End of molecule record
            pose_id += 1
            capture_score = False
        elif line == ">  <minimizedAffinity>":
            capture_score = True
        elif capture_score and line:
            # This line should be the score
            try:
                score = float(line)
                scores.append((pose_id, score))
            except ValueError:
                pass
            capture_score = False

with open(out_csv, "w") as fo:
    fo.write("pose_id,score\n")
    for pid, sc in scores:
        fo.write(f"{pid},{sc}\n")
