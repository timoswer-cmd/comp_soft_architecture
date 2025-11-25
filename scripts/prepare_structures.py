"""
Split PDB into protein.pdb (ATOM) and ligand.pdb (largest HETATM group).

snakemake input: raw PDB file
snakemake output: protein.pdb and ligand.pdb
"""

import os
from collections import defaultdict

in_pdb = snakemake.input[0]
out_protein = snakemake.output[0]
out_ligand = snakemake.output[1]

os.makedirs(os.path.dirname(out_protein), exist_ok=True)

atom_lines = []
hetatm_lines = []

with open(in_pdb) as fh:
    for line in fh:
        if line.startswith("ATOM"):
            atom_lines.append(line)
        elif line.startswith("HETATM"):
            hetatm_lines.append(line)

with open(out_protein, "w") as fo:
    fo.write("".join(atom_lines))

groups = defaultdict(list)
for ln in hetatm_lines:
    key = (ln[17:20].strip(), ln[21].strip(), ln[22:26].strip())
    groups[key].append(ln)

if not groups:
    with open(out_ligand, "w") as fo:
        fo.write("")
else:
    best = max(groups.items(), key=lambda kv: len(kv[1]))
    with open(out_ligand, "w") as fo:
        fo.write("".join(best[1]))
