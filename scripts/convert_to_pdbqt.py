"""
Convert protein and ligand PDB files to PDBQT format using Open Babel.

snakemake input: protein.pdb and ligand.pdb
snakemake output: protein.pdbqt and ligand.pdbqt
"""

import os
import subprocess

in_prot = snakemake.input[0]
in_lig = snakemake.input[1]
out_prot_qt = snakemake.output[0]
out_lig_qt = snakemake.output[1]

os.makedirs(os.path.dirname(out_prot_qt), exist_ok=True)

subprocess.check_call(
    ["obabel", "-i", "pdb", in_prot, "-o", "pdbqt", "-O", out_prot_qt, "--addhyd"]
)

subprocess.check_call(
    ["obabel", "-i", "pdb", in_lig, "-o", "pdbqt", "-O", out_lig_qt, "--addhyd"]
)
