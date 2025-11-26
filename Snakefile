# Snakefile
"""Snakemake workflow for molecular docking using smina.

Steps:
1. Prepare structures: Split PDB into protein and ligand files.
2. Convert to PDBQT format.
3. Run docking with smina.
4. Post-process results to extract scores.
5. Split docking poses into individual SDF files.
6. Visualize poses using PyMOL or RDKit (if PyMOL isn't working).
"""
# import statements
import os
configfile: "config.yaml"

# Define directories
RAW_DIR = config["raw_dir"]
PROC_DIR = config["proc_dir"]
DOCK_DIR = config["dock_dir"]
SCRIPTS = config["scripts_dir"]

PDB_FILES = config["pdb_files"]


rule all:
    input:
        expand(os.path.join(DOCK_DIR, "{pdb}", "{pdb}_docked.sdf"), pdb=PDB_FILES),
        expand(os.path.join(DOCK_DIR, "{pdb}", "{pdb}_scores.csv"), pdb=PDB_FILES),
        expand(os.path.join(DOCK_DIR, "{pdb}", "poses", ".done"), pdb=PDB_FILES),
        expand(os.path.join(DOCK_DIR, "{pdb}", "visualize", "{pose}.png"), pdb=PDB_FILES, pose=[f"pose_{i:02d}" for i in range(5)])


rule prepare_structures:
    input:
        os.path.join(RAW_DIR, "{pdb}.pdb")
    output:
        protein=os.path.join(PROC_DIR, "{pdb}", "{pdb}_protein.pdb"),
        ligand=os.path.join(PROC_DIR, "{pdb}", "{pdb}_ligand.pdb")
    conda:
        "environment.yml"
    log:
        os.path.join(DOCK_DIR, "{pdb}", "logs", "prepare_structures.log")
    benchmark:
        os.path.join(DOCK_DIR, "{pdb}", "benchmarks", "prepare_structures.bench.txt")
    script:
        os.path.join(SCRIPTS, "prepare_structures.py")


rule convert_to_pdbqt:
    input:
        protein=os.path.join(PROC_DIR, "{pdb}", "{pdb}_protein.pdb"),
        ligand=os.path.join(PROC_DIR, "{pdb}", "{pdb}_ligand.pdb")
    output:
        protein_qt=os.path.join(PROC_DIR, "{pdb}", "{pdb}_protein.pdbqt"),
        ligand_qt=os.path.join(PROC_DIR, "{pdb}", "{pdb}_ligand.pdbqt")
    conda:
        "environment.yml"
    log:
        os.path.join(DOCK_DIR, "{pdb}", "logs", "convert_to_pdbqt.log")
    benchmark:
        os.path.join(DOCK_DIR, "{pdb}", "benchmarks", "convert_to_pdbqt.bench.txt")
    script:
        os.path.join(SCRIPTS, "convert_to_pdbqt.py")


rule run_docking:
    input:
        protein_qt=os.path.join(PROC_DIR, "{pdb}", "{pdb}_protein.pdbqt"),
        ligand_qt=os.path.join(PROC_DIR, "{pdb}", "{pdb}_ligand.pdbqt")
    output:
        os.path.join(DOCK_DIR, "{pdb}", "{pdb}_docked.sdf")
    params:
        outdir=lambda wildcards: os.path.join(DOCK_DIR, wildcards.pdb),
        num_modes=5,
        autobox_add=4
    conda:
        "environment.yml"
    log:
        os.path.join(DOCK_DIR, "{pdb}", "logs", "run_docking.log")
    benchmark:
        os.path.join(DOCK_DIR, "{pdb}", "benchmarks", "run_docking.bench.txt")
    shell:
        "mkdir -p {params.outdir} {params.outdir}/logs {params.outdir}/benchmarks && "
        "( /usr/bin/time -v -o {params.outdir}/logs/run_docking.time.txt -- bash -lc \"smina -r {input.protein_qt} -l {input.ligand_qt} --autobox_ligand {input.ligand_qt} --autobox_add {params.autobox_add} --num_modes {params.num_modes} -o {output}\" ) > {log} 2>&1"


rule postprocess:
    input:
        os.path.join(DOCK_DIR, "{pdb}", "{pdb}_docked.sdf")
    output:
        os.path.join(DOCK_DIR, "{pdb}", "{pdb}_scores.csv")
    conda:
        "environment.yml"
    log:
        os.path.join(DOCK_DIR, "{pdb}", "logs", "postprocess.log")
    benchmark:
        os.path.join(DOCK_DIR, "{pdb}", "benchmarks", "postprocess.bench.txt")
    script:
        os.path.join(SCRIPTS, "run_analysis.py")


rule split_docking_poses:
    input:
        os.path.join(DOCK_DIR, "{pdb}", "{pdb}_docked.sdf")
    output:
        os.path.join(DOCK_DIR, "{pdb}", "poses", "pose_00.sdf"),
        os.path.join(DOCK_DIR, "{pdb}", "poses", "pose_01.sdf"),
        os.path.join(DOCK_DIR, "{pdb}", "poses", "pose_02.sdf"),
        os.path.join(DOCK_DIR, "{pdb}", "poses", "pose_03.sdf"),
        os.path.join(DOCK_DIR, "{pdb}", "poses", "pose_04.sdf"),
        os.path.join(DOCK_DIR, "{pdb}", "poses", ".done")
    conda:
        "environment.yml"
    log:
        os.path.join(DOCK_DIR, "{pdb}", "poses", "split_docking_poses.log")
    benchmark:
        os.path.join(DOCK_DIR, "{pdb}", "benchmarks", "split_docking_poses.bench.txt")
    params:
        outdir=lambda wildcards: os.path.join(DOCK_DIR, wildcards.pdb, "poses")
    shell:
        "mkdir -p {params.outdir} {params.outdir}/../benchmarks {params.outdir}/../logs && "
        "( /usr/bin/time -v -o {params.outdir}/../logs/split_docking_poses.time.txt -- bash -lc \"python3 - << 'SPLIT_SCRIPT'\nimport os\nfrom openbabel import pybel\n\nsdf_file = \"{input}\"\nout_dir = \"{params.outdir}\"\nmolecules = pybel.readfile(\"sdf\", sdf_file)\nfor i, mol in enumerate(molecules):\n    out_path = os.path.join(out_dir, f\"pose_{{i:02d}}.sdf\")\n    mol.write(\"sdf\", out_path, overwrite=True)\n\ndone_path = os.path.join(out_dir, \".done\")\nwith open(done_path, \"w\") as f:\n    f.write(\"done\")\nSPLIT_SCRIPT\" ) > {log} 2>&1"


rule visualize_poses:
    input:
        split_done=os.path.join(DOCK_DIR, "{pdb}", "poses", ".done"),
        pose_sdf=os.path.join(DOCK_DIR, "{pdb}", "poses", "{pose}.sdf"),
        protein_pdb=os.path.join(PROC_DIR, "{pdb}", "{pdb}_protein.pdb")
    output:
        os.path.join(DOCK_DIR, "{pdb}", "visualize", "{pose}.png")
    conda:
        "environment.yml"
    log:
        os.path.join(DOCK_DIR, "{pdb}", "logs", "visualize_poses.{pose}.log")
    benchmark:
        os.path.join(DOCK_DIR, "{pdb}", "benchmarks", "visualize_poses.{pose}.bench.txt")
    script:
        os.path.join(SCRIPTS, "visualize_pose.py")
