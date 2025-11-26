<<<<<<< HEAD
# comp_soft_architecture
MSC Life science 1. Semester computer and software architecture
=======
# Molecular Docking Workflow with Snakemake

A complete workflow for protein-ligand docking using Snakemake, based on the TeachOpenCADD tutorial T015 (https://projects.volkamerlab.org/teachopencadd/talktorials/T015_protein_ligand_docking.html).

## Overview

This pipeline automates the molecular docking process:
1. **Prepare structures** - Extract protein and ligand from PDB
2. **Convert to PDBQT** - Format conversion for docking
3. **Molecular docking** - Run smina (AutoDock Vina fork) with 5 poses per ligand
4. **Post-processing** - Extract docking scores
5. **Split poses** - Separate 5 individual docking poses
6. **Visualization** - Generate 3D protein-ligand complex images using PyMOL for each pose found in the PDB file data

## Requirements

- Python 3.10+
- Conda (Miniforge or Anaconda)

## Setup

### 1. Clone the repository
```bash
git clone <your-repo-url>
cd dock_project
```

### 2. Create the conda environment
```bash
conda env create -f environment.yml
conda activate dock_env
```

This installs all dependencies:
- snakemake (workflow engine)
- openbabel (structure format conversion)
- pymol-open-source (3D visualization)
- rdkit (molecular toolkit)
- smina (docking engine)
- And supporting libraries

### 3. Add your input PDB files

Place PDB files in `data/raw/`:
```
data/raw/2ito.pdb
data/raw/1A2C.pdb (example)
...
```
### 4. Add your input PDB files to config.yaml file

Update `config.yaml` with the PDB file names (without `.pdb` extension). Only include files that are actually in `data/raw/`:
```yaml
pdb_files:
  - "2ito"
  - "1A2C"
```

## Running the Workflow

### 1. Dry-run (optional)
```bash
snakemake -n
```

### 2. Execute the workflow
```bash
snakemake -j 1 #(or other number can be used)
```
(Use `-j x` to run x "jobs" in parallel)

### Force re-run all steps
```bash
snakemake -F -j 1
```

## Output Structure with two example PDB Files:

```
data/
├── processed/
│   ├── 2ito/
│   │   ├── 2ito_protein.pdb
│   │   ├── 2ito_ligand.pdb
│   │   ├── 2ito_protein.pdbqt
│   │   └── 2ito_ligand.pdbqt
│   └── 1A2C/
│       ├── 1A2C_protein.pdb
│       ├── 1A2C_ligand.pdb
│       ├── 1A2C_protein.pdbqt
│       └── 1A2C_ligand.pdbqt
└── docking/
    ├── 2ito/
    │   ├── 2ito_docked.sdf          # All poses that were in PDB File
    │   ├── 2ito_scores.csv          # Affinity scores
    │   ├── poses/
    │   │   ├── pose_00.sdf
    │   │   ├── pose_01.sdf
    │   │   └── ...
    │   └── visualize/
    │       ├── pose_00.png
    │       ├── pose_01.png
    │       └── ...
    └── 1A2C/
        ├── 1A2C_docked.sdf
        ├── 1A2C_scores.csv
        ├── poses/
        │   ├── pose_00.sdf
        │   ├── pose_01.sdf
        │   └── ...
        └── visualize/
            ├── pose_00.png
            ├── pose_01.png
            └── ...
```

## Workflow Visualization

A DAG (Directed Acyclic Graph) visualization is included as `dag_protein_docking.svg`. This shows how the 6 rules are connected and their dependencies.

To regenerate the DAG after modifying the Snakefile:
```bash
snakemake --dag --snakefile Snakefile | dot -Tsvg > dag_protein_docking.svg
```
(Requires graphviz: `conda install -c conda-forge graphviz`)

## File Structure

- **Snakefile** - Workflow definition with all 6 rules
- **environment.yml** - Conda dependencies specification
- **config.yaml** - Project configuration (paths, PDB files)
- **scripts/** - Python scripts for each workflow step
  - `prepare_structures.py` - Extract protein/ligand from PDB
  - `convert_to_pdbqt.py` - Convert to PDBQT format
  - `run_analysis.py` - Extract docking scores
  - `visualize_pose.py` - Generate 3D visualizations


## Citation

Based on TeachOpenCADD T015: Protein-Ligand Docking
https://projects.volkamerlab.org/teachopencadd/talktorials/T015_protein_ligand_docking.html

PDB Files used in creation:
https://www.rcsb.org/structure/2ITO
https://www.rcsb.org/structure/2XNI

## Contact

FHNW

## Known limitations

- **Unsupported atom types (e.g. Boron)**: The pipeline converts ligand and protein files to the PDBQT format and runs `smina` for docking. The AutoDock/Vina-style PDBQT format (used by `smina`) supports a limited set of atom types. Ligands that contain some elements (for example, boron `B`) may fail during PDBQT parsing and cause the workflow to stop with a parse error.

  Detection (quick): run this to list element symbols found in a ligand PDB:

  ```bash
  awk '{print substr($0,77,2)}' data/processed/<PDB>/<PDB>_ligand.pdb | sort | uniq -c
  ```
>>>>>>> 59d5283 (Initial project: Snakemake docking pipeline)
