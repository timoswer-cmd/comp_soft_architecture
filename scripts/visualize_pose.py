"""
Visualize a protein-ligand pose using PyMOL if available, otherwise RDKit 2D rendering.

snakemake input: protein_pdb (PDB file), pose_sdf (SDF file)
snakemake output: PNG image file
"""

import os
import subprocess
import sys
import tempfile

# Resolve output path (support list/tuple or single string)
out = (
    snakemake.output[0]
    if isinstance(snakemake.output, (list, tuple))
    else str(snakemake.output)
)
os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

# Expect inputs named `protein_pdb` and `pose_sdf` from the Snakefile
protein_path = (
    snakemake.input.protein_pdb if hasattr(snakemake.input, "protein_pdb") else None
)
pose_path = snakemake.input.pose_sdf if hasattr(snakemake.input, "pose_sdf") else None

if protein_path is None:
    raise ValueError("visualize_pose.py requires a protein_pdb input")
if pose_path is None:
    raise ValueError("visualize_pose.py requires a pose_sdf input")

# First, try to use PyMOL for 3D visualization
pymol_success = False
try:
    # Check if pymol is available
    result = subprocess.run(
        ["pymol", "--version"], capture_output=True, text=True, timeout=5
    )
    pymol_available = result.returncode == 0

    if pymol_available:
        print(f"PyMOL detected, attempting 3D visualization", file=sys.stderr)

        # Create a PyMOL script for visualization
        with tempfile.NamedTemporaryFile(mode="w", suffix=".pml", delete=False) as f:
            pml_script = f.name
            f.write(
                f"""
# Load structures
load {protein_path}, protein
load {pose_path}, ligand

# Set representations
hide everything
show cartoon, protein
color spectrum, protein
show sticks, ligand
color cyan, ligand

# Zoom and orient
zoom
set bg_rgb, white

# Save image
png {out}, width=1200, height=800, dpi=100
quit
"""
            )

        # Run PyMOL with the script
        pymol_cmd = ["pymol", "-c", pml_script]
        result = subprocess.run(pymol_cmd, capture_output=True, text=True, timeout=30)

        if (
            result.returncode == 0
            and os.path.exists(out)
            and os.path.getsize(out) > 1000
        ):
            pymol_success = True
            print(f"3D visualization with PyMOL saved to: {out}", file=sys.stderr)
        else:
            print(f"PyMOL rendering failed, falling back to RDKit", file=sys.stderr)

        # Clean up temp script
        try:
            os.remove(pml_script)
        except:
            pass

except Exception as e:
    print(f"PyMOL attempt failed: {e}", file=sys.stderr)

# If PyMOL didn't work, use RDKit fallback
if not pymol_success:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import FancyBboxPatch
        from rdkit import Chem
        from rdkit.Chem import AllChem, Crippen, Descriptors, Draw

        print(f"Using RDKit 2D fallback visualization", file=sys.stderr)

        # Load ligand from SDF
        suppl = Chem.SDMolSupplier(pose_path, removeHs=False)
        if suppl is None or len(suppl) == 0:
            raise ValueError(f"Could not read SDF file: {pose_path}")
        mol = suppl[0]
        if mol is None:
            raise ValueError(f"First molecule in {pose_path} is invalid")

        # Generate 2D coordinates for ligand
        AllChem.Compute2DCoords(mol)

        # Draw the ligand structure
        img_ligand = Draw.MolToImage(mol, size=(400, 400))

        # Load protein info
        protein_atoms = []
        try:
            with open(protein_path, "r") as f:
                for line in f:
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        try:
                            protein_atoms.append(line[12:16].strip())
                        except:
                            pass
        except Exception as e:
            print(f"Warning: Could not parse protein: {e}", file=sys.stderr)

        # Create figure with matplotlib
        fig = plt.figure(figsize=(12, 8), dpi=100)
        ax = fig.add_subplot(111)

        # Display title
        ax.text(
            0.5,
            0.95,
            "Protein-Ligand Complex (2D Fallback)",
            ha="center",
            va="top",
            fontsize=14,
            fontweight="bold",
            transform=ax.transAxes,
        )

        # Display protein info
        ax.text(
            0.05,
            0.88,
            f"Protein: {os.path.basename(protein_path)}",
            ha="left",
            va="top",
            fontsize=11,
            transform=ax.transAxes,
            family="monospace",
        )
        ax.text(
            0.05,
            0.82,
            f"Atoms: {len(protein_atoms)}",
            ha="left",
            va="top",
            fontsize=10,
            transform=ax.transAxes,
        )

        # Display ligand info
        ax.text(
            0.05,
            0.74,
            f"Ligand: {os.path.basename(pose_path)}",
            ha="left",
            va="top",
            fontsize=11,
            transform=ax.transAxes,
            family="monospace",
        )
        ax.text(
            0.05,
            0.68,
            f"Atoms: {mol.GetNumAtoms()}, Bonds: {mol.GetNumBonds()}",
            ha="left",
            va="top",
            fontsize=10,
            transform=ax.transAxes,
        )

        # Molecular properties
        try:
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            ax.text(
                0.05,
                0.60,
                f"Formula: {formula}",
                ha="left",
                va="top",
                fontsize=10,
                transform=ax.transAxes,
            )
            ax.text(
                0.05,
                0.54,
                f"MW: {mw:.2f} g/mol, LogP: {logp:.2f}",
                ha="left",
                va="top",
                fontsize=10,
                transform=ax.transAxes,
            )
        except:
            pass

        # Display the 2D ligand structure
        ax_img = fig.add_axes([0.55, 0.35, 0.4, 0.45])
        ax_img.imshow(img_ligand)
        ax_img.axis("off")
        ax_img.set_title("Ligand 2D Structure", fontsize=11, fontweight="bold")

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

        plt.savefig(out, dpi=100, bbox_inches="tight", facecolor="white")
        plt.close()

        print(f"2D visualization with RDKit saved to: {out}", file=sys.stderr)

    except Exception as e:
        print(f"Error during visualization: {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()

        # Create error image
        from PIL import Image, ImageDraw

        img = Image.new("RGB", (800, 600), color="white")
        draw = ImageDraw.Draw(img)
        draw.text((10, 10), f"Rendering Error: {str(e)[:80]}", fill="red")
        img.save(out)
