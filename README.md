# MDAnalysis Automated Pipeline

A Python script to automate the analysis of molecular dynamics simulations (GROMACS format). It automatically detects protein complexes, distributes ligands, and calculates RMSD and Atom Distances for specific binding sites.

## Features
- **Auto-Detection:** Automatically separates protein chains into complexes (e.g., Dimer, Tetramer) without needing a `.tpr` topology.
- **Ligand Sorting:** Cyclically distributes ligands (e.g., SAM, NAS) to their respective protein complexes.
- **Analysis Modules:**
  - **RMSD (Protein-Fitted):** Tracks ligand movement relative to the protein backbone.
  - **RMSD (Ligand-Fitted):** Tracks the internal flexibility of the ligand itself.
  - **Atom Distances:** Monitors catalytic distances over time.
  - **Radius of Gyration:** Measures ligand compactness.

## Dependencies
Ensure you have the following Python libraries installed:

```bash
pip install MDAnalysis pandas matplotlib
```
