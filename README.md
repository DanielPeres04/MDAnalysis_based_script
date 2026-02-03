# MDAnalysis Automated Pipeline

A Python tool for automated Molecular Dynamics analysis. It handles complex detection, ligand distribution, and trajectory analysis (RMSD, Distances, Rgyr) with automatic plotting.

## Features
- **Auto-Detection:** Identifies protein complexes and distributes ligands (e.g., SAM, NAS) automatically.
- **Smart Fitting:** Calculates RMSD by automatically identifying the "Static Anchor" of the protein (Whole Protein minus the Mobile Loop you define).
- **Visualization:** Generates `.csv` data and `.png` plots with customizable Y-axis scaling.

## Setup & Configuration
Open `MD_Analysis_final.py` and adjust the **constants** at the top of the file to match your system.

### 1. Analysis Targets
* **`MOBILE_RESIDUES_KEY`**: The selection string for the **flexible part** of your protein.
  * *Example:* `"resid 120-135"` (The script will automatically use the *rest* of the protein as the static anchor for alignment).
* **`LIGANDS_KEYS`**: List of ligand residue names.
  * *Example:* `["resname SAM", "resname NAS"]`
* **`ATOM_PAIRS`**: List of atom ID tuples (GROMACS 1-based indexing) for distance monitoring.
  * *Example:* `[(450, 500), (1020, 1025)]`
  * *Tip:* Open your `.gro` in VMD to find these numbers.

### 2. Plot Scaling (Y-Limits)
Control the Y-axis limits in the `main()` function to ensure consistent graphs:
* **`PROTEIN_Y_LIM`**: `[0, 15]` (Forces protein RMSD plots to 0-15 Å).
* **`LIGAND_Y_LIM`**: `[0, 10]` (Forces ligand RMSD plots to 0-10 Å).
* **`DISTANCE_YLIM`**: `None` (Auto-scales each distance plot individually).
* **`RGYR_YLIM`**: `None` (Auto-scales to show subtle fluctuations).

## Dependencies
Ensure you have the following Python libraries installed:

```bash
pip install MDAnalysis pandas matplotlib
```

## How to Run
Run the script from the terminal passing your topology and trajectory files:

```bash
python MD_Analysis_final.py system.gro trajectory.xtc
```


