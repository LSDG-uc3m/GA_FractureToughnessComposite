# GA_FractureToughnessComposite
===============================

## Purpose
----------
This repository provides a Python script that calibrates (G, Gm, n) parameters of a cohesive law
using a Genetic Algorithm (GA) and Abaqus CAE simulations (noGUI).

The typical use case is matching an experimental Force–Displacement curve with a numerical curve
produced by Abaqus, by iteratively updating cohesive parameters.

Intended audience:
- Mid-level engineering students needing a clear workflow
- University researchers in:
  * Continuum Mechanics
  * Structural Mechanics
  * Materials Science / Damage mechanics
  * Cohesive-zone modeling (delamination, fracture, etc.)

## What the script does
-----------------------
For a specific case defined by `(MODE, MATERIAL, LAMINATE)`, the script executes the following pipeline:

1.  **Load Experimental Data** Imports the reference curve from:  
    `./Experimental/Exp_{MODE}_{MATERIAL}_{LAMINATE}.txt`

2.  **Population Initialization** Initializes a population of cohesive parameters $(G, G_m, n)$ within the user-defined bounds.

3.  **Individual Evaluation (Loop)** For each individual in the generation:
    * **Data Injection:** Writes $(G, G_m, n)$ into the corresponding row of `Mechanical_properties.csv` using the **READ_FLAG** protocol.
    * **Simulation:** Runs Abaqus CAE in `noGUI` mode using `Ct_Shell.py`.
    * **Data Extraction:** `Ct_Shell.py` must output a numerical curve file to the root:  
        `NUMERICAL_{MODE}_{MATERIAL}_{LAMINATE}*.txt`  
        *(The script automatically selects the newest matching file if suffixes exist).*
    * **Error Calculation:** Computes the fit error: $raw\_err = 1 - R^2$, evaluated over a common displacement window.

4.  **Live Plotting (Optional)** * Controlled by the `LIVE_PLOTS` flag (`True`/`False`).
    * If `False`, no windows appear during execution, but all figures are still generated and saved at the end.

5.  **Finalization & Export** Once the optimization is complete, the script:
    * **Best Solution Re-run:** Executes Abaqus one final time with the best parameters found.
    * **File Organization:**
        * Saves exactly **one** numerical "best" `.txt` file in: `./{MODE}_{MATERIAL}_{LAMINATE}/`
        * Saves all plots in: `./{MODE}_{MATERIAL}_{LAMINATE}/Graficas/`
        * Saves the full optimization log in: `./{MODE}_{MATERIAL}_{LAMINATE}/individual_errors.csv`

## 📂 Recommended Repository Structure

Place the following core files in the **repository root**:


├── ga_abaqus_cohesive_fit.py      # Main Genetic Algorithm script
├── Ct_Shell.py                   # Abaqus CAE noGUI script
├── Mechanical_properties.csv      # Material/Case database
└── Experimental/                 # Folder for reference data
    └── Exp_{MODE}_{MATERIAL}_{LAMINATE}.txt

## 📊 Experimental File Format

The experimental data must be provided as a **plain text file** using a semicolon `;` as the separator. It should contain exactly two columns:

x;y
x;y
...
x;y
x;y
...

- x: displacement (e.g., mm)
- y: force (e.g., N)
Empty lines and lines starting with '#' are ignored.

Example:
0.00;0.0
0.10;120.5
0.20;250.1



## Citation

If you use this software in your research or found it interesting for your develpment, please cite the following paper:

**Reference:**
> Cimadevilla-Díez, A., Vaz-Romero, A., Artero-Guerrero, J. A., Pernas-Sánchez, J., Maimí, P., González, E. V., Jacques, V., & De Blanpre, E. (2026) "Genetic algorithm-based optimization for deriving traction–separation laws of CFRPs translaminar fracture at high strain rate". Composite Structures, Vol(No), pages. https://doi.org/10.xxxx/xxxx

**BibTeX:**

@article{lab_code_2024,
  author = {Cimadevilla-Díez, A. and Vaz-Romero, A. and Artero-Guerrero, J. A. and Pernas-Sánchez, J. and Maimí, P. and González, E. V. and Jacques, V. and De Blanpre, E.},
  title = {Genetic algorithm-based optimization for deriving traction–separation laws of CFRPs translaminar fracture at high strain rate},
  journal = {Composite Structures},
  year = {2026},
  doi = {10.xxxx/xxxx},
  url = {[https://doi.org/10.xxxx/xxxx](https://doi.org/10.xxxx/xxxx)}
}


## 🛠️ Requirements

To ensure the calibration script runs correctly, please verify the following prerequisites:

### 🐍 Python Environment
* **Version:** Python 3.9+ (Recommended: **3.10** or **3.11**).
* **Packages:** Install all necessary dependencies using the provided requirements file:
  ```bash
  pip install -r requirements.txt
## Requirements
---------------
- Python 3.9+ (recommended 3.10/3.11)
- Abaqus installed and callable as "abaqus" from terminal or Abaqus Command Prompt
- Python packages:
  pip install -r requirements.txt

IMPORTANT (Windows):
- Most reliable: run from "Abaqus Command Prompt"
  OR set environment variable:
    ABAQUS_BAT = full\path\to\abaqus.bat
- numpy
- matplotlib
- scikit-learn
- deap

## How to run
-------------
1. Verify these exist:
    - Mechanical_properties.csv (repo root)
    - Ct_Shell.py (repo root)
    - ./Experimental/Exp_{MODE}_{MATERIAL}_{LAMINATE}.txt

2. Edit ga_abaqus_cohesive_fit.py:
    - MODE, MATERIAL, LAMINATE
    - POP_SIZE, N_GENERATIONS, bounds, etc.
    - LIVE_PLOTS = False if you want terminal-only execution

3. Run:
    python ga_abaqus_cohesive_fit.py

## Outputs
----------
Inside ./{MODE}_{MATERIAL}_{LAMINATE}/Graficas:
- Final_ForceDisp_vs_Experimental.(png/pdf)
- Final_CohesiveLaw.(png/pdf)
- Evolution_ForceDisp_Live.(png/pdf)
- Evolution_Cohesive_Live.(png/pdf)
- Convergence_*. (png/pdf)
- DIAG_PairPlot_*. (png/pdf)

Inside ./{MODE}_{MATERIAL}_{LAMINATE}/:
- individual_errors.csv
- exactly one NUMERICAL_... file (final best rerun)

## Troubleshooting
------------------
1. "Abaqus launcher not found"
- Run from Abaqus Command Prompt
- Or set ABAQUS_BAT to the full path of abaqus.bat
- Or set abaqus_launcher in the script to an absolute path

2. Abaqus runs but no NUMERICAL_*.txt appears
- Check Ct_Shell.py:
  * must write NUMERICAL_{MODE}_{MATERIAL}_{LAMINATE}*.txt into outdir (repo root)
  * must read the correct CSV row based on READ_FLAG

3. Very large errors / R2 = -inf
- Often caused by insufficient overlap between numerical and experimental displacement ranges.
- Verify disp_min and disp_max in Mechanical_properties.csv.

4. Runtime
- Each evaluation launches Abaqus, so runtime is dominated by simulation cost.
- Increase POP_SIZE / N_GENERATIONS gradually based on your compute budget.

### Good research practice
-------------------------
- Track individual_errors.csv for traceability and reproducibility.
- Document Abaqus version, mesh, boundary conditions, and step definitions as they affect calibration.



