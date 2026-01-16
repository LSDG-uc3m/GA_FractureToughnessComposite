GA + ABAQUS COHESIVE LAW FITTING
===============================

1) Purpose
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

2) What the script does
-----------------------
For a case defined by (MODE, MATERIAL, LAMINATE), the script:

(1) Loads the experimental curve:
    ./Experimental/Exp_{MODE}_{MATERIAL}_{LAMINATE}.txt

(2) Initializes a population of (G, Gm, n) within user bounds.

(3) For each individual:
    - Writes (G, Gm, n) into the selected row of Mechanical_properties.csv using READ_FLAG protocol
    - Runs Abaqus CAE in noGUI mode using Ct_Shell.py
    - Ct_Shell.py must output a numerical curve file into repo root:
      NUMERICAL_{MODE}_{MATERIAL}_{LAMINATE}*.txt
      (suffixes allowed; Python selects the newest matching file)
    - Computes the fit error:
      raw_err = 1 - R², evaluated over a common displacement window

(4) Live plots (optional)
    - Controlled by LIVE_PLOTS = True/False in the script
    - If False, no live windows appear, but all figures are still saved at the end

(5) End of run:
    - Re-runs Abaqus one final time with the best solution
    - Saves EXACTLY ONE numerical “best” .txt file inside:
      ./{MODE}_{MATERIAL}_{LAMINATE}/
    - Saves figures inside:
      ./{MODE}_{MATERIAL}_{LAMINATE}/Graficas/
    - Saves a CSV log:
      ./{MODE}_{MATERIAL}_{LAMINATE}/individual_errors.csv

3) Recommended repository structure
-----------------------------------
Place these in the repository root:

./ga_abaqus_cohesive_fit.py      (main Python script)
./Ct_Shell.py                   (your Abaqus CAE noGUI script)
./Mechanical_properties.csv      (material/case database)
./Experimental/                 (folder containing experimental curves)
    Exp_{MODE}_{MATERIAL}_{LAMINATE}.txt

When you run the code, a case folder is created:
./{MODE}_{MATERIAL}_{LAMINATE}/
    Graficas/
    individual_errors.csv
    NUMERICAL_{...best...}.txt

4) Experimental file format
---------------------------
The experimental file must be plain text with semicolon ';' separator and 2 columns:

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

5) Requirements
---------------
- Python 3.9+ (recommended 3.10/3.11)
- Abaqus installed and callable as "abaqus" from terminal or Abaqus Command Prompt
- Python packages:
  pip install -r requirements.txt

IMPORTANT (Windows):
- Most reliable: run from "Abaqus Command Prompt"
  OR set environment variable:
    ABAQUS_BAT = full\path\to\abaqus.bat

6) How to run
-------------
(1) Verify these exist:
    - Mechanical_properties.csv (repo root)
    - Ct_Shell.py (repo root)
    - ./Experimental/Exp_{MODE}_{MATERIAL}_{LAMINATE}.txt

(2) Edit ga_abaqus_cohesive_fit.py:
    - MODE, MATERIAL, LAMINATE
    - POP_SIZE, N_GENERATIONS, bounds, etc.
    - LIVE_PLOTS = False if you want terminal-only execution

(3) Run:
    python ga_abaqus_cohesive_fit.py

7) Outputs
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

8) Troubleshooting
------------------
A) "Abaqus launcher not found"
- Run from Abaqus Command Prompt
- Or set ABAQUS_BAT to the full path of abaqus.bat
- Or set abaqus_launcher in the script to an absolute path

B) Abaqus runs but no NUMERICAL_*.txt appears
- Check Ct_Shell.py:
  * must write NUMERICAL_{MODE}_{MATERIAL}_{LAMINATE}*.txt into outdir (repo root)
  * must read the correct CSV row based on READ_FLAG

C) Very large errors / R2 = -inf
- Often caused by insufficient overlap between numerical and experimental displacement ranges.
- Verify disp_min and disp_max in Mechanical_properties.csv.

D) Runtime
- Each evaluation launches Abaqus, so runtime is dominated by simulation cost.
- Increase POP_SIZE / N_GENERATIONS gradually based on your compute budget.

9) Good research practice
-------------------------
- Track individual_errors.csv for traceability and reproducibility.
- Document Abaqus version, mesh, boundary conditions, and step definitions as they affect calibration.

REQUIREMENTS:

numpy
matplotlib
scikit-learn
deap
