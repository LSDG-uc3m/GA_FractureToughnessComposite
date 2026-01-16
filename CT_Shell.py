# =============================================================================
# CT cohesive Abaqus/CAE script (noGUI) — Model build + run + postprocess export
# =============================================================================
#
# Purpose
# -------
# This script is designed to be executed inside the Abaqus/CAE Python environment
# (typically via `abaqus cae noGUI=<this_file>.py -- <args>`). It automates a full
# simulation pipeline for a CT-like 2D specimen with a cohesive interface:
#
#   1) Reads the *active* case (READ_FLAG=1) from Mechanical_properties.csv
#   2) Builds the Abaqus model (geometry, materials, cohesive law, assembly)
#   3) Runs an Abaqus/Explicit analysis
#   4) Postprocesses the generated ODB
#   5) Exports a numerical curve NUMERICAL_*.txt with (DISP;FORCE)
#
# This script is intended to be called repeatedly by an external optimizer (e.g.
# the DEAP genetic algorithm driver). Each run reads updated cohesive parameters
# (G, Gm, n) from the CSV and produces a corresponding numerical Force–Displacement curve.
#
#
# Expected repository layout
# --------------------------
# Place this file in the repository root. The script expects:
#
#   <repo_root>/
#     ├─ Mechanical_properties.csv                 # Material DB + GA parameters + READ_FLAG
#     ├─ Experimental/
#     │    └─ Exp_<ensayo>_<MATERIAL>_<Laminado>.txt
#     ├─ SHPB_disp/                               # Only required for dynamic tests (HR)
#     │    └─ Disp_<MATERIAL>_<Laminado>.txt       # or supported name variants (see below)
#     └─ Ct_Shell.py (this file)
#
# Output files are written in the working directory (repo root by default):
#   - Abaqus job files: *.inp, *.odb, *.dat, *.sta, etc.
#   - Exported curve:   NUMERICAL_<ensayo>_<MATERIAL>_<Laminado>[_Gt_...].txt
#
#
# How the active case is selected (READ_FLAG protocol)
# ----------------------------------------------------
# Mechanical_properties.csv contains multiple rows (cases). The external driver sets
# one row as active by writing READ_FLAG = 1 (all others 0).
#
# This script:
#   - scans the CSV
#   - finds the first row where the last field is "1"
#   - uses that row to define:
#       ensayo, MATERIAL, Laminado and all model constants
#   - reads the cohesive parameters:
#       Gg (G), m1 (Gm), n1 (n)
#   - reads control flags from the tail of the row:
#       q   = evaluation index (optional, used for logging)
#       end = 0/1 (controls job naming and output naming)
#
# Note: Some repositories store each CSV line as a single quoted field; helper
# functions normalize that into ';'-separated tokens.
#
#
# Supported test modes
# --------------------
# ensayo identifies the test mode:
#
#   - "QS" : quasi-static
#   - "HR" : dynamic high-rate (requested update: "HR" replaces previous "dyn")
#
# QS behavior:
#   - uses a fixed simulation time ti = 0.005 s
#   - uses a SmoothStep amplitude from 0 to 1 over [0, ti]
#   - boundary conditions are symmetric:
#       RP1: u1 = -Ui
#       RP2: u1 = +Ui  (Ut = Ui)
#
# HR behavior:
#   - uses the SHPB displacement-time history to build the loading amplitude
#   - simulation time ti is set to max(time) in the SHPB file
#   - requested constraint: Ut = 0.0 ALWAYS (only Ui is used)
#       RP1: u1 = -Ui (with time-dependent amplitude)
#       RP2: u1 = 0.0 (with the same amplitude, but zero magnitude)
#
#
# Experimental displacement definition (Ui)
# ----------------------------------------
# Ui is derived from the experimental curve file:
#
#   Experimental/Exp_<ensayo>_<MATERIAL>_<Laminado>.txt
#
# The file is expected to contain displacement in the first column (disp;force).
# The script computes:
#   Ui = max(disp) / 2
#
# For HR, Ui is then doubled (Ui = Ui*2.0) to match the requested kinematics:
#   - RP1 moves -Ui while RP2 stays at 0 (Ut=0), effectively applying the full opening.
#
# Non-numeric lines, empty lines, and comment lines (#...) are ignored.
#
#
# Dynamic loading from SHPB (HR)
# ------------------------------
# For HR runs, the script reads the measured displacement-time history and converts it
# into a normalized TabularAmplitude (0 → 1).
#
# SHPB input folder/file:
#   SHPB_disp/
#     - Disp_<MATERIAL>_<Laminado>.txt
#     - Disp<MATERIAL><Laminado>.txt
#     - Disp<MATERIAL><Laminado>.txt (robust concatenation variants)
#
# Expected SHPB file format:
#   - two columns: time, displacement
#   - delimiter can be whitespace/tab OR ';'
#   - optional header lines containing letters are ignored
#
# Processing steps:
#   - reads (t, d), sorts by time, removes repeated times (keeps last value)
#   - shifts time so that t0 = 0
#   - defines ti = max(t)
#   - normalizes displacement:
#       amp(t) = d(t) / max(d)
#     and clamps amp to [0,1], enforcing:
#       amp[0]  = 0
#       amp[-1] = 1
#
# The resulting tuple Velo = ((t0, a0), (t1, a1), ...) is used in:
#   TabularAmplitude(name="Amp-V5", data=Velo)
#
#
# Cohesive law definition
# -----------------------
# The cohesive behavior is defined using:
#   - MaxsDamageInitiation
#   - DamageEvolution (TABULAR, DISPLACEMENT)
#
# The script computes a damage evolution table V3 = [(D, separation), ...] based on:
#   - G  (total energy, from CSV)
#   - Gm (intermediate energy, from CSV)
#   - n  (shape parameter, from CSV)
#   - material constants (E, E_lam, Tt, etc.)
#
# This table is passed to:
#   materials["Material-1"].maxsDamageInitiation.DamageEvolution(..., table=V3)
#
# Notes:
# - The algorithm creates two segments (or one if Ti=0), following the original model logic.
# - Some checks can trigger sys.exit() if invalid geometry/ordering is detected.
#
#
# Model summary (what is built)
# -----------------------------
# - 2D planar parts: CT_0a, CT_0b (two arms) and Cohesivo_1 (cohesive layer)
# - tie constraints connect each arm to one side of the cohesive part
# - an ExplicitDynamicsStep is created with timePeriod = ti
# - history output requests are created on RP1 and RP2 for U1 and RF1
# - boundary conditions are applied at RP1 and RP2 using the configured amplitude
# - meshing is generated (cohesive uses COH2D4 elements, arms use CPS4R/CPS3)
#
#
# Job execution
# -------------
# A job is created and an input file is written. The analysis is launched via:
#   os.system("abq<version> job=<jobname>.inp cpu=12 inter double ask=off")
#
# Ensure that:
# - the correct Abaqus command is available in PATH inside the Abaqus environment
# - abaqus_version is consistent with the installed release
#
#
# Postprocessing and exported numerical curve
# -------------------------------------------
# After the job completes, the script opens the ODB and extracts:
#   - U1 at RP1 and RP2
#   - RF1 at RP1 and RP2
#
# It resolves history region keys robustly by:
#   - reading nodeSets "RP1" and "RP2" from odb.rootAssembly
#   - extracting a node label and instance name
#   - searching the Step-1 historyRegions for matching keys (e.g., "Node <inst>.<label>")
#
# Exported displacement:
#   DISP = -U1(RP1) + U1(RP2)
#
# Exported force:
#   - QS: FORCE = 0.5 * (-RF1(RP1) + RF1(RP2))  # symmetric reaction
#   - HR: FORCE = RF1(RP2)                      # requested behavior in this implementation
#
# The output is saved as a semicolon-separated file:
#   NUMERICAL_<ensayo>_<MATERIAL>_<Laminado>.txt                (end=0)
#   NUMERICAL_<ensayo>_<MATERIAL>_<Laminado>_Gt_<...>_Gm<...>_n<...>.txt  (end=1)
#
# Format:
#   DISP;FORCE
#   with float formatting '%.6f'
#
#
# How to run (typical usage)
# --------------------------
# This script is normally launched by an external driver (GA script) that:
#   - sets READ_FLAG in Mechanical_properties.csv
#   - updates G, Gm, n and control fields in the active row
#   - calls Abaqus:
#       abaqus cae noGUI=Ct_Shell.py -- --mode ... (if args are used)
#
# If running manually, ensure:
#   1) Mechanical_properties.csv has exactly one row with READ_FLAG=1
#   2) Experimental file exists:
#        Experimental/Exp_QS_<MATERIAL>_<Laminado>.txt  or Exp_HR_...
#   3) For HR, SHPB file exists:
#        SHPB_disp/Disp_<MATERIAL>_<Laminado>.txt
#   4) Run from an Abaqus/CAE environment:
#        abaqus cae noGUI=Ct_Shell.py
#
#
# Common failure modes / troubleshooting
# --------------------------------------
# - "No row with READ_FLAG=1":
#     The driver did not set READ_FLAG correctly, or the CSV format is unexpected.
#
# - "Experimental file not found":
#     Verify Experimental/Exp_<ensayo>_<MATERIAL>_<Laminado>.txt exists and naming matches CSV.
#
# - "SHPB displacement file not found" (HR only):
#     Place the correct Disp_<MATERIAL>_<Laminado>.txt in SHPB_disp/ (or match one of the supported names).
#
# - ODB not found:
#     The job may have failed. Inspect *.msg, *.sta, *.dat, or run interactively to diagnose.
#
# - History variables missing (U1/RF1):
#     Ensure HistoryOutputRequests were created and that RP1/RP2 sets exist in the assembly.
#
# =============================================================================

# =========================
# Abaqus imports (CAE macro style)
# =========================
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *

from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup

import section
import regionToolset
import displayGroupMdbToolset as dgm
import displayGroupOdbToolset as dgo
import connectorBehavior
import odbAccess
from odbAccess import *
from abaqus import session

# =========================
# Python / third-party imports
# =========================
import os
import sys
import time
import csv
import inspect
import numpy as np
from math import pi

executeOnCaeStartup()

# =========================
# Paths / configuration
# =========================
name = "CT"
abaqus_version = "2021"

filename = inspect.getframeinfo(inspect.currentframe()).filename
ruta = os.path.dirname(os.path.abspath(filename))
ruta2 = ruta
ruta3 = ruta + "/Experimental"
ruta_excel = ruta + "/Mechanical_properties.csv"


def _to_tokens_from_csv_row(row):
    """
    In some repos, Mechanical_properties.csv lines are stored as a single quoted field;
    csv.reader may return a single column with all ';' inside. This normalizes into tokens.
    """
    if row is None:
        return []

    # Typical case: one "cell" with the entire line
    if len(row) == 1:
        raw = row[0].strip()
        raw = raw.strip('"').strip()
        return [t.strip() for t in raw.split(";")]

    # Alternative case: reader already split by ';'
    return [t.strip().strip('"').strip() for t in row]


def _parse_q_end_from_tokens(tokens):
    """
    With READ_FLAG appended at the end, the tail may be:
        ... ; q ; end ; (other field) ; READ_FLAG
    or:
        ... ; q ; end ; READ_FLAG

    Robust parsing:
    - READ_FLAG: last token if 0/1
    - end: search from the right for the first 0/1 (before the flag)
    - q: token immediately before 'end' (if numeric)
    """
    read_flag = None
    end_val = None
    q_val = None

    if not tokens:
        return read_flag, end_val, q_val

    # READ_FLAG if last token is 0/1
    last = tokens[-1].strip()
    if last in ("0", "1"):
        read_flag = int(last)
        right_limit = len(tokens) - 2
    else:
        right_limit = len(tokens) - 1

    # Find end from the right: first 0/1
    end_idx = None
    for j in range(right_limit, -1, -1):
        tj = tokens[j].strip()
        if tj in ("0", "1"):
            end_val = float(tj)
            end_idx = j
            break

    # q is immediately before end
    if end_idx is not None and end_idx - 1 >= 0:
        try:
            q_val = float(tokens[end_idx - 1].split(",")[0])
        except:
            q_val = 0.0

    return read_flag, end_val, q_val


def load_active_case_from_csv(csv_path):
    """
    Returns (tokens, line_number) for the first row with READ_FLAG == 1.
    """
    f = open(csv_path, "r")
    try:
        reader = csv.reader(f, delimiter=";")
        for i, row in enumerate(reader, start=1):
            tokens = _to_tokens_from_csv_row(row)
            if not tokens:
                continue

            # Skip header
            if tokens[0].lower().startswith("ensayo"):
                continue

            # Active row if last field == 1
            if tokens[-1].strip() == "1":
                return tokens, i
    finally:
        f.close()

    return None, None


def _load_exp_disp_max(exp_path):
    """
    Reads an experimental file 'disp;force' and returns max(disp).
    Skips empty lines, comments (#...), and non-numeric lines.
    """
    if not os.path.exists(exp_path):
        raise IOError("Experimental file not found: %s" % exp_path)

    disp = []
    f = open(exp_path, "r")
    try:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = [p.strip() for p in s.split(";")]
            if len(parts) < 1:
                continue
            try:
                x = float(parts[0].replace(",", "."))
            except:
                continue
            if x == x:  # not NaN
                disp.append(x)
    finally:
        f.close()

    if len(disp) < 3:
        raise ValueError("Not enough valid displacement points in: %s" % exp_path)

    return max(disp)


def _load_shpb_disp_as_tabular_amplitude(base_dir, material, laminate):
    """
    Reads SHPB_disp displacement-time file and returns:
      - Velo: tuple((t, amp_norm), ...) for TabularAmplitude (amp_norm in [0,1])
      - ti  : simulation time = max(t) - t0

    Expected format (robust):
      - 2 columns: time, displacement
      - delimiter can be whitespace/tab or ';'
      - optional header line with text (ignored)

    Files searched (in order):
      1) SHPB_disp/Disp_<material>_<laminate>.txt
      2) SHPB_disp/Disp<material><laminate>.txt
      3) SHPB_disp/Disp + material + laminate + .txt
    """
    folder = os.path.join(base_dir, "SHPB_disp")

    candidates = [
        os.path.join(folder, "Disp_%s_%s.txt" % (material, laminate)),
        os.path.join(folder, "Disp%s%s.txt" % (material, laminate)),
        os.path.join(folder, "Disp" + str(material) + str(laminate) + ".txt"),
    ]

    shpb_path = None
    for p in candidates:
        if os.path.exists(p):
            shpb_path = p
            break

    if shpb_path is None:
        raise IOError(
            "SHPB displacement file not found.\n"
            "Searched:\n  - %s\n  - %s\n  - %s"
            % (candidates[0], candidates[1], candidates[2])
        )

    t_list = []
    d_list = []

    f = open(shpb_path, "r")
    try:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue

            # Skip header lines containing letters
            if any(c.isalpha() for c in s):
                continue

            # Parse delimiter
            if ";" in s:
                parts = [p.strip() for p in s.split(";") if p.strip() != ""]
            else:
                parts = s.split()  # whitespace/tab

            if len(parts) < 2:
                continue

            try:
                t = float(parts[0].replace(",", "."))
                d = float(parts[1].replace(",", "."))
            except:
                continue

            if not (t == t and d == d):
                continue

            t_list.append(t)
            d_list.append(d)
    finally:
        f.close()

    if len(t_list) < 2:
        raise ValueError("Not enough valid points in SHPB file: %s" % shpb_path)

    # Sort by time
    pairs = sorted(zip(t_list, d_list), key=lambda x: x[0])

    # Unique times (keep last displacement for repeated time)
    t_u = []
    d_u = []
    last_t = None
    for t, d in pairs:
        if last_t is None or t != last_t:
            t_u.append(t)
            d_u.append(d)
            last_t = t
        else:
            d_u[-1] = d

    # Shift time to start at 0
    t0 = t_u[0]
    t_u = [t - t0 for t in t_u]
    ti = float(t_u[-1])

    d_max = max(d_u)
    if d_max <= 0.0:
        raise ValueError("SHPB displacement max <= 0 in: %s" % shpb_path)

    # Normalize displacement to amplitude 0..1
    amp = [d / d_max for d in d_u]
    amp = [0.0 if a < 0.0 else (1.0 if a > 1.0 else a) for a in amp]

    # Enforce exact endpoints
    amp[0] = 0.0
    amp[-1] = 1.0

    Velo = tuple((float(t), float(a)) for t, a in zip(t_u, amp))

    print("[SHPB] Loaded amplitude from:", shpb_path)
    print("[SHPB] ti (simulation time) = %.9f s | disp_max = %.6f | npts = %d" % (ti, d_max, len(Velo)))

    return Velo, ti


# --------- Read active case (READ_FLAG=1) ---------
valores_fila, line_no = load_active_case_from_csv(ruta_excel)

if valores_fila is None:
    print("ERROR: No row with READ_FLAG=1 was found in:", ruta_excel)
    sys.exit()

# These variables come from the CSV
ensayo = valores_fila[0]
MATERIAL = valores_fila[1]
Laminado = valores_fila[2]

print("\n[CSV] Active row (READ_FLAG=1) detected:")
print("      File  :", ruta_excel)
print("      Line  :", line_no)
print("      ensayo=", ensayo, "| MATERIAL=", MATERIAL, "| Laminado=", Laminado)
print("      Row   :", valores_fila)
print("")

# --------- Parse remaining constants ---------
Gg = float(valores_fila[3].split(",")[0])
m1 = float(valores_fila[4].split(",")[0])
n1 = float(valores_fila[5].split(",")[0])
E = float(valores_fila[6])
E2 = float(valores_fila[7])
Ec = float(valores_fila[8])
Tt = float(valores_fila[9].split(",")[0])
G12 = float(valores_fila[10])
G23 = float(valores_fila[11])
NU = float(valores_fila[12])
E_lam = float(valores_fila[13])
n_lam = float(valores_fila[14])
n_lam0 = float(valores_fila[15])
compl = float(valores_fila[16])
desp_min = float(valores_fila[17])
desp_max = float(valores_fila[18])

# q and end robust parse with READ_FLAG at the end
_read_flag, end, q = _parse_q_end_from_tokens(valores_fila)

if end is None:
    end = 0.0
if q is None:
    q = 0.0

V1 = 13500.0
mass = 0.0
mesh1 = E_lam * 2.0

mesh2 = 2
eb = 5
D_titanio = 4.51e-09
ec1 = E_lam * n_lam0
ec2 = E_lam * n_lam0
ef = E_lam * n_lam

# === Pre-parameters (must exist): Tt, n1, m1, E, E_lam, Gg ===
Ti = Tt * (n1)
Gm = m1
Ec = E / E_lam
T = Tt
G = Gg

# =========================
# Compute cohesive damage evolution table (as in your script)
# =========================
b = 10000.0
bb = int(b)

x = [0.0] * bb
y = [0.0] * bb
D1 = [0.0] * bb
x2 = [0.0] * bb
y2 = [0.0] * bb
D2 = [0.0] * bb

x3 = [0.0] * bb
x4 = [0.0] * bb

eps = 1e-15

rho_0 = (2.0 * Gm + (T / Ec) * Ti) / T
if Ti == 0.0:
    rho_f = rho_0
else:
    rho_f = 2.0 * (G - Gm) / Ti

den1 = rho_0 - T / Ec
if abs(den1) < eps:
    den1 = eps if den1 >= 0.0 else -eps

if Ti != 0.0:
    den2 = rho_f - rho_0
    if abs(den2) < eps:
        den2 = eps if den2 >= 0.0 else -eps

for a in range(bb):
    # Segment 1
    x[a] = T / Ec + (den1 / b) * (a + 1)
    y[a] = (x[a] - T / Ec) * (Ti - T) / den1 + T
    x3[a] = (den1 / b) * (a + 1)

    denom = Ec * x[a]
    if abs(denom) < eps:
        denom = eps
    D1[a] = 1.0 - y[a] / denom
    if D1[a] < 0.0:
        D1[a] = 0.0
    if D1[a] > 1.0:
        D1[a] = 1.0

    # Segment 2 (only if Ti != 0)
    if Ti != 0.0:
        x2[a] = rho_0 + ((rho_f - rho_0) / b) * (a + 1)
        y2[a] = (x2[a] - rho_0) * (-Ti) / den2 + Ti
        x4[a] = x2[a] - T / Ec

        denom2 = Ec * x2[a]
        if abs(denom2) < eps:
            denom2 = eps
        D2[a] = 1.0 - y2[a] / denom2
        if D2[a] < 0.0:
            D2[a] = 0.0
        if D2[a] > 1.0:
            D2[a] = 1.0
    else:
        x2[a] = rho_0
        y2[a] = 0.0
        x4[a] = 0.0
        D2[a] = 1.0

xx = [0.0]
if Ti == 0.0:
    x5 = xx + x3
    D3 = xx + D1
else:
    x5 = xx + x3 + x4
    D3 = xx + D1 + D2

V3 = list(zip(D3, x5))

# =========================
# Points / areas (same logic)
# =========================
p0a = (0.0, 0.0)
p5a = (T / Ec, T)
p6a = (x3[-1], y[-1])
if Ti == 0.0:
    p7a = (p6a[0], p6a[1])
else:
    p7a = (x4[-1], y2[-1])

area_triangular1 = 0.5 * abs(p5a[0] * p6a[1] - p6a[0] * p5a[1])
area_triangular2 = 0.5 * abs(p7a[0] * p6a[1] - p6a[0] * p7a[1])

p66a = "{0:.3f}".format(p6a[0])
xc = [p0a[0], p5a[0], p6a[0]]
yc = [p0a[1], p5a[1], p6a[1]]

if any(x4_val < x3_val for x3_val in x3 for x4_val in x4):
    sys.exit()

# =========================
# Geometry / general settings
# =========================
w = 51.0
a = 26.0
frames = 1

# =========================
# Load experimental Ui and configure Ut/ti
# =========================
exp_file = os.path.join(ruta3, "Exp_%s_%s_%s.txt" % (ensayo, MATERIAL, Laminado))

# Ui = maximum experimental displacement / 2
Ui = _load_exp_disp_max(exp_file) / 2.0

# QS keeps original fixed time; HR time comes from SHPB displacement file
Velo = None

if ensayo == "QS":
    ti = 0.005
    Ut = Ui

elif ensayo == "HR":
    # Requested behavior: Ut = 0 ALWAYS in HR (only Ui is used)
    Ut = 0.0
    Ui=Ui*2.0
    # Requested behavior: simulation time is max time from SHPB displacement file
    Velo, ti = _load_shpb_disp_as_tabular_amplitude(ruta, MATERIAL, Laminado)

else:
    raise ValueError("Unknown ensayo='%s' (expected 'QS' or 'HR')" % ensayo)

print("[EXP] Using experimental displacement max/2 as Ui: %.6f (file: %s)" % (Ui, exp_file))
print("[CFG] ensayo=%s | Ui=%.6f | Ut=%.6f | ti=%.9f" % (ensayo, Ui, Ut, ti))

# -------------------------------------------------------------------------

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
Mdb()
session.viewports["Viewport: 1"].setValues(displayedObject=None)

mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=200.0)
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(0.0, 0.0), point2=(0.0, 14 + w))
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(0.0, 14 + w), point2=(30.0, 14 + w))
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(30.0, 14 + w), point2=(30.0, 14 + a))
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(30.0, 14 + a), point2=(29.9, 14 + a - 10))
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(29.9, 14 + a - 10), point2=(28.0, 14 + a - 13))
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(28.0, 14 + a - 13), point2=(28.0, 0.0))
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(28.0, 0.0), point2=(0.0, 0.0))

mdb.models["Model-1"].Part(dimensionality=TWO_D_PLANAR, name="CT_0a", type=DEFORMABLE_BODY)
mdb.models["Model-1"].parts["CT_0a"].BaseShell(sketch=mdb.models["Model-1"].sketches["__profile__"])
del mdb.models["Model-1"].sketches["__profile__"]

mdb.models["Model-1"].parts["CT_0a"].Surface(
    name="CT_0a",
    side1Edges=mdb.models["Model-1"].parts["CT_0a"].edges.findAt(((30.0, 57.8125, 0.0),)),
)

mdb.models["Model-1"].Part(name="CT_0b", objectToCopy=mdb.models["Model-1"].parts["CT_0a"])
mdb.models["Model-1"].ConstrainedSketch(
    name="__edit__", objectToCopy=mdb.models["Model-1"].parts["CT_0b"].features["Shell planar-1"].sketch
)
mdb.models["Model-1"].parts["CT_0b"].projectReferencesOntoSketch(
    filter=COPLANAR_EDGES,
    sketch=mdb.models["Model-1"].sketches["__edit__"],
    upToFeature=mdb.models["Model-1"].parts["CT_0b"].features["Shell planar-1"],
)
mdb.models["Model-1"].sketches["__edit__"].mirror(
    mirrorLine=mdb.models["Model-1"].sketches["__edit__"].geometry[4],
    objectList=(
        mdb.models["Model-1"].sketches["__edit__"].geometry[2],
        mdb.models["Model-1"].sketches["__edit__"].geometry[3],
        mdb.models["Model-1"].sketches["__edit__"].geometry[4],
        mdb.models["Model-1"].sketches["__edit__"].geometry[5],
        mdb.models["Model-1"].sketches["__edit__"].geometry[6],
        mdb.models["Model-1"].sketches["__edit__"].geometry[7],
        mdb.models["Model-1"].sketches["__edit__"].geometry[8],
    ),
)
mdb.models["Model-1"].parts["CT_0b"].features["Shell planar-1"].setValues(sketch=mdb.models["Model-1"].sketches["__edit__"])
del mdb.models["Model-1"].sketches["__edit__"]

mdb.models["Model-1"].parts["CT_0b"].regenerate()
mdb.models["Model-1"].parts["CT_0b"].surfaces.changeKey(fromName="CT_0a", toName="CT_0b")

mdb.models["Model-1"].ConstrainedSketch(
    gridSpacing=3.57,
    name="__profile__",
    sheetSize=143.17,
    transform=mdb.models["Model-1"].parts["CT_0a"].MakeSketchTransform(
        sketchPlane=mdb.models["Model-1"].parts["CT_0a"].faces.findAt((18.666667, 9.0, 0.0), (0.0, 0.0, 1.0)),
        sketchPlaneSide=SIDE1,
        sketchOrientation=RIGHT,
        origin=(14.572686, 33.049374, 0.0),
    ),
)
mdb.models["Model-1"].parts["CT_0a"].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models["Model-1"].sketches["__profile__"])
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(3.57, -33.0493739999683 + 14.0), point2=(3.57, -33.0493739999683))
mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((3.57, -28.127187))
mdb.models["Model-1"].sketches["__profile__"].VerticalConstraint(addUndoState=False, entity=mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((3.57, -28.127187)))
mdb.models["Model-1"].sketches["__profile__"].vertices.findAt((3.57, -33.049374))
mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((-0.572686, -33.049374))
mdb.models["Model-1"].sketches["__profile__"].CoincidentConstraint(
    addUndoState=False,
    entity1=mdb.models["Model-1"].sketches["__profile__"].vertices.findAt((3.57, -33.0493739999683)),
    entity2=mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((-0.572686, -33.049374)),
)
mdb.models["Model-1"].parts["CT_0a"].PartitionFaceBySketch(
    faces=mdb.models["Model-1"].parts["CT_0a"].faces.findAt(((18.666667, 9.0, 0.0),)),
    sketch=mdb.models["Model-1"].sketches["__profile__"],
)

mdb.models["Model-1"].ConstrainedSketch(
    gridSpacing=3.57,
    name="__profile__",
    sheetSize=143.17,
    transform=mdb.models["Model-1"].parts["CT_0a"].MakeSketchTransform(
        sketchPlane=mdb.models["Model-1"].parts["CT_0a"].faces.findAt((18.666667, 9.0, 0.0), (0.0, 0.0, 1.0)),
        sketchPlaneSide=SIDE1,
        sketchOrientation=RIGHT,
        origin=(14.572686, 33.049374, 0.0),
    ),
)
mdb.models["Model-1"].parts["CT_0a"].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models["Model-1"].sketches["__profile__"])
mdb.models["Model-1"].sketches["__profile__"].CircleByCenterPerimeter(center=(3.57, -19.0493739999683), point1=(3.57, -19.0493739999683 - 5.0))
mdb.models["Model-1"].parts["CT_0a"].PartitionFaceBySketch(
    faces=mdb.models["Model-1"].parts["CT_0a"].faces.findAt(((18.666667, 9.0, 0.0),)),
    sketch=mdb.models["Model-1"].sketches["__profile__"],
)
del mdb.models["Model-1"].sketches["__profile__"]

mdb.models["Model-1"].ConstrainedSketch(
    gridSpacing=4.42,
    name="__profile__",
    sheetSize=176.91,
    transform=mdb.models["Model-1"].parts["CT_0b"].MakeSketchTransform(
        sketchPlane=mdb.models["Model-1"].parts["CT_0b"].faces.findAt((41.333333, 9.0, 0.0), (0.0, 0.0, 1.0)),
        sketchPlaneSide=SIDE1,
        sketchOrientation=RIGHT,
        origin=(45.427314, 33.049374, 0.0),
    ),
)
mdb.models["Model-1"].parts["CT_0b"].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models["Model-1"].sketches["__profile__"])
mdb.models["Model-1"].sketches["__profile__"].Line(point1=(-2.21, -33.0493739999683 + 14.0), point2=(-2.21, -33.0493739999683))
mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((-2.21, -27.574687))
mdb.models["Model-1"].sketches["__profile__"].VerticalConstraint(addUndoState=False, entity=mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((-2.21, -27.574687)))
mdb.models["Model-1"].sketches["__profile__"].vertices.findAt((-2.21, -33.049374))
mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((0.572686, -33.049374))
mdb.models["Model-1"].sketches["__profile__"].CoincidentConstraint(
    addUndoState=False,
    entity1=mdb.models["Model-1"].sketches["__profile__"].vertices.findAt((-2.21, -33.0493739999683)),
    entity2=mdb.models["Model-1"].sketches["__profile__"].geometry.findAt((0.572686, -33.049374)),
)
mdb.models["Model-1"].parts["CT_0b"].PartitionFaceBySketch(
    faces=mdb.models["Model-1"].parts["CT_0b"].faces.findAt(((41.333333, 9.0, 0.0),)),
    sketch=mdb.models["Model-1"].sketches["__profile__"],
)
del mdb.models["Model-1"].sketches["__profile__"]

mdb.models["Model-1"].ConstrainedSketch(
    gridSpacing=3.57,
    name="__profile__",
    sheetSize=143.17,
    transform=mdb.models["Model-1"].parts["CT_0b"].MakeSketchTransform(
        sketchPlane=mdb.models["Model-1"].parts["CT_0b"].faces.findAt((41.333333, 9.0, 0.0), (0.0, 0.0, 1.0)),
        sketchPlaneSide=SIDE1,
        sketchOrientation=RIGHT,
        origin=(45.427314, 33.049374, 0.0),
    ),
)
mdb.models["Model-1"].parts["CT_0b"].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models["Model-1"].sketches["__profile__"])
mdb.models["Model-1"].sketches["__profile__"].CircleByCenterPerimeter(center=(-2.21, -19.0493739999683), point1=(-2.21, -19.0493739999683 - 5.0))

mdb.models["Model-1"].parts["CT_0b"].PartitionFaceBySketch(
    faces=mdb.models["Model-1"].parts["CT_0b"].faces.findAt(((41.333333, 9.0, 0.0),)),
    sketch=mdb.models["Model-1"].sketches["__profile__"],
)
del mdb.models["Model-1"].sketches["__profile__"]

mdb.models["Model-1"].ConstrainedSketch(name="__profile__", sheetSize=200.0)
mdb.models["Model-1"].sketches["__profile__"].rectangle(point1=(30.0, 14 + a), point2=(30.0 + E_lam, 14 + w))
mdb.models["Model-1"].Part(dimensionality=TWO_D_PLANAR, name="Cohesivo_1", type=DEFORMABLE_BODY)
mdb.models["Model-1"].parts["Cohesivo_1"].BaseShell(sketch=mdb.models["Model-1"].sketches["__profile__"])
del mdb.models["Model-1"].sketches["__profile__"]
mdb.models["Model-1"].parts["Cohesivo_1"].Surface(
    name="C1_a", side1Edges=mdb.models["Model-1"].parts["Cohesivo_1"].edges.findAt(((30.0, 45.9375, 0.0),))
)
mdb.models["Model-1"].parts["Cohesivo_1"].Surface(
    name="C1_b", side1Edges=mdb.models["Model-1"].parts["Cohesivo_1"].edges.findAt(((30.0 + E_lam, 57.8125, 0.0),))
)

# =========================
# Materials
# =========================
mdb.models["Model-1"].Material(name="AS4/8552")
mdb.models["Model-1"].materials["AS4/8552"].Density(table=((1.59e-09,),))
mdb.models["Model-1"].materials["AS4/8552"].Elastic(table=((E, E2, NU, G12, G12, G23),), type=LAMINA)

mdb.models["Model-1"].Material(name="Acero")
mdb.models["Model-1"].materials["Acero"].Density(table=((D_titanio,),))
mdb.models["Model-1"].materials["Acero"].Elastic(table=((E, NU),), type=ISOTROPIC)

mdb.models["Model-1"].Material(name="Material-1")
mdb.models["Model-1"].materials["Material-1"].Density(table=((1.59e-09,),))
mdb.models["Model-1"].materials["Material-1"].Elastic(table=((Ec, G12 / E_lam, G12 / E_lam),), type=TRACTION)
mdb.models["Model-1"].materials["Material-1"].MaxsDamageInitiation(table=((T, T, T),))
mdb.models["Model-1"].materials["Material-1"].maxsDamageInitiation.DamageEvolution(
    softening=TABULAR, table=(V3), type=DISPLACEMENT
)

# =========================
# Sections
# =========================
mdb.models["Model-1"].HomogeneousSolidSection(material="AS4/8552", name="Fibra", thickness=ef)
mdb.models["Model-1"].HomogeneousSolidSection(material="Acero", name="Bulon", thickness=ef)
mdb.models["Model-1"].sections["Fibra"].setValues(material="AS4/8552", thickness=ef)
mdb.models["Model-1"].sections["Bulon"].setValues(material="Acero", thickness=eb)
mdb.models["Model-1"].CohesiveSection(material="Material-1", name="Cohesive-1", outOfPlaneThickness=ec1, response=TRACTION_SEPARATION)

mdb.models["Model-1"].parts["CT_0a"].SectionAssignment(
    offset=0.0,
    offsetField="",
    offsetType=MIDDLE_SURFACE,
    region=Region(faces=mdb.models["Model-1"].parts["CT_0a"].faces.findAt(((1.0, 1.0, 0.0), (0.0, 0.0, 1.0)),)),
    sectionName="Fibra",
    thicknessAssignment=FROM_SECTION,
)
mdb.models["Model-1"].parts["CT_0b"].SectionAssignment(
    offset=0.0,
    offsetField="",
    offsetType=MIDDLE_SURFACE,
    region=Region(faces=mdb.models["Model-1"].parts["CT_0b"].faces.findAt(((39.489755, 6.768161, 0.0), (0.0, 0.0, 1.0)),)),
    sectionName="Fibra",
    thicknessAssignment=FROM_SECTION,
)
mdb.models["Model-1"].parts["Cohesivo_1"].SectionAssignment(
    offset=0.0,
    offsetField="",
    offsetType=MIDDLE_SURFACE,
    region=Region(faces=mdb.models["Model-1"].parts["Cohesivo_1"].faces.findAt(((30.041667, 47.916667, 0.0), (0.0, 0.0, 1.0)),)),
    sectionName="Cohesive-1",
    thicknessAssignment=FROM_SECTION,
)

mdb.models["Model-1"].parts["CT_0a"].MaterialOrientation(
    additionalRotationType=ROTATION_NONE,
    axis=AXIS_3,
    fieldName="",
    localCsys=None,
    orientationType=GLOBAL,
    region=Region(faces=mdb.models["Model-1"].parts["CT_0a"].faces.findAt(((10.156422, 6.768161, 0.0), (0.0, 0.0, 1.0)),)),
    stackDirection=STACK_3,
)
mdb.models["Model-1"].parts["CT_0b"].MaterialOrientation(
    additionalRotationType=ROTATION_NONE,
    axis=AXIS_3,
    fieldName="",
    localCsys=None,
    orientationType=GLOBAL,
    region=Region(faces=mdb.models["Model-1"].parts["CT_0b"].faces.findAt(((39.489755, 6.768161, 0.0), (0.0, 0.0, 1.0)),)),
    stackDirection=STACK_3,
)
mdb.models["Model-1"].parts["CT_0a"].SectionAssignment(
    offset=0.0,
    offsetField="",
    offsetType=MIDDLE_SURFACE,
    region=Region(faces=mdb.models["Model-1"].parts["CT_0a"].faces.findAt(((18.780492, 10.793534, 0.0), (0.0, 0.0, 1.0)),)),
    sectionName="Bulon",
    thicknessAssignment=FROM_SECTION,
)
mdb.models["Model-1"].parts["CT_0b"].SectionAssignment(
    offset=0.0,
    offsetField="",
    offsetType=MIDDLE_SURFACE,
    region=Region(faces=mdb.models["Model-1"].parts["CT_0b"].faces.findAt(((43.85512, 10.793534, 0.0), (0.0, 0.0, 1.0)),)),
    sectionName="Bulon",
    thicknessAssignment=FROM_SECTION,
)

# =========================
# Assembly
# =========================
mdb.models["Model-1"].rootAssembly.DatumCsysByDefault(CARTESIAN)
mdb.models["Model-1"].rootAssembly.Instance(dependent=ON, name="CT_0a-1", part=mdb.models["Model-1"].parts["CT_0a"])
mdb.models["Model-1"].rootAssembly.Instance(dependent=ON, name="CT_0b-1", part=mdb.models["Model-1"].parts["CT_0b"])
mdb.models["Model-1"].rootAssembly.Instance(dependent=ON, name="Cohesivo_1-1", part=mdb.models["Model-1"].parts["Cohesivo_1"])

mdb.models["Model-1"].rootAssembly.translate(instanceList=("CT_0b-1",), vector=(E_lam, 0.0, 0.0))

mdb.models["Model-1"].rootAssembly.Set(
    faces=mdb.models["Model-1"].rootAssembly.instances["CT_0a-1"].faces.findAt(((17.504881, 6.126867, 0.0), (0.0, 0.0, 1.0)),)
    + mdb.models["Model-1"].rootAssembly.instances["CT_0b-1"].faces.findAt(((42.77151, 6.126867, 0.0), (0.0, 0.0, 1.0)),),
    name="Lam",
)

# =========================
# Step
# =========================
mdb.models["Model-1"].rootAssembly.regenerate()

# QS uses fixed ti; HR uses ti from SHPB displacement file
mdb.models["Model-1"].ExplicitDynamicsStep(improvedDtMethod=ON, name="Step-1", previous="Initial", timePeriod=ti)
mdb.models["Model-1"].steps["Step-1"].setValues(improvedDtMethod=ON, scaleFactor=1)

mdb.models["Model-1"].fieldOutputRequests["F-Output-1"].setValues(
    variables=(
        "S",
        "MISES",
        "E",
        "LE",
        "U",
        "UT",
        "UR",
        "V",
        "VT",
        "VR",
        "A",
        "AT",
        "AR",
        "RBANG",
        "RBROT",
        "RF",
        "RT",
        "RM",
        "CF",
        "SF",
        "NFORC",
        "NFORCSO",
        "RBFOR",
        "BF",
        "GRAV",
        "P",
        "HP",
        "IWCONWEP",
        "TRSHR",
        "TRNOR",
        "VP",
        "STAGP",
        "SBF",
        "SDEG",
        "SDV",
        "STATUS",
        "CFAILURE",
        "DMICRT",
        "CSDMG",
        "CSQUADSCRT",
        "STATUS",
    )
)
mdb.models["Model-1"].fieldOutputRequests["F-Output-1"].setValues(numIntervals=frames)

mdb.models["Model-1"].Tie(
    adjust=ON,
    master=mdb.models["Model-1"].rootAssembly.instances["CT_0a-1"].surfaces["CT_0a"],
    name="Constraint-1",
    positionToleranceMethod=COMPUTED,
    slave=mdb.models["Model-1"].rootAssembly.instances["Cohesivo_1-1"].surfaces["C1_a"],
    thickness=ON,
    tieRotations=ON,
)

mdb.models["Model-1"].Tie(
    adjust=ON,
    master=mdb.models["Model-1"].rootAssembly.instances["CT_0b-1"].surfaces["CT_0b"],
    name="Constraint-3",
    positionToleranceMethod=COMPUTED,
    slave=mdb.models["Model-1"].rootAssembly.instances["Cohesivo_1-1"].surfaces["C1_b"],
    thickness=ON,
    tieRotations=ON,
)

mdb.models["Model-1"].rootAssembly.Set(
    name="RP1",
    vertices=mdb.models["Model-1"].rootAssembly.instances["CT_0a-1"].vertices.findAt(((18.142686, 14.0, 0.0),)),
)
mdb.models["Model-1"].rootAssembly.Set(
    name="RP2",
    vertices=mdb.models["Model-1"].rootAssembly.instances["CT_0b-1"].vertices.findAt(((43.217314 + E_lam, 14.0, 0.0),)),
)

# =========================
# Amplitudes (QS smooth step; HR tabular from SHPB)
# =========================
if ensayo == "QS":
    mdb.models["Model-1"].SmoothStepAmplitude(data=((0.0, 0.0), (ti, 1.0)), name="Amp-2", timeSpan=STEP)
    AMPLITUD = "Amp-2"
elif ensayo == "HR":
    mdb.models["Model-1"].TabularAmplitude(data=Velo, name="Amp-V5", smooth=SOLVER_DEFAULT, timeSpan=STEP)
    AMPLITUD = "Amp-V5"
else:
    raise ValueError("Unknown ensayo: %s (expected 'QS' or 'HR')" % ensayo)

# =========================
# Boundary conditions
# =========================
mdb.models["Model-1"].DisplacementBC(
    amplitude=AMPLITUD,
    createStepName="Step-1",
    distributionType=UNIFORM,
    fieldName="",
    fixed=OFF,
    localCsys=None,
    name="Apertura_I",
    region=mdb.models["Model-1"].rootAssembly.sets["RP1"],
    u1=-Ui,
    u2=0,
    u3=0,
    ur1=UNSET,
    ur2=UNSET,
    ur3=UNSET,
)

# QS: Ut = Ui, HR: Ut = 0.0 (RP2 stays at 0)
mdb.models["Model-1"].DisplacementBC(
    amplitude=AMPLITUD,
    createStepName="Step-1",
    distributionType=UNIFORM,
    fieldName="",
    fixed=OFF,
    localCsys=None,
    name="Apertura_i2",
    region=mdb.models["Model-1"].rootAssembly.sets["RP2"],
    u1=Ut,
    u2=0,
    u3=0,
    ur1=UNSET,
    ur2=UNSET,
    ur3=UNSET,
)

# =========================
# Mesh
# =========================
mdb.models["Model-1"].parts["CT_0a"].seedEdgeBySize(
    constraint=FINER,
    deviationFactor=0.1,
    edges=mdb.models["Model-1"].parts["CT_0a"].edges.findAt(((30.0, 57.8125, 0.0),)),
    size=mesh1,
)

mdb.models["Model-1"].parts["CT_0a"].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=mesh2)
mdb.models["Model-1"].parts["CT_0a"].generateMesh()
mdb.models["Model-1"].parts["CT_0b"].seedEdgeBySize(
    constraint=FINER,
    deviationFactor=0.1,
    edges=mdb.models["Model-1"].parts["CT_0b"].edges.findAt(((30.0, 57.8125, 0.0),)),
    size=mesh1,
)
mdb.models["Model-1"].parts["CT_0b"].seedEdgeBySize(
    constraint=FINER,
    deviationFactor=0.1,
    edges=mdb.models["Model-1"].parts["CT_0b"].edges.findAt(((30.0, 57.8125, 0.0),)),
    size=mesh1,
)
mdb.models["Model-1"].parts["CT_0b"].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=mesh2)
mdb.models["Model-1"].parts["CT_0b"].generateMesh()

mdb.models["Model-1"].parts["Cohesivo_1"].seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=mesh1)
mdb.models["Model-1"].parts["Cohesivo_1"].setMeshControls(
    regions=mdb.models["Model-1"].parts["Cohesivo_1"].faces.findAt(((30.041667, 47.916667, 0.0),)), technique=SWEEP
)
mdb.models["Model-1"].parts["Cohesivo_1"].setSweepPath(
    edge=mdb.models["Model-1"].parts["Cohesivo_1"].edges.findAt((30.09375, 40.0, 0.0),),
    region=mdb.models["Model-1"].parts["Cohesivo_1"].faces.findAt((30.041667, 47.916667, 0.0),),
    sense=REVERSE,
)
mdb.models["Model-1"].parts["Cohesivo_1"].generateMesh()

mdb.models["Model-1"].parts["Cohesivo_1"].deleteMesh(regions=mdb.models["Model-1"].parts["Cohesivo_1"].faces.findAt(((30.041667, 47.916667, 0.0),)))
mdb.models["Model-1"].parts["Cohesivo_1"].setMeshControls(elemShape=QUAD, regions=mdb.models["Model-1"].parts["Cohesivo_1"].faces.findAt(((30.041667, 47.916667, 0.0),)))
mdb.models["Model-1"].parts["Cohesivo_1"].generateMesh()
mdb.models["Model-1"].parts["Cohesivo_1"].setElementType(
    elemTypes=(ElemType(elemCode=COH2D4, elemLibrary=EXPLICIT, elemDeletion=ON), ElemType(elemCode=UNKNOWN_TRI, elemLibrary=EXPLICIT)),
    regions=(mdb.models["Model-1"].parts["Cohesivo_1"].faces.findAt(((30.041667, 47.916667, 0.0),)),),
)
mdb.models["Model-1"].parts["CT_0a"].setElementType(
    elemTypes=(
        ElemType(elemCode=CPS4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED, distortionControl=ON, lengthRatio=0.100000001490116),
        ElemType(elemCode=CPS3, elemLibrary=EXPLICIT),
    ),
    regions=(mdb.models["Model-1"].parts["CT_0a"].faces.findAt(((10.156422, 6.768161, 0.0),)),),
)
mdb.models["Model-1"].parts["CT_0b"].setElementType(
    elemTypes=(
        ElemType(elemCode=CPS4R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED, distortionControl=ON, lengthRatio=0.100000001490116),
        ElemType(elemCode=CPS3, elemLibrary=EXPLICIT),
    ),
    regions=(mdb.models["Model-1"].parts["CT_0b"].faces.findAt(((39.489755, 6.768161, 0.0),)),),
)

mdb.models["Model-1"].rootAssembly.regenerate()
mdb.models["Model-1"].HistoryOutputRequest(
    createStepName="Step-1",
    name="H-Output-9",
    numIntervals=1200,
    rebar=EXCLUDE,
    region=mdb.models["Model-1"].rootAssembly.sets["RP1"],
    sectionPoints=DEFAULT,
    variables=("U1", "U2", "U3", "RF1", "RF2", "RF3", "RT"),
)
mdb.models["Model-1"].HistoryOutputRequest(
    createStepName="Step-1",
    name="H-Output-10",
    numIntervals=1200,
    rebar=EXCLUDE,
    region=mdb.models["Model-1"].rootAssembly.sets["RP2"],
    sectionPoints=DEFAULT,
    variables=("U1", "U2", "U3", "RF1", "RF2", "RF3", "RT"),
)

filename = inspect.getframeinfo(inspect.currentframe()).filename
ruta = os.path.dirname(os.path.abspath(filename))
ruta2 = ruta

m11 = int(m1 * 100)
n = int(n1 * 100)
EE1 = int(E)
NU1 = int(100 * NU)
T1 = int(T)
G1 = int(G)
ppi = int(pi * 100)
G12i = int(G12)

if end == 0:
    NOMBRE = "" + name + "_" + ensayo + "_" + MATERIAL + "_" + Laminado
if end == 1:
    NOMBRE = "" + name + "_" + ensayo + "_" + MATERIAL + "_" + Laminado + "Gt_" + str(int(G)) + "_Gm" + str(int(Gm)) + "_n" + str(int(n1 * 100))

os.chdir(r"" + ruta2)
mdb.Job(
    activateLoadBalancing=False,
    atTime=None,
    contactPrint=OFF,
    description="",
    echoPrint=OFF,
    explicitPrecision=DOUBLE,
    historyPrint=OFF,
    memory=90,
    memoryUnits=PERCENTAGE,
    model="Model-1",
    modelPrint=OFF,
    multiprocessingMode=DEFAULT,
    name=NOMBRE,
    nodalOutputPrecision=FULL,
    numCpus=4,
    numDomains=4,
    parallelizationMethodExplicit=DOMAIN,
    queue=None,
    resultsFormat=ODB,
    scratch="",
    type=ANALYSIS,
    userSubroutine="",
    waitHours=0,
    waitMinutes=0,
)
mdb.jobs[NOMBRE].writeInput(consistencyChecking=OFF)

sistem = "abq" + str(abaqus_version) + " job=" + NOMBRE + ".inp cpu=12 inter double ask=off"
os.system(sistem)

# =========================
# POSTPROCESS ODB -> NUMERICAL_*.txt
# Uses RP1/RP2 nodeSets from ODB to build historyRegion keys robustly.
# =========================
odb_path = NOMBRE + ".odb"
if not os.path.exists(odb_path):
    raise IOError("ODB not found: %s" % odb_path)

odb = session.openOdb(name=odb_path)

step_name = "Step-1"
if step_name not in odb.steps.keys():
    odb.close()
    raise RuntimeError("Step not found in ODB: %s" % step_name)


def _first_node_from_nodeset(odb, set_name):
    """
    Returns one OdbMeshNode from odb.rootAssembly.nodeSets[set_name],
    robust across Abaqus versions where nodeSet.nodes may be:
      - a flat list/sequence of OdbMeshNode
      - an OdbMeshNodeArray
      - a sequence of OdbMeshNodeArray, typically per instance
    """
    ra = odb.rootAssembly
    if set_name not in ra.nodeSets.keys():
        raise RuntimeError("NodeSet '%s' not found in ODB. Available: %s" % (set_name, ra.nodeSets.keys()))

    ns = ra.nodeSets[set_name]
    nodes = ns.nodes
    if nodes is None:
        raise RuntimeError("NodeSet '%s' exists but contains no nodes (nodes is None)." % set_name)

    # Case A: flat sequence
    try:
        n0 = nodes[0]
        if hasattr(n0, "label"):
            return n0
    except:
        pass

    # Case B: nested
    try:
        n00 = nodes[0][0]
        if hasattr(n00, "label"):
            return n00
    except:
        pass

    # Case C: iterate
    try:
        for blk in nodes:
            if hasattr(blk, "label"):
                return blk
            try:
                for n in blk:
                    if hasattr(n, "label"):
                        return n
            except:
                pass
    except:
        pass

    raise RuntimeError("Could not extract an OdbMeshNode with attribute 'label' from nodeSet '%s'." % set_name)


def _node_instance_name(node):
    """
    Returns the instance name for an OdbMeshNode across Abaqus versions.
    """
    if hasattr(node, "instanceName"):
        return node.instanceName
    if hasattr(node, "instance") and hasattr(node.instance, "name"):
        return node.instance.name
    return None


def _history_region_key_for_nodeset(odb, set_name, step_name):
    """
    Builds the correct historyRegion key for the node contained in set_name.
    Typical key format:
      "Node CT_0A-1.10"
      "Node CT_0B-1.10"
    """
    node = _first_node_from_nodeset(odb, set_name)
    inst = _node_instance_name(node)
    lab = int(node.label)

    keys = odb.steps[step_name].historyRegions.keys()

    if inst is not None:
        key = "Node %s.%d" % (inst, lab)
        if key in keys:
            return key

        suffix = ".%d" % lab
        for k in keys:
            if k.startswith("Node ") and (inst in k) and k.endswith(suffix):
                return k

    suffix = ".%d" % lab
    cand = [k for k in keys if k.startswith("Node ") and k.endswith(suffix)]
    if len(cand) == 1:
        return cand[0]

    raise RuntimeError(
        "Could not build historyRegion key from nodeSet '%s'.\n"
        "  Node label detected: %d\n"
        "  Instance detected  : %s\n"
        "  Example ODB keys   : %s"
        % (set_name, lab, str(inst), list(keys)[:20])
    )


def _series_from_history(odb, hr_key, var_name, step_name):
    """
    Extracts (values, times) from historyOutputs[var_name] in a given historyRegion.
    """
    step = odb.steps[step_name]
    hr = step.historyRegions[hr_key]

    if var_name not in hr.historyOutputs.keys():
        raise RuntimeError(
            "Variable '%s' not found in historyOutputs for '%s'. Available: %s"
            % (var_name, hr_key, hr.historyOutputs.keys())
        )

    data = hr.historyOutputs[var_name].data  # list of (time, value)
    times = np.array([p[0] for p in data], dtype=float)
    vals = np.array([p[1] for p in data], dtype=float)
    return vals, times


# Resolve historyRegion keys using ODB nodeSets RP1/RP2
hr_key_rp1 = _history_region_key_for_nodeset(odb, "RP1", step_name)
hr_key_rp2 = _history_region_key_for_nodeset(odb, "RP2", step_name)

print("[ODB] RP1 historyRegion key:", hr_key_rp1)
print("[ODB] RP2 historyRegion key:", hr_key_rp2)

# Extract U1 and RF1 from each
U1_rp1, t1 = _series_from_history(odb, hr_key_rp1, "U1", step_name)
U1_rp2, t2 = _series_from_history(odb, hr_key_rp2, "U1", step_name)

RF1_rp1, t3 = _series_from_history(odb, hr_key_rp1, "RF1", step_name)
RF1_rp2, t4 = _series_from_history(odb, hr_key_rp2, "RF1", step_name)

# Truncate to common length
nmin = min(len(U1_rp1), len(U1_rp2), len(RF1_rp1), len(RF1_rp2))
U1_rp1 = U1_rp1[:nmin]
U1_rp2 = U1_rp2[:nmin]
RF1_rp1 = RF1_rp1[:nmin]
RF1_rp2 = RF1_rp2[:nmin]

# DISPLACEMENT = -U1(RP1) + U1(RP2)
DISP = (-U1_rp1 + U1_rp2)

# FORCE
if ensayo == "QS":
    FORCE = 0.5 * (-RF1_rp1 + RF1_rp2)
elif ensayo == "HR":
    FORCE = RF1_rp2
else:
    odb.close()
    raise ValueError("Unknown ensayo: %s (expected 'QS' or 'HR')" % ensayo)

matriz_combinada = np.column_stack((DISP, FORCE))

# Save output
if int(end) == 0:
    outname = "NUMERICAL_%s_%s_%s.txt" % (ensayo, MATERIAL, Laminado)
else:
    outname = "NUMERICAL_%s_%s_%s_Gt_%d_Gm%d_n%d.txt" % (ensayo, MATERIAL, Laminado, int(G), int(Gm), int(n1 * 100))

np.savetxt(outname, matriz_combinada, fmt="%.6f", delimiter=";")
print("[POST] Saved:", outname, "| npts=", len(matriz_combinada))

odb.close()

# Safe cleanup (optional)
for _k in list(session.xyDataObjects.keys()):
    if _k.startswith("Fx_") or _k in ("DISPLACEMENT", "FORCE"):
        try:
            del session.xyDataObjects[_k]
        except:
            pass
 