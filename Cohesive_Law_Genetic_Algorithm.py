```python
# =============================================================================
# Genetic Algorithm (DEAP) + Abaqus/CAE curve-fitting for a cohesive law (QS/dynamic)
# =============================================================================
#
# Overview
# --------
# This repository provides an automated inverse-identification workflow to fit a
# cohesive traction–separation law (parameterized by G, Gm, n) to an experimental
# Force–Displacement curve. The optimization is performed with a Genetic Algorithm
# (DEAP), and each candidate parameter set is evaluated by launching Abaqus/CAE
# in no-GUI mode. Abaqus produces a numerical Force–Displacement curve, which is
# compared against the experimental curve through a (densified-only) R² metric.
#
# The objective minimized by the GA is:
#     fitness = 1 - R²
#
# Key features
# ------------
# 1) Fully automated loop:
#    - GA proposes (G, Gm, n)
#    - parameters are written into the material database CSV
#    - Abaqus runs (noGUI) and exports a NUMERICAL_*.txt curve
#    - R² is computed versus the experimental curve
#    - GA evolves the population until the generation limit or early-stop
#
# 2) Physical/consistency constraints:
#    - parameter bounds are enforced
#    - infeasible (G, Gm, n) combinations are penalized (energy consistency, snapback,
#      triangular consistency when n=0, etc.)
#    - a "repair" step clamps and nudges individuals back into feasible regions
#
# 3) Live visualization:
#    - optional real-time plots during evolution (Force–Displacement and cohesive law)
#    - summary plots, convergence plots, and a diagnostic pair-plot are exported at the end
#
# IMPORTANT: R² computation policy (no preprocessing, densification only)
# ----------------------------------------------------------------------
# The experimental curve used for R² is NOT preprocessed:
#   - no sorting
#   - no clipping to a window
#   - no smoothing/filtering
#   - no overlap filtering
#
# The ONLY operation applied is densification of the experimental polyline:
# points are inserted along the straight segments that connect consecutive points
# from the Experimental .txt file.
#
# Then, the numerical curve is interpolated at those densified experimental x-values
# using numpy.interp. Values outside the numerical x-range are clamped to the nearest
# endpoint (left/right).
#
# NOTE: numpy.interp assumes the numerical x-array is monotonic increasing.
# If your Abaqus numerical output is not monotonic in x, the interpolation (and therefore
# the R²) can be invalid. In that case, revise the interpolation strategy or ensure that
# Abaqus exports a monotonic displacement vector.
#
#
# Repository layout (expected)
# ----------------------------
# Place this script in the repository root. The script expects the following files/folders:
#
#   <repo_root>/
#     ├─ Mechanical_properties.csv          # Material DB used by Abaqus script
#     ├─ Ct_Shell.py                        # Abaqus CAE noGUI script (model + export)
#     ├─ Experimental/
#     │    └─ Exp_<MODE>_<MATERIAL>_<LAMINATE>.txt
#     └─ <this_script>.py
#
# During execution, results are written to:
#
#   <repo_root>/<MODE>_<MATERIAL>_<LAMINATE>/
#     ├─ individual_errors.csv              # log of all evaluations
#     ├─ Graficas/                          # PNG/PDF figures (auto-saved)
#     └─ NUMERICAL_...best....txt           # final best curve (only one is kept)
#
#
# Input file formats
# ------------------
# 1) Experimental curve file:
#    - File: Experimental/Exp_<MODE>_<MATERIAL>_<LAMINATE>.txt
#    - Format: 2 columns separated by semicolon ';'
#      Example:
#         0.0000;0.0
#         0.0100;12.3
#         ...
#    - Lines starting with '#' are ignored.
#    - Non-parsable lines are ignored (required to proceed).
#    - IMPORTANT: The point order is preserved (no sorting).
#
# 2) Mechanical_properties.csv:
#    - Must contain a row matching the selected (MODE, MATERIAL, LAMINATE)
#    - The script sets/clears a READ_FLAG column to select the active row for Abaqus.
#    - The script writes G, Gm, n and control flags to the matching row before each Abaqus run.
#
# 3) Abaqus numerical output:
#    - Abaqus is expected to export a file named:
#        NUMERICAL_<MODE>_<MATERIAL>_<LAMINATE>*.txt
#      in the repository root after each run.
#    - The file must be the same ';' two-column format (x;y).
#
#
# Requirements
# ------------
# - Abaqus installed and callable from command line:
#     * Option A (recommended): run this script from an "Abaqus Command Prompt"
#     * Option B: set environment variable ABAQUS_BAT / ABAQUS_CMD / ABAQUS_LAUNCHER
#     * Option C: set `abaqus_launcher` in the config to a full path to abaqus.bat/exe
#
# - Python packages:
#     numpy, matplotlib, scikit-learn, deap
#
#   Example (pip):
#     pip install numpy matplotlib scikit-learn deap
#
# - OS:
#     Windows is explicitly supported for Abaqus launching (cmd.exe wrapping for .bat/.cmd),
#     but the subprocess logic also supports non-Windows launchers.
#
#
# How to run
# ----------
# 1) Configure the case and GA parameters at the top of the script:
#       MODE = "QS"                   # e.g., "QS" or "DYNAMIC" (must match filenames/CSV rows)
#       MATERIAL = "IM7_8552"
#       LAMINATE = "XP"
#       POP_SIZE = 20
#       N_GENERATIONS = 15
#       G_BOUNDS, Gm_BOUNDS, n_BOUNDS
#
# 2) Ensure required files exist:
#       Mechanical_properties.csv
#       Ct_Shell.py
#       Experimental/Exp_QS_IM7_8552_XP.txt   (example)
#
# 3) Launch from a terminal where Abaqus is available:
#       python <this_script>.py
#
# 4) Outputs:
#    - The case folder is created automatically:
#         <MODE>_<MATERIAL>_<LAMINATE>/
#    - Figures are saved into:
#         <case_folder>/Graficas/  as PNG and PDF (dpi=300)
#    - The final best numerical curve is saved in the case folder.
#
#
# Troubleshooting
# ---------------
# - "Abaqus launcher not found (WinError 2)":
#     Run from Abaqus Command Prompt or set ABAQUS_BAT / ABAQUS_CMD / ABAQUS_LAUNCHER,
#     or provide a full path in `abaqus_launcher`.
#
# - No NUMERICAL_*.txt generated:
#     Check Ct_Shell.py paths and export logic. The optimizer will penalize the evaluation
#     when Abaqus fails or when the file is not found.
#
# - R² looks wrong / unstable:
#     Confirm the numerical x-values are monotonic increasing. If not, numpy.interp is not valid.
#     Consider modifying the Abaqus export to output monotonic displacement or implement a different
#     interpolation method for non-monotonic curves.
#
# - CSV READ_FLAG left active:
#     The script clears READ_FLAG markers in a finally-block. If execution is forcibly killed,
#     you may need to manually reset READ_FLAG values to 0.
#
# License / citation
# ------------------
# Add your chosen license and (optionally) a citation request here.
#
# =============================================================================
```


from __future__ import annotations

import ast
import csv
import copy
import os
import random
import shutil
import subprocess
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LogNorm
from sklearn.metrics import r2_score
from deap import base, creator, tools, algorithms

warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
#                           USER CONFIGURATION
# =============================================================================
MODE = "QS"
MATERIAL = "IM7_8552"
LAMINATE = "XP"

POP_SIZE = 20
N_GENERATIONS = 15

G_BOUNDS = (50.0, 250.0)
Gm_BOUNDS = (25.0, 150.0)
n_BOUNDS = (0.0, 1.0)

EARLY_STOP_ERR = 0.01
COHESIVE_XLIM_MAX = 1.5

LEGEND_NCOL = 3
ROUND_G_TO = 0.1
ROUND_N_TO = 0.1

IMMIGRANTS_FRAC_EARLY = 0.20
IMMIGRANTS_FRAC_LATE = 0.05
IMMIGRANTS_SWITCH_GEN_FRAC = 0.40

TRUST_INJECT_EVERY = 5
TRUST_INJECT_START_GEN = 5
TRUST_INJECT_COUNT_FRAC = 0.20
TRUST_SIGMA_FRAC = (0.15, 0.15, 0.08)

DIAG_ERR_CLIP_MAX = 10.0
DIAG_ERR_EPS = 1e-4

SAVE_ALL_NUMERICAL = False  # forced false
LIVE_PLOTS = True

# Plot densification (for aesthetics only)
DENSIFY_EXPERIMENTAL_FOR_PLOTS = True
PLOT_DENSIFY_MAX_POINTS = 1200

# NEW: Densification used for R² (THIS IS THE KEY CHANGE)
# "Solo añadir puntos en la recta que une los puntos ya definidos"
R2_DENSIFY_TARGET_POINTS = 1200  # total points after densification (approx.)
R2_DENSIFY_MIN_POINTS = 50       # minimum densified points

# ------------------------- Units (edit if needed) -----------------------------
UNIT_FORCE = "N"
UNIT_DISP = "mm"
UNIT_SEP = "mm"
UNIT_TRAC = "MPa"
UNIT_G = r"J/m$^2$"
UNIT_DIMLESS = "-"
# =============================================================================


# =============================================================================
#                             PLOTTING STYLE (USER DEFAULTS)
# =============================================================================
def apply_plot_style() -> Tuple[float, float]:
    """
    Consistent plotting style without external LaTeX installation.
    """
    rcParams["mathtext.fontset"] = "stix"
    rcParams["font.family"] = "STIXGeneral"
    rcParams["font.family"] = "serif"

    width = 12
    aurea = 1.618033

    label_size = 24
    rcParams["xtick.labelsize"] = label_size + 5
    rcParams["ytick.labelsize"] = label_size + 5
    rcParams["font.size"] = label_size + 5

    rcParams["figure.facecolor"] = "white"
    rcParams["axes.facecolor"] = "white"
    rcParams["savefig.facecolor"] = "white"

    plt.rcParams.update({
        "text.usetex": False,
        "font.family": "serif",
    })

    return float(width), float(aurea)


def apply_plot_style_compact(font_size: int = 20) -> Tuple[float, float]:
    width, aurea = apply_plot_style()
    rcParams["xtick.labelsize"] = font_size
    rcParams["ytick.labelsize"] = font_size
    rcParams["font.size"] = font_size
    return width, aurea


def setup_axes(ax: plt.Axes) -> None:
    ax.grid(True)
    ax.grid(which="major", linestyle=":", linewidth=0.5)
    ax.minorticks_on()
    ax.grid(which="minor", linestyle=":", linewidth=0.125)


def save_figure(fig: plt.Figure, outbase: Path) -> None:
    outbase = Path(outbase)
    outbase.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(outbase.with_suffix(".png")), dpi=300, bbox_inches="tight")
    fig.savefig(str(outbase.with_suffix(".pdf")), dpi=300, bbox_inches="tight")


def remove_axes_frame(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# =============================================================================
#                               FILE I/O
# =============================================================================
def load_curve_semicolon_2col(path: Path) -> Tuple[np.ndarray, np.ndarray]:
    """
    Minimal parsing of a 2-col ';' file.
    (We must ignore non-parsable lines; otherwise computation is impossible.)
    NO sorting. NO filtering beyond "can be parsed and finite".
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Curve file not found: {path}")

    xs: List[float] = []
    ys: List[float] = []

    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                continue
            parts = [p.strip() for p in s.split(";")]
            if len(parts) < 2:
                continue
            try:
                x = float(parts[0])
                y = float(parts[1])
            except ValueError:
                continue
            if np.isfinite(x) and np.isfinite(y):
                xs.append(x)
                ys.append(y)

    if len(xs) < 3:
        raise ValueError(f"Not enough valid points in: {path}")

    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def densify_polyline_linear(
    x: np.ndarray,
    y: np.ndarray,
    *,
    target_n: int,
    min_n: int = 50,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Densify a polyline by adding points on straight segments between consecutive points.
    NO sorting, NO filtering, preserves the original point order.

    Strategy:
    - Allocate an integer number of points per segment proportional to segment length in x
      (fallback: uniform allocation if total length is ~0).
    - Always includes original vertices.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n_in = int(x.size)

    if n_in < 2:
        return x, y

    target_n = int(max(min_n, target_n))
    # At least original points
    target_n = int(max(target_n, n_in))

    dx = np.diff(x)
    seg_len = np.abs(dx)  # 1D measure; consistent with "add points along x"
    total = float(np.sum(seg_len))

    n_segments = n_in - 1
    # We will add points per segment (including the segment's start point; end points handled by chaining)
    # Total output points is approx target_n.
    n_add_total = target_n - n_in  # number of extra points besides original vertices
    if n_add_total <= 0:
        return x, y

    if total <= 1e-15 or not np.isfinite(total):
        # Uniform allocation if lengths are degenerate
        base = n_add_total // n_segments
        rem = n_add_total % n_segments
        adds = np.array([base + (1 if i < rem else 0) for i in range(n_segments)], dtype=int)
    else:
        w = seg_len / total
        raw = w * float(n_add_total)
        adds = np.floor(raw).astype(int)
        # distribute remainder
        rem = int(n_add_total - int(np.sum(adds)))
        if rem > 0:
            frac = raw - np.floor(raw)
            order = np.argsort(frac)[::-1]
            for k in range(rem):
                adds[order[k % n_segments]] += 1

    x_out: List[float] = [float(x[0])]
    y_out: List[float] = [float(y[0])]

    for i in range(n_segments):
        x0, y0 = float(x[i]), float(y[i])
        x1, y1 = float(x[i + 1]), float(y[i + 1])

        n_add = int(adds[i])
        # Add interior points (excluding x0 and x1). Then append x1,y1 at end of segment.
        if n_add > 0:
            for j in range(1, n_add + 1):
                t = j / float(n_add + 1)
                x_out.append(x0 + t * (x1 - x0))
                y_out.append(y0 + t * (y1 - y0))

        x_out.append(x1)
        y_out.append(y1)

    return np.asarray(x_out, dtype=float), np.asarray(y_out, dtype=float)


def densify_curve_for_plot(
    x: np.ndarray,
    y: np.ndarray,
    *,
    target_n: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Plot-only densification (independent of R²). Uses linear interpolation on a uniform grid.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if x.size < 2:
        return x, y

    n = int(max(50, min(int(target_n), int(PLOT_DENSIFY_MAX_POINTS))))
    # For np.interp, x should be increasing; if not, we fall back to polyline densification for plots
    if np.any(np.diff(x) < 0):
        return densify_polyline_linear(x, y, target_n=n, min_n=50)

    g = np.linspace(float(x.min()), float(x.max()), n)
    yg = np.interp(g, x, y)
    return g, yg


# =============================================================================
#                       CSV HANDLING (READ_FLAG PROTOCOL)
# =============================================================================
READ_FLAG_COL = "READ_FLAG"


@dataclass(frozen=True)
class CaseConfig:
    mode: str
    material: str
    laminate: str


@dataclass(frozen=True)
class MaterialRow:
    E: float
    E_lam: float
    Tt: float
    compl: float
    disp_min: float
    disp_max: float


def _norm_line(line: str) -> str:
    return line.strip().strip('"').strip()


def _is_header_line(line: str) -> bool:
    return _norm_line(line).lower().startswith("ensayo;")


def _read_csv_lines(csv_path: Path) -> List[str]:
    return [_norm_line(l) for l in csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()]


def _write_csv_lines_quoted(csv_path: Path, lines: List[str]) -> None:
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        for line in lines:
            f.write(f"\"{line}\"\n")


def _infer_has_read_flag(data_tokens: List[List[str]], header_tokens: Optional[List[str]]) -> bool:
    if header_tokens and header_tokens[-1].strip().upper() == READ_FLAG_COL:
        return True

    if not data_tokens:
        return False
    lens = [len(t) for t in data_tokens]
    if len(set(lens)) != 1:
        return False
    last_vals = [t[-1].strip() for t in data_tokens if t]
    return bool(last_vals) and all(v in ("0", "1") for v in last_vals)


def ensure_read_flag_column(
    csv_path: Path,
    case: CaseConfig,
    flag_value: int = 1,
    reset_others: bool = True,
) -> None:
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Material DB not found: {csv_path}")

    lines = _read_csv_lines(csv_path)

    header_idx: Optional[int] = None
    header_tokens: Optional[List[str]] = None
    data_rows: List[Tuple[int, List[str]]] = []

    for i, line in enumerate(lines):
        if not line:
            continue
        if _is_header_line(line):
            header_idx = i
            header_tokens = line.split(";")
            continue
        tokens = line.split(";")
        if len(tokens) >= 3:
            data_rows.append((i, tokens))

    data_tokens_only = [t for _, t in data_rows]
    has_flag = _infer_has_read_flag(data_tokens_only, header_tokens)

    if header_tokens is not None and not has_flag:
        header_tokens = header_tokens + [READ_FLAG_COL]
        assert header_idx is not None
        lines[header_idx] = ";".join(header_tokens)

    updated = list(lines)
    found_target = False

    for idx, tokens in data_rows:
        is_target = (tokens[0] == case.mode and tokens[1] == case.material and tokens[2] == case.laminate)
        if is_target:
            found_target = True

        if not has_flag:
            tokens = tokens + ["0"]

        if is_target:
            tokens[-1] = str(int(flag_value))
        elif reset_others:
            tokens[-1] = "0"

        updated[idx] = ";".join(tokens)

    _write_csv_lines_quoted(csv_path, updated)

    if not found_target:
        raise RuntimeError(f"No row found for {case} in {csv_path} while setting READ_FLAG")


def clear_all_read_flags(csv_path: Path) -> None:
    csv_path = Path(csv_path)
    if not csv_path.exists():
        return

    lines = _read_csv_lines(csv_path)
    header_tokens: Optional[List[str]] = None
    data_rows: List[Tuple[int, List[str]]] = []

    for i, line in enumerate(lines):
        if not line:
            continue
        if _is_header_line(line):
            header_tokens = line.split(";")
            continue
        tokens = line.split(";")
        if len(tokens) >= 3:
            data_rows.append((i, tokens))

    data_tokens_only = [t for _, t in data_rows]
    has_flag = _infer_has_read_flag(data_tokens_only, header_tokens)
    if not has_flag:
        return

    updated = list(lines)
    for idx, tokens in data_rows:
        if tokens and tokens[-1].strip() in ("0", "1"):
            tokens[-1] = "0"
            updated[idx] = ";".join(tokens)

    _write_csv_lines_quoted(csv_path, updated)
    print(f"[CSV] Cleared READ_FLAG markers in: {csv_path}")


def read_material_constants(csv_path: Path, case: CaseConfig) -> MaterialRow:
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Material DB not found: {csv_path}")

    ensure_read_flag_column(csv_path, case, flag_value=1, reset_others=True)

    selected_tokens: Optional[List[str]] = None
    with csv_path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            line = _norm_line(raw)
            if not line or _is_header_line(line):
                continue
            tokens = line.split(";")
            if len(tokens) < 19:
                continue
            if tokens[0] == case.mode and tokens[1] == case.material and tokens[2] == case.laminate:
                selected_tokens = tokens
                break

    if selected_tokens is None:
        raise RuntimeError(f"No row found for {case} in {csv_path}")

    return MaterialRow(
        E=float(selected_tokens[6]),
        E_lam=float(selected_tokens[13]),
        Tt=float(selected_tokens[9].split(",")[0]),
        compl=float(selected_tokens[16]),
        disp_min=float(selected_tokens[17]),
        disp_max=float(selected_tokens[18]),
    )


def write_ga_parameters(
    csv_path: Path,
    case: CaseConfig,
    G: float,
    Gm: float,
    n: float,
    eval_index: int,
    end_flag: int,
) -> None:
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Material DB not found: {csv_path}")

    ensure_read_flag_column(csv_path, case, flag_value=1, reset_others=True)
    lines = _read_csv_lines(csv_path)

    updated: List[str] = []
    found = False

    for line in lines:
        if not line or _is_header_line(line):
            updated.append(line)
            continue

        tokens = line.split(";")
        if len(tokens) >= 19 and tokens[0] == case.mode and tokens[1] == case.material and tokens[2] == case.laminate:
            if tokens[-1] not in ("0", "1"):
                tokens.append("1")

            tokens[3] = str(float(G))
            tokens[4] = str(float(Gm))
            tokens[5] = str(float(n))

            tokens[-3] = str(int(eval_index))  # q
            tokens[-2] = str(int(end_flag))    # end
            tokens[-1] = "1"                   # keep marker while running

            updated.append(";".join(tokens))
            found = True
        else:
            updated.append(line)

    if not found:
        raise RuntimeError(f"Could not update GA parameters: row not found for {case}")

    _write_csv_lines_quoted(csv_path, updated)


# =============================================================================
#                     COHESIVE LAW: MODEL + CONSTRAINTS
# =============================================================================
@dataclass(frozen=True)
class Bounds:
    G: Tuple[float, float]
    Gm: Tuple[float, float]
    n: Tuple[float, float]


@dataclass(frozen=True)
class Violations:
    total: float
    v_energy: float
    v_rho0_start: float
    v_snapback: float
    v_triangular_consistency: float


def compute_rho0_rhof_x0(G: float, Gm: float, n: float, row: MaterialRow) -> Tuple[float, float, float]:
    if row.E_lam <= 0.0 or row.Tt <= 0.0:
        raise ValueError("E_lam and Tt must be > 0.")
    Ec = row.E / row.E_lam
    if Ec <= 0.0:
        raise ValueError("Ec must be > 0.")

    T = row.Tt
    Ti = T * n
    x0 = T / Ec

    if Ti <= 1e-15:
        rho0 = x0 + 2.0 * G / T
        rhof = rho0
    else:
        rho0 = x0 + 2.0 * Gm / (T + Ti)
        rhof = rho0 + 2.0 * (G - Gm) / Ti

    return float(rho0), float(rhof), float(x0)


def violations(G: float, Gm: float, n: float, row: MaterialRow) -> Violations:
    eps = 1e-12

    v_energy = 0.0
    if G < -eps or Gm < -eps:
        v_energy += max(0.0, -G) + max(0.0, -Gm) + 1.0
    if Gm > G:
        v_energy += (Gm - G)

    try:
        rho0, rhof, x0 = compute_rho0_rhof_x0(G, Gm, n, row)
    except Exception:
        return Violations(1e6, 1e6, 0.0, 0.0, 0.0)

    v_rho0_start = max(0.0, x0 - rho0)
    v_snapback = max(0.0, rho0 - rhof) if n > 0.0 else 0.0
    v_tri = abs(G - Gm) if n == 0.0 else 0.0

    total = v_energy + v_rho0_start + v_snapback + v_tri
    return Violations(float(total), float(v_energy), float(v_rho0_start), float(v_snapback), float(v_tri))


def sanitize_in_bounds(x: Any, lo: float, hi: float) -> float:
    if isinstance(x, complex):
        x = x.real
    try:
        xf = float(x)
    except Exception:
        xf = float("nan")
    if not np.isfinite(xf):
        xf = random.uniform(lo, hi)
    return float(min(max(xf, lo), hi))


def repair(G: Any, Gm: Any, n: Any, bounds: Bounds, row: MaterialRow) -> Tuple[float, float, float]:
    G = sanitize_in_bounds(G, bounds.G[0], bounds.G[1])
    Gm = sanitize_in_bounds(Gm, bounds.Gm[0], bounds.Gm[1])
    n = sanitize_in_bounds(n, bounds.n[0], bounds.n[1])

    if Gm > G:
        Gm = min(G - 1e-9, max(1e-9, 0.95 * G))

    try:
        rho0, _, x0 = compute_rho0_rhof_x0(G, Gm, n, row)
        if rho0 < x0:
            n = max(bounds.n[0], n - 0.1 * (x0 - rho0) / max(x0, 1e-6))
            n = sanitize_in_bounds(n, bounds.n[0], bounds.n[1])
    except Exception:
        pass

    return float(G), float(Gm), float(n)


def cohesive_curve(G: float, Gm: float, n: float, row: MaterialRow, npts: int = 300) -> Tuple[np.ndarray, np.ndarray]:
    rho0, rhof, x0 = compute_rho0_rhof_x0(G, Gm, n, row)
    T = float(row.Tt)
    Ti = T * float(n)

    x0_ = max(0.0, float(x0))
    rho0_ = max(x0_, float(rho0))
    rhof_ = max(rho0_, float(rhof))

    s_pts = np.array([0.0, x0_, rho0_, rhof_], dtype=float)
    t_pts = np.array([0.0, T,  Ti,   0.0], dtype=float)

    s_grid = np.linspace(0.0, rhof_, int(npts))
    t_grid = np.interp(s_grid, s_pts, t_pts)
    t_grid[t_grid < 0.0] = 0.0
    return s_grid, t_grid


# =============================================================================
#                               SCORING (R²)
# =============================================================================
def score_r2_densified_only(
    exp_x: np.ndarray,
    exp_y: np.ndarray,
    num_x: np.ndarray,
    num_y: np.ndarray,
) -> Tuple[float, float]:
    """
    R² with ONLY densification of experimental polyline.
    - No sorting
    - No window clipping
    - No filtering by overlap
    - Numerical is interpolated at densified experimental x values
      and extrapolated by clamping endpoints (np.interp left/right).
    """
    exp_x = np.asarray(exp_x, dtype=float)
    exp_y = np.asarray(exp_y, dtype=float)
    num_x = np.asarray(num_x, dtype=float)
    num_y = np.asarray(num_y, dtype=float)

    if exp_x.size < 2 or num_x.size < 2:
        return 9.99999999e5, float("-inf")

    # Densify experimental along its consecutive points (polyline).
    xD, yD = densify_polyline_linear(
        exp_x, exp_y,
        target_n=int(R2_DENSIFY_TARGET_POINTS),
        min_n=int(R2_DENSIFY_MIN_POINTS),
    )

    # Interpolate numerical on densified experimental x values.
    # This assumes num_x is increasing. If not, np.interp behavior is not valid.
    # (No sorting here, per user request.)
    yN = np.interp(xD, num_x, num_y, left=float(num_y[0]), right=float(num_y[-1]))

    r2 = float(r2_score(yD, yN))
    one_minus_r2 = float(1.0 - r2)
    return one_minus_r2, r2


# =============================================================================
#                          LABELS / FORMATTING
# =============================================================================
def round_to(x: float, step: float) -> float:
    if step <= 0:
        return float(x)
    return float(step * round(float(x) / step))


def fmt_best_label(G: float, Gm: float, n: float, one_minus_r2: float) -> str:
    Gs = round_to(G, ROUND_G_TO)
    Gms = round_to(Gm, ROUND_G_TO)
    ns = round_to(n, ROUND_N_TO)
    return (
        "Best\n"
        f"G={Gs:.2f}\n"
        f"Gm={Gms:.2f}\n"
        f"n={ns:.2f}\n"
        f"1-R^2={one_minus_r2:.4f}"
    )


# =============================================================================
#                          LIVE PLOTS (EVOLUTION)
# =============================================================================
class LivePlotHistory:
    def __init__(self, width: float, aurea: float, title_prefix: str):
        self.enabled = bool(LIVE_PLOTS)
        if self.enabled:
            plt.ion()

        self.width = width
        self.aurea = aurea
        self.title_prefix = title_prefix
        self.fig: Optional[plt.Figure] = None
        self.ax: Optional[plt.Axes] = None

        self.history: List[Dict[str, Any]] = []
        self.best_idx: Optional[int] = None

    @staticmethod
    def _downsample_xy(x: np.ndarray, y: np.ndarray, max_points: int = 600) -> Tuple[np.ndarray, np.ndarray]:
        if x.size <= max_points:
            return x, y
        idx = np.linspace(0, x.size - 1, max_points).astype(int)
        return x[idx], y[idx]

    def add_case(self, x: np.ndarray, y: np.ndarray, err: float, params: Tuple[float, float, float]) -> None:
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        # For plotting only, keep user request about no preprocessing confined to R².
        # Here we lightly sort to get a clean plot if needed; if you want strictly no sorting anywhere, remove this.
        idx = np.argsort(x)
        x = x[idx]
        y = y[idx]
        x, y = self._downsample_xy(x, y)

        rec = dict(x=x, y=y, err=float(err), params=(float(params[0]), float(params[1]), float(params[2])))
        self.history.append(rec)

        if self.best_idx is None:
            self.best_idx = int(np.argmin([h["err"] for h in self.history]))
        else:
            if rec["err"] < self.history[self.best_idx]["err"]:
                self.best_idx = len(self.history) - 1

    def _ensure_axes(self) -> plt.Axes:
        if self.fig is None or self.ax is None or not plt.fignum_exists(self.fig.number):
            self.fig, self.ax = plt.subplots(figsize=(self.width, self.width / self.aurea))
        return self.ax

    @staticmethod
    def _hist_style(err: float) -> Tuple[float, float]:
        if not np.isfinite(err):
            err = 1.0
        err = float(min(max(err, 0.0), 1.0))
        vis = 1.0 - err
        alpha = 0.06 + 0.10 * vis
        lw = 0.25 + 0.30 * vis
        return alpha, lw


class LivePlotForceDisp(LivePlotHistory):
    def update(
        self,
        x_exp: np.ndarray,
        y_exp: np.ndarray,
        x_cur: Optional[np.ndarray],
        y_cur: Optional[np.ndarray],
        *,
        disp_min: float,
        disp_max: float,
        cur_err: Optional[float],
        highlight_current: bool = True,
    ) -> None:
        if not self.enabled:
            return

        ax = self._ensure_axes()
        ax.clear()
        setup_axes(ax)

        for rec in self.history:
            a, lw = self._hist_style(rec["err"])
            ax.plot(rec["x"], rec["y"], "-", color="0.45", alpha=a, linewidth=lw, label="_nolegend_", zorder=1)

        if self.history and self.best_idx is not None:
            b = self.history[self.best_idx]
            berr = float(b["err"])
            Gb, Gmb, nb = b["params"]
            ax.plot(
                b["x"], b["y"],
                "-", linewidth=1.6, alpha=0.95,
                label=fmt_best_label(Gb, Gmb, nb, berr),
                zorder=5
            )

        if highlight_current and x_cur is not None and y_cur is not None and cur_err is not None:
            idx = np.argsort(x_cur)
            ax.plot(np.asarray(x_cur)[idx], np.asarray(y_cur)[idx], "-", linewidth=1.0, alpha=0.92, label="Current", zorder=6)

        # Plot densification is independent and optional
        if DENSIFY_EXPERIMENTAL_FOR_PLOTS and (x_cur is not None) and (np.asarray(x_cur).size >= 2):
            xE, yE = densify_curve_for_plot(
                x_exp, y_exp,
                target_n=int(min(max(int(np.asarray(x_cur).size), 200), PLOT_DENSIFY_MAX_POINTS)),
            )
            ax.plot(xE, yE, "-", linewidth=1.4, label="Experimental", zorder=7)
        else:
            ax.plot(x_exp, y_exp, "-", linewidth=1.4, label="Experimental", zorder=7)

        ax.axvline(float(disp_min), linestyle="--", linewidth=1.0, label="_nolegend_", zorder=2)
        ax.axvline(float(disp_max), linestyle="--", linewidth=1.0, label="_nolegend_", zorder=2)

        ax.set_title(self.title_prefix)
        ax.set_xlabel(f"Displacement [{UNIT_DISP}]")
        ax.set_ylabel(f"Force [{UNIT_FORCE}]")
        ax.legend(loc="best", ncol=int(LEGEND_NCOL))
        plt.draw()
        plt.pause(0.01)

    def render_summary(self, x_exp: np.ndarray, y_exp: np.ndarray, disp_min: float, disp_max: float) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(self.width, self.width / self.aurea))
        setup_axes(ax)

        for rec in self.history:
            a, lw = self._hist_style(rec["err"])
            ax.plot(rec["x"], rec["y"], "-", color="0.45", alpha=a, linewidth=lw, label="_nolegend_", zorder=1)

        if self.history and self.best_idx is not None:
            b = self.history[self.best_idx]
            Gb, Gmb, nb = b["params"]
            ax.plot(b["x"], b["y"], "-", linewidth=1.8, alpha=0.95, label=fmt_best_label(Gb, Gmb, nb, float(b["err"])), zorder=5)

        if DENSIFY_EXPERIMENTAL_FOR_PLOTS and self.history and self.best_idx is not None:
            xb = np.asarray(self.history[self.best_idx]["x"], dtype=float)
            xE, yE = densify_curve_for_plot(
                x_exp, y_exp,
                target_n=int(min(max(int(xb.size), 200), PLOT_DENSIFY_MAX_POINTS))
            )
            ax.plot(xE, yE, "-", linewidth=1.4, label="Experimental", zorder=6)
        else:
            ax.plot(x_exp, y_exp, "-", linewidth=1.4, label="Experimental", zorder=6)

        ax.axvline(float(disp_min), linestyle="--", linewidth=1.0, label="_nolegend_", zorder=2)
        ax.axvline(float(disp_max), linestyle="--", linewidth=1.0, label="_nolegend_", zorder=2)

        ax.set_xlabel(f"Displacement [{UNIT_DISP}]")
        ax.set_ylabel(f"Force [{UNIT_FORCE}]")
        ax.set_title(self.title_prefix)
        ax.legend(loc="best", ncol=int(LEGEND_NCOL))
        return fig


class LivePlotCohesive(LivePlotHistory):
    def update(
        self,
        s_cur: Optional[np.ndarray],
        t_cur: Optional[np.ndarray],
        *,
        cur_err: Optional[float],
        highlight_current: bool = True,
    ) -> None:
        if not self.enabled:
            return

        ax = self._ensure_axes()
        ax.clear()
        setup_axes(ax)

        for rec in self.history:
            a, lw = self._hist_style(rec["err"])
            ax.plot(rec["x"], rec["y"], "-", color="0.45", alpha=a, linewidth=lw, label="_nolegend_", zorder=1)

        if self.history and self.best_idx is not None:
            b = self.history[self.best_idx]
            berr = float(b["err"])
            Gb, Gmb, nb = b["params"]
            ax.plot(b["x"], b["y"], "-", linewidth=1.6, alpha=0.95, label=fmt_best_label(Gb, Gmb, nb, berr), zorder=5)

        if highlight_current and s_cur is not None and t_cur is not None and cur_err is not None:
            ax.plot(np.asarray(s_cur, dtype=float), np.asarray(t_cur, dtype=float), "-", linewidth=1.0, alpha=0.92, label="Current", zorder=6)

        ax.set_title(self.title_prefix)
        ax.set_xlabel(f"Separation [{UNIT_SEP}]")
        ax.set_ylabel(f"Traction [{UNIT_TRAC}]")
        ax.set_xlim(0.0, float(COHESIVE_XLIM_MAX))
        ax.legend(loc="upper right", ncol=int(LEGEND_NCOL))
        plt.draw()
        plt.pause(0.01)

    def render_summary(self) -> plt.Figure:
        fig, ax = plt.subplots(figsize=(self.width, self.width / self.aurea))
        setup_axes(ax)

        for rec in self.history:
            a, lw = self._hist_style(rec["err"])
            ax.plot(rec["x"], rec["y"], "-", color="0.45", alpha=a, linewidth=lw, label="_nolegend_", zorder=1)

        if self.history and self.best_idx is not None:
            b = self.history[self.best_idx]
            Gb, Gmb, nb = b["params"]
            ax.plot(b["x"], b["y"], "-", linewidth=1.8, alpha=0.95, label=fmt_best_label(Gb, Gmb, nb, float(b["err"])), zorder=5)

        ax.set_xlabel(f"Separation [{UNIT_SEP}]")
        ax.set_ylabel(f"Traction [{UNIT_TRAC}]")
        ax.set_title(self.title_prefix)
        ax.set_xlim(0.0, float(COHESIVE_XLIM_MAX))
        ax.legend(loc="upper right", ncol=int(LEGEND_NCOL))
        return fig


# =============================================================================
#                        ABAQUS LAUNCHING / NUMERICAL I/O
# =============================================================================
@dataclass(frozen=True)
class AbaqusConfig:
    abaqus_launcher: str = "abaqus"
    cae_mode: str = "cae"
    no_gui_flag: str = "noGUI"
    step_name: str = "Step-1"


@dataclass(frozen=True)
class PathsConfig:
    repo_root: Path
    mechanical_csv: Path
    ensayos_dir: Path
    abaqus_script: Path
    case_dir: Path
    figures_dir: Path
    errors_csv: Path


@dataclass(frozen=True)
class GAConfig:
    pop_size: int
    generations: int


@dataclass(frozen=True)
class AppConfig:
    case: CaseConfig
    bounds: Bounds
    ga: GAConfig
    paths: PathsConfig
    abaqus: AbaqusConfig


def _which(cmd: str) -> Optional[str]:
    from shutil import which
    return which(cmd)


def resolve_abaqus_launcher(cfg: AppConfig) -> str:
    cand = cfg.abaqus.abaqus_launcher.strip().strip('"')
    p = Path(cand)
    if p.is_absolute() and p.exists():
        return str(p)

    for envk in ("ABAQUS_BAT", "ABAQUS_CMD", "ABAQUS_LAUNCHER"):
        v = os.environ.get(envk, "").strip().strip('"')
        if v and Path(v).exists():
            return str(Path(v))

    w = _which(cand)
    if w:
        return w

    guesses = [
        r"C:\SIMULIA\Commands\abaqus.bat",
        r"C:\SIMULIA\Commands\abaqus.exe",
        r"C:\Program Files\SIMULIA\Commands\abaqus.bat",
    ]
    for g in guesses:
        if Path(g).exists():
            return g

    raise FileNotFoundError(
        "Abaqus launcher not found (WinError 2).\n"
        "Fix options:\n"
        "1) Run from 'Abaqus Command Prompt', or\n"
        "2) define ABAQUS_BAT with the full path to abaqus.bat, or\n"
        "3) set abaqus_launcher to a full path in the config."
    )


def run_subprocess_windows_aware(args: List[str], cwd: Path) -> subprocess.CompletedProcess:
    launcher = args[0].strip().strip('"')
    ext = Path(launcher).suffix.lower()

    if os.name == "nt" and ext in (".bat", ".cmd"):
        def q(s: str) -> str:
            s = str(s)
            return f"\"{s}\"" if any(c in s for c in (" ", "(", ")", "&", "^", "%", "!")) else s
        cmdline = " ".join(q(a) for a in args)
        proc = subprocess.Popen(["cmd.exe", "/c", cmdline], cwd=str(cwd))
    else:
        proc = subprocess.Popen(args, cwd=str(cwd))

    while True:
        ret = proc.poll()
        if LIVE_PLOTS:
            try:
                plt.pause(0.05)
            except Exception:
                pass
        if ret is not None:
            return subprocess.CompletedProcess(args=args, returncode=int(ret))


def validate_abaqus_available(cfg: AppConfig) -> str:
    launcher = resolve_abaqus_launcher(cfg)
    try:
        _ = run_subprocess_windows_aware([launcher, "information=release"], cfg.paths.repo_root)
    except Exception as e:
        raise RuntimeError(f"Unable to execute Abaqus launcher: {launcher}\nError: {e}")
    return launcher


def call_abaqus(cfg: AppConfig, *, end_flag: int) -> int:
    launcher = resolve_abaqus_launcher(cfg)
    cmd = [
        launcher,
        cfg.abaqus.cae_mode,
        f"{cfg.abaqus.no_gui_flag}={str(cfg.paths.abaqus_script)}",
        "--",
        "--mode", cfg.case.mode,
        "--material", cfg.case.material,
        "--laminate", cfg.case.laminate,
        "--step", cfg.abaqus.step_name,
        "--csv", str(cfg.paths.mechanical_csv),
        "--outdir", str(cfg.paths.repo_root),
        "--end", str(int(end_flag)),
    ]
    proc = run_subprocess_windows_aware(cmd, cfg.paths.repo_root)
    return int(proc.returncode)


def numerical_prefix(cfg: AppConfig) -> str:
    return f"NUMERICAL_{cfg.case.mode}_{cfg.case.material}_{cfg.case.laminate}"


def numerical_name_param(cfg: AppConfig, G: float, Gm: float, n: float) -> str:
    return (
        f"{numerical_prefix(cfg)}"
        f"_Gt_{int(round(G))}_Gm{int(round(Gm))}_n{int(round(n * 100))}.txt"
    )


def safe_copy(src: Path, dst: Path) -> None:
    src = Path(src)
    dst = Path(dst)
    if not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(src), str(dst))


def find_latest_numerical(repo_root: Path, prefix: str) -> Optional[Path]:
    repo_root = Path(repo_root)
    candidates = [p for p in repo_root.glob(f"{prefix}*.txt") if p.is_file()]
    if not candidates:
        return None
    candidates.sort(key=lambda x: x.stat().st_mtime, reverse=True)
    return candidates[0]


def purge_numerical_repo_root(repo_root: Path, prefix: str) -> None:
    repo_root = Path(repo_root)
    for p in repo_root.glob(f"{prefix}*.txt"):
        try:
            if p.is_file():
                p.unlink()
        except Exception:
            pass


def cleanup_repo_root(repo_root: Path) -> None:
    repo_root = Path(repo_root)
    keep_ext = {".py", ".csv", ".txt"}

    deleted = 0
    kept = 0
    for p in repo_root.iterdir():
        try:
            if p.is_dir():
                kept += 1
                continue
            if p.suffix.lower() in keep_ext:
                kept += 1
                continue
            p.unlink()
            deleted += 1
        except Exception as e:
            print(f"[CLEANUP] Could not delete: {p.name} | {e}")

    print(f"[CLEANUP] repo_root={repo_root} | deleted_files={deleted} | kept_entries={kept}")


# =============================================================================
#                      DIAGNOSTICS: CLIPPED LOG COLOR
# =============================================================================
def err_for_diag_color(err: np.ndarray) -> Tuple[np.ndarray, LogNorm]:
    e = np.asarray(err, dtype=float)
    e = np.where(np.isfinite(e) & (e > 0.0), e, np.nan)
    e_col = np.where(np.isfinite(e), np.minimum(e, float(DIAG_ERR_CLIP_MAX)), np.nan)

    finite = e_col[np.isfinite(e_col)]
    if finite.size == 0:
        vmin = float(DIAG_ERR_EPS)
        vmax = float(DIAG_ERR_CLIP_MAX)
    else:
        vmin = float(max(DIAG_ERR_EPS, np.nanpercentile(finite, 5)))
        vmax = float(max(vmin * 10.0, np.nanpercentile(finite, 95)))
        vmax = float(min(vmax, float(DIAG_ERR_CLIP_MAX)))

    norm = LogNorm(vmin=vmin, vmax=vmax)
    return e_col, norm


def diag_pairplot_3vars_compact_3x2(
    X: np.ndarray,
    one_minus_r2: np.ndarray,
    labels: List[str],
    title: str,
    outbase: Path,
    best_idx: int,
) -> None:
    width, aurea = apply_plot_style_compact(font_size=20)

    fig, axes = plt.subplots(3, 2, figsize=(width, width / aurea))
    fig.suptitle(title)

    fig.subplots_adjust(left=0.08, right=0.90, bottom=0.08, top=0.88, wspace=0.28, hspace=0.32)

    err_col, norm = err_for_diag_color(one_minus_r2)

    for i in range(3):
        ax = axes[i, 0]
        setup_axes(ax)
        ax.hist(X[:, i], bins=18)
        ax.axvline(X[best_idx, i], linewidth=2.0, linestyle="--")
        ax.set_xlabel(labels[i])
        ax.set_ylabel("Count [-]")

    scat_axes = [axes[0, 1], axes[1, 1], axes[2, 1]]
    pairs = [(0, 1), (1, 2), (0, 2)]
    sc_last = None

    for ax, (ix, iy) in zip(scat_axes, pairs):
        setup_axes(ax)
        sc_last = ax.scatter(X[:, ix], X[:, iy], c=err_col, s=18, norm=norm)
        ax.scatter(
            X[best_idx, ix], X[best_idx, iy],
            marker="*", s=220, edgecolors="k", linewidths=1.0
        )
        ax.set_xlabel(labels[ix])
        ax.set_ylabel(labels[iy])

    if sc_last is not None:
        cbar = fig.colorbar(sc_last, ax=scat_axes, shrink=0.90, pad=0.03)
        cbar.set_label(f"1-R$^2$ [{UNIT_DIMLESS}] (log scale, clipped)")

    save_figure(fig, outbase)
    plt.close(fig)


# =============================================================================
#                          CONVERGENCE PLOTS
# =============================================================================
def plot_convergence_best_per_gen(best_per_gen: List[float], outbase: Path) -> None:
    width, aurea = apply_plot_style()
    fig, ax = plt.subplots(figsize=(width, width / aurea))
    setup_axes(ax)
    remove_axes_frame(ax)

    x = np.arange(1, len(best_per_gen) + 1, dtype=float)
    y = np.asarray(best_per_gen, dtype=float)
    ax.plot(x, y, "-o")

    y_max = float(np.nanmax(y)) if y.size else 1.0
    if not np.isfinite(y_max) or y_max <= 0.0:
        y_max = 1.0
    ax.set_ylim(0.0, 1.05 * y_max)

    ax.set_xlabel("Generation [-]")
    ax.set_ylabel(f"1-R$^2$ [{UNIT_DIMLESS}]")
    ax.set_title("Convergence per generation")

    save_figure(fig, outbase)
    plt.close(fig)


def plot_convergence_all_from_csv(
    errors_csv: Path,
    outbase: Path,
    *,
    keep_status: Tuple[str, ...] = ("OK", "CACHE"),
    max_reasonable_err: float = 50.0,
) -> None:
    errors_csv = Path(errors_csv)
    if not errors_csv.exists():
        return

    xs: List[int] = []
    ys: List[float] = []

    with errors_csv.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f)
        for rr in r:
            try:
                status = str(rr.get("Status", "")).strip()
                if status not in keep_status:
                    continue

                val = float(rr["one_minus_r2"])
                if (not np.isfinite(val)) or (val < 0.0):
                    continue
                if val > float(max_reasonable_err):
                    continue

                x = int(rr["EvalIndex"])
                xs.append(x)
                ys.append(val)
            except Exception:
                continue

    if len(xs) < 3:
        return

    order = np.argsort(np.asarray(xs, dtype=int))
    xs = list(np.asarray(xs, dtype=int)[order])
    ys = list(np.asarray(ys, dtype=float)[order])

    width, aurea = apply_plot_style()
    fig, ax = plt.subplots(figsize=(width, width / aurea))
    setup_axes(ax)

    ax.plot(xs, ys, "o", markersize=3.0)
    ax.set_xlabel("EvalIndex [-]")
    ax.set_ylabel(f"1-R$^2$ [{UNIT_DIMLESS}]")
    ax.set_title("Convergence of evaluated individuals (penalized excluded)")

    save_figure(fig, outbase)
    plt.close(fig)


# =============================================================================
#                               FINAL PLOTS
# =============================================================================
def plot_final_force_disp(
    exp_x: np.ndarray, exp_y: np.ndarray,
    num_x: np.ndarray, num_y: np.ndarray,
    disp_min: float, disp_max: float,
    G: float, Gm: float, n: float, one_minus_r2: float,
    outbase: Path,
) -> None:
    width, aurea = apply_plot_style()
    fig, ax = plt.subplots(figsize=(width, width / aurea))
    setup_axes(ax)

    idx = np.argsort(num_x)
    num_xs = np.asarray(num_x)[idx]
    num_ys = np.asarray(num_y)[idx]

    if DENSIFY_EXPERIMENTAL_FOR_PLOTS and num_xs.size >= 2:
        xE, yE = densify_curve_for_plot(
            exp_x, exp_y,
            target_n=int(min(max(int(num_xs.size), 200), PLOT_DENSIFY_MAX_POINTS))
        )
        ax.plot(xE, yE, "-", linewidth=1.6, label="Experimental", zorder=3)
    else:
        ax.plot(exp_x, exp_y, "-", linewidth=1.6, label="Experimental", zorder=3)

    ax.plot(num_xs, num_ys, "-", label=fmt_best_label(G, Gm, n, one_minus_r2), zorder=4)
    ax.axvline(float(disp_min), linestyle="--", linewidth=1.0, label="_nolegend_", zorder=2)
    ax.axvline(float(disp_max), linestyle="--", linewidth=1.0, label="_nolegend_", zorder=2)

    ax.set_xlabel(f"Displacement [{UNIT_DISP}]")
    ax.set_ylabel(f"Force [{UNIT_FORCE}]")
    ax.set_title("Final Force-Displacement vs Experimental")
    ax.legend(loc="best", ncol=int(LEGEND_NCOL))
    save_figure(fig, outbase)
    plt.close(fig)


def plot_final_cohesive(row: MaterialRow, G: float, Gm: float, n: float, one_minus_r2: float, outbase: Path) -> None:
    width, aurea = apply_plot_style()
    fig, ax = plt.subplots(figsize=(width, width / aurea))
    setup_axes(ax)

    s, t = cohesive_curve(G, Gm, n, row, npts=500)
    ax.plot(s, t, "-", label=fmt_best_label(G, Gm, n, one_minus_r2), zorder=3)

    ax.set_xlabel(f"Separation [{UNIT_SEP}]")
    ax.set_ylabel(f"Traction [{UNIT_TRAC}]")
    ax.set_title("Final Cohesive Law")
    ax.set_xlim(0.0, float(COHESIVE_XLIM_MAX))
    ax.legend(loc="upper right", ncol=int(LEGEND_NCOL))
    save_figure(fig, outbase)
    plt.close(fig)


# =============================================================================
#                     FINAL: RERUN BEST AND SAVE ONE TXT
# =============================================================================
def rerun_best_and_save_txt(
    cfg: AppConfig,
    row: MaterialRow,
    best_ind: List[float],
    *,
    eval_index_final: int,
) -> Optional[Path]:
    G, Gm, n = map(float, best_ind)
    G, Gm, n = repair(G, Gm, n, cfg.bounds, row)

    write_ga_parameters(cfg.paths.mechanical_csv, cfg.case, G, Gm, n, eval_index_final, end_flag=1)
    ret = call_abaqus(cfg, end_flag=1)

    prefix = numerical_prefix(cfg)
    num_latest = find_latest_numerical(cfg.paths.repo_root, prefix)
    if ret != 0 or num_latest is None or not num_latest.exists():
        print(f"[FINAL BEST RERUN] Abaqus failed or no numerical file. code={ret} | found={num_latest}")
        purge_numerical_repo_root(cfg.paths.repo_root, prefix)
        return None

    best_name = numerical_name_param(cfg, G, Gm, n)
    dst = cfg.paths.case_dir / best_name

    try:
        if dst.exists():
            dst.unlink()
    except Exception:
        pass

    safe_copy(num_latest, dst)
    purge_numerical_repo_root(cfg.paths.repo_root, prefix)

    print(f"[FINAL BEST RERUN] Saved best numerical txt in CASE folder: {dst}")
    return dst


# =============================================================================
#                      GA INITIALIZATION (LHS) + DIVERSITY
# =============================================================================
def init_population_lhs(cfg: AppConfig, row: MaterialRow) -> List[List[float]]:
    N = cfg.ga.pop_size
    (G_lo, G_hi) = cfg.bounds.G
    (Gm_lo, Gm_hi) = cfg.bounds.Gm
    (n_lo, n_hi) = cfg.bounds.n

    bins = np.linspace(0.0, 1.0, N + 1)
    uG = np.array([random.uniform(bins[i], bins[i + 1]) for i in range(N)])
    uGm = np.array([random.uniform(bins[i], bins[i + 1]) for i in range(N)])
    un = np.array([random.uniform(bins[i], bins[i + 1]) for i in range(N)])
    np.random.shuffle(uG)
    np.random.shuffle(uGm)
    np.random.shuffle(un)

    pop: List[List[float]] = []
    for i in range(N):
        G = G_lo + (G_hi - G_lo) * float(uG[i])
        Gm = Gm_lo + (Gm_hi - Gm_lo) * float(uGm[i])
        n = n_lo + (n_hi - n_lo) * float(un[i])

        if Gm > G:
            Gm = 0.90 * G

        G, Gm, n = repair(G, Gm, n, cfg.bounds, row)
        pop.append([G, Gm, n])

    return pop


def diversity(pop: List[List[float]]) -> float:
    if not pop:
        return 0.0
    dim = len(pop[0])
    centroid = [sum(ind[i] for ind in pop) / float(len(pop)) for i in range(dim)]
    d = 0.0
    for ind in pop:
        for i in range(dim):
            d += (ind[i] - centroid[i]) ** 2
    return float(d / len(pop))


# =============================================================================
#                                   MAIN
# =============================================================================
def main() -> None:
    repo_root = Path(__file__).resolve().parent

    cleanup_repo_root(repo_root)

    global SAVE_ALL_NUMERICAL
    if SAVE_ALL_NUMERICAL:
        print("[WARN] SAVE_ALL_NUMERICAL is forced to False (only final best NUMERICAL is kept).")
        SAVE_ALL_NUMERICAL = False

    case = CaseConfig(mode=MODE, material=MATERIAL, laminate=LAMINATE)
    case_name = f"{case.mode}_{case.material}_{case.laminate}"
    case_dir = repo_root / case_name
    figures_dir = case_dir / "Graficas"
    case_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    cfg = AppConfig(
        case=case,
        bounds=Bounds(G=G_BOUNDS, Gm=Gm_BOUNDS, n=n_BOUNDS),
        ga=GAConfig(pop_size=int(POP_SIZE), generations=int(N_GENERATIONS)),
        paths=PathsConfig(
            repo_root=repo_root,
            mechanical_csv=repo_root / "Mechanical_properties.csv",
            ensayos_dir=repo_root / "Experimental",
            abaqus_script=repo_root / "Ct_Shell.py",
            case_dir=case_dir,
            figures_dir=figures_dir,
            errors_csv=case_dir / "individual_errors.csv",
        ),
        abaqus=AbaqusConfig(
            abaqus_launcher="abaqus",
            step_name="Step-1",
        ),
    )

    launcher = validate_abaqus_available(cfg)
    print(f"[INFO] Abaqus launcher resolved as: {launcher}")

    row = read_material_constants(cfg.paths.mechanical_csv, cfg.case)
    exp_file = cfg.paths.ensayos_dir / f"Exp_{cfg.case.mode}_{cfg.case.material}_{cfg.case.laminate}.txt"
    exp_x, exp_y = load_curve_semicolon_2col(exp_file)

    width, aurea = apply_plot_style()
    live_fd = LivePlotForceDisp(width, aurea, title_prefix=f"Evolution Force-Displacement — {case_name}")
    live_coh = LivePlotCohesive(width, aurea, title_prefix=f"Evolution Cohesive Law — {case_name}")

    if not hasattr(creator, "FitnessMin"):
        creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    if not hasattr(creator, "Individual"):
        creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    toolbox.register("clone", copy.deepcopy)

    def mate_safe(ind1, ind2, low, up, eta):
        tools.cxSimulatedBinaryBounded(ind1, ind2, low=low, up=up, eta=eta)
        for i in range(3):
            ind1[i] = sanitize_in_bounds(ind1[i], float(low[i]), float(up[i]))
            ind2[i] = sanitize_in_bounds(ind2[i], float(low[i]), float(up[i]))
        return ind1, ind2

    def mutate_safe(ind, low, up, eta, indpb):
        tools.mutPolynomialBounded(ind, low=low, up=up, eta=eta, indpb=indpb)
        for i in range(3):
            ind[i] = sanitize_in_bounds(ind[i], float(low[i]), float(up[i]))
        return (ind,)

    evaluated_cache: Dict[Tuple[float, float, float], Dict[str, Any]] = {}
    eval_index = 0
    current_generation = 0

    best_global: Optional[creator.Individual] = None
    best_global_one_minus_r2 = float("inf")
    best_global_r2 = float("-inf")

    best_txt_path: Optional[Path] = None
    best_per_gen: List[float] = []

    cfg.paths.errors_csv.parent.mkdir(parents=True, exist_ok=True)
    with cfg.paths.errors_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Individual", "one_minus_r2", "R2", "Fitness", "Generation", "EvalIndex", "Status"])

    def k(ind: List[float]) -> Tuple[float, float, float]:
        return tuple(round(float(x), 12) for x in ind)

    def log_individual(ind: List[float], one_minus_r2: float, r2: float, fitness: float, gen: int, idx: int, status: str) -> None:
        with cfg.paths.errors_csv.open("a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow([list(map(float, ind)), float(one_minus_r2), float(r2), float(fitness), int(gen), int(idx), str(status)])

    def evaluate(individual: creator.Individual) -> Tuple[float]:
        nonlocal eval_index, best_global, best_global_one_minus_r2, best_global_r2

        G, Gm, n = repair(individual[0], individual[1], individual[2], cfg.bounds, row)
        individual[:] = [G, Gm, n]
        key = k(individual)

        if key in evaluated_cache:
            rec = evaluated_cache[key]
            one_minus_r2 = float(rec["one_minus_r2"])
            r2 = float(rec["r2"])
            fitness = float(rec["fitness"])
            log_individual(list(individual), one_minus_r2, r2, fitness, current_generation, eval_index, "CACHE")
            eval_index += 1
            return (fitness,)

        v = violations(G, Gm, n, row)
        if v.total > 1e-9:
            one_minus_r2 = 1.0
            r2 = float("-inf")
            fitness = 1e6 + v.total
            evaluated_cache[key] = dict(one_minus_r2=one_minus_r2, r2=r2, fitness=fitness, status="INFEASIBLE")
            log_individual(list(individual), one_minus_r2, r2, fitness, current_generation, eval_index, "INFEASIBLE")

            live_fd.update(
                exp_x, exp_y, None, None,
                disp_min=row.disp_min, disp_max=row.disp_max,
                cur_err=None,
                highlight_current=False,
            )
            sC, tC = cohesive_curve(G, Gm, n, row, npts=300)
            live_coh.update(sC, tC, cur_err=None, highlight_current=False)

            eval_index += 1
            return (fitness,)

        write_ga_parameters(cfg.paths.mechanical_csv, cfg.case, G, Gm, n, eval_index, end_flag=0)
        ret = call_abaqus(cfg, end_flag=0)

        prefix = numerical_prefix(cfg)
        num_latest = find_latest_numerical(cfg.paths.repo_root, prefix)

        if ret != 0 or num_latest is None or not num_latest.exists():
            one_minus_r2 = 9.99999999e5
            r2 = float("-inf")
            fitness = one_minus_r2
            status = f"ABAQUS_FAIL(code={ret})"
            evaluated_cache[key] = dict(one_minus_r2=one_minus_r2, r2=r2, fitness=fitness, status=status)
            log_individual(list(individual), one_minus_r2, r2, fitness, current_generation, eval_index, status)

            live_fd.update(
                exp_x, exp_y, None, None,
                disp_min=row.disp_min, disp_max=row.disp_max,
                cur_err=None,
                highlight_current=False,
            )
            sC, tC = cohesive_curve(G, Gm, n, row, npts=300)
            live_coh.update(sC, tC, cur_err=None, highlight_current=False)

            purge_numerical_repo_root(cfg.paths.repo_root, prefix)
            eval_index += 1
            return (fitness,)

        num_x, num_y = load_curve_semicolon_2col(num_latest)

        # >>> KEY CHANGE: R² computed using ONLY densification of experimental curve
        one_minus_r2, r2 = score_r2_densified_only(exp_x, exp_y, num_x, num_y)
        fitness = float(one_minus_r2)
        status = "OK"

        evaluated_cache[key] = dict(one_minus_r2=one_minus_r2, r2=r2, fitness=fitness, status=status)
        log_individual(list(individual), one_minus_r2, r2, fitness, current_generation, eval_index, status)

        if np.isfinite(one_minus_r2) and one_minus_r2 < best_global_one_minus_r2:
            best_global_one_minus_r2 = float(one_minus_r2)
            best_global_r2 = float(r2)
            best_global = toolbox.clone(individual)

        live_fd.add_case(num_x, num_y, err=float(one_minus_r2), params=(G, Gm, n))
        live_fd.update(
            exp_x, exp_y, num_x, num_y,
            disp_min=row.disp_min, disp_max=row.disp_max,
            cur_err=one_minus_r2,
            highlight_current=True,
        )

        sC, tC = cohesive_curve(G, Gm, n, row, npts=320)
        live_coh.add_case(sC, tC, err=float(one_minus_r2), params=(G, Gm, n))
        live_coh.update(sC, tC, cur_err=one_minus_r2, highlight_current=True)

        print(
            f"[EVAL {eval_index:05d}] gen={current_generation:02d} | "
            f"G={G:9.4f}  Gm={Gm:9.4f}  n={n:7.4f} | "
            f"1-R^2={one_minus_r2:10.6f}  R2={r2:10.6f} | {status}"
        )

        purge_numerical_repo_root(cfg.paths.repo_root, prefix)

        eval_index += 1
        return (fitness,)

    toolbox.register("evaluate", evaluate)

    pop0 = init_population_lhs(cfg, row)
    population = [creator.Individual(ind) for ind in pop0]

    for ind in population:
        ind.fitness.values = toolbox.evaluate(ind)

    try:
        for gen in range(cfg.ga.generations):
            current_generation = gen

            if best_global is not None and np.isfinite(best_global_one_minus_r2) and best_global_one_minus_r2 <= float(EARLY_STOP_ERR):
                print(f"Early stop: best 1-R^2={best_global_one_minus_r2:.6f} <= {EARLY_STOP_ERR:.6f}")
                break

            factor = gen / float(cfg.ga.generations if cfg.ga.generations > 0 else 1)
            tourn = 2 if factor < 0.40 else (3 if factor < 0.75 else 4)

            cxpb = 0.85 - 0.25 * factor
            mutpb = 0.65 - 0.40 * factor
            indpb = 0.55 - 0.45 * factor
            eta_c = int(8 + (30 - 8) * factor)
            eta_m = int(15 + (60 - 15) * factor)

            op_low = [cfg.bounds.G[0], cfg.bounds.Gm[0], cfg.bounds.n[0]]
            op_up = [cfg.bounds.G[1], cfg.bounds.Gm[1], cfg.bounds.n[1]]

            for op in ("mate", "mutate", "select"):
                if op in toolbox.__dict__:
                    toolbox.unregister(op)

            toolbox.register("mate", mate_safe, low=op_low, up=op_up, eta=eta_c)
            toolbox.register("mutate", mutate_safe, low=op_low, up=op_up, eta=eta_m, indpb=indpb)
            toolbox.register("select", tools.selTournament, tournsize=tourn)

            offspring = algorithms.varAnd(population, toolbox, cxpb, mutpb)
            for ind in offspring:
                ind[:] = list(repair(ind[0], ind[1], ind[2], cfg.bounds, row))

            frac_switch = int(max(1, round(IMMIGRANTS_SWITCH_GEN_FRAC * cfg.ga.generations)))
            imm_frac = float(IMMIGRANTS_FRAC_EARLY) if gen < frac_switch else float(IMMIGRANTS_FRAC_LATE)
            n_imm = max(0, int(round(cfg.ga.pop_size * imm_frac)))
            immigrants: List[creator.Individual] = []
            if n_imm > 0:
                immigrants = [creator.Individual(ind) for ind in init_population_lhs(cfg, row)[:n_imm]]
                for ind in immigrants:
                    ind.fitness.values = toolbox.evaluate(ind)

            injected: List[creator.Individual] = []
            if (
                best_global is not None
                and gen >= int(TRUST_INJECT_START_GEN)
                and TRUST_INJECT_EVERY > 0
                and (gen % int(TRUST_INJECT_EVERY) == 0)
            ):
                n_inj = max(1, int(round(cfg.ga.pop_size * float(TRUST_INJECT_COUNT_FRAC))))
                bg = best_global
                sigG = float(TRUST_SIGMA_FRAC[0]) * (cfg.bounds.G[1] - cfg.bounds.G[0])
                sigGm = float(TRUST_SIGMA_FRAC[1]) * (cfg.bounds.Gm[1] - cfg.bounds.Gm[0])
                sign = float(TRUST_SIGMA_FRAC[2]) * (cfg.bounds.n[1] - cfg.bounds.n[0])

                for _ in range(n_inj):
                    g = np.random.normal(float(bg[0]), sigG)
                    gm = np.random.normal(float(bg[1]), sigGm)
                    nn = np.random.normal(float(bg[2]), sign)
                    g, gm, nn = repair(g, gm, nn, cfg.bounds, row)
                    injected.append(creator.Individual([g, gm, nn]))

                for ind in injected:
                    ind.fitness.values = toolbox.evaluate(ind)

            for ind in offspring:
                if not ind.fitness.valid:
                    ind.fitness.values = toolbox.evaluate(ind)

            if n_imm > 0:
                worst = tools.selWorst(offspring, n_imm)
                for w, im in zip(worst, immigrants):
                    w[:] = im[:]
                    w.fitness.values = im.fitness.values

            if injected:
                worst = tools.selWorst(offspring, len(injected))
                for w, inj in zip(worst, injected):
                    w[:] = inj[:]
                    w.fitness.values = inj.fitness.values

            elit_frac = 0.08 + 0.17 * factor
            elit_n = max(1, int(round(cfg.ga.pop_size * elit_frac)))
            elites = tools.selBest(population, elit_n)

            new_pop = tools.selBest(offspring, cfg.ga.pop_size - elit_n)
            new_pop.extend(elites)
            population[:] = new_pop

            best_per_gen.append(float(best_global_one_minus_r2))

            div = diversity([list(ind) for ind in population])
            print(f"\n--- Generation {gen+1}/{cfg.ga.generations} ---")
            print(f"Diversity: {div:.6f} | Best 1-R^2: {best_global_one_minus_r2:.6f}")

        if best_global is not None:
            best_final_path = rerun_best_and_save_txt(cfg, row, list(best_global), eval_index_final=eval_index)
            if best_final_path is not None:
                best_txt_path = best_final_path

        if best_txt_path is not None and best_txt_path.exists():
            num_x, num_y = load_curve_semicolon_2col(best_txt_path)

            # Final score also uses the same densified-only method (consistent)
            final_one_minus_r2, _final_r2 = score_r2_densified_only(exp_x, exp_y, num_x, num_y)

            G_best, Gm_best, n_best = map(float, best_global) if best_global is not None else (np.nan, np.nan, np.nan)

            plot_final_force_disp(
                exp_x, exp_y, num_x, num_y,
                row.disp_min, row.disp_max,
                G_best, Gm_best, n_best, final_one_minus_r2,
                cfg.paths.figures_dir / "Final_ForceDisp_vs_Experimental"
            )
            plot_final_cohesive(
                row, G_best, Gm_best, n_best, final_one_minus_r2,
                cfg.paths.figures_dir / "Final_CohesiveLaw"
            )

        fig_fd = live_fd.render_summary(exp_x, exp_y, row.disp_min, row.disp_max)
        save_figure(fig_fd, cfg.paths.figures_dir / "Evolution_ForceDisp_Live")
        plt.close(fig_fd)

        fig_coh = live_coh.render_summary()
        save_figure(fig_coh, cfg.paths.figures_dir / "Evolution_Cohesive_Live")
        plt.close(fig_coh)

        plot_convergence_all_from_csv(
            cfg.paths.errors_csv,
            cfg.paths.figures_dir / "Convergence_AllIndividuals",
            keep_status=("OK", "CACHE"),
            max_reasonable_err=50.0,
        )
        plot_convergence_best_per_gen(best_per_gen, cfg.paths.figures_dir / "Convergence_BestPerGeneration")

        X_list: List[List[float]] = []
        E_list: List[float] = []
        with cfg.paths.errors_csv.open("r", encoding="utf-8", newline="") as f:
            r = csv.DictReader(f)
            for rr in r:
                try:
                    ind = ast.literal_eval(rr["Individual"])
                    if len(ind) == 3:
                        X_list.append([float(ind[0]), float(ind[1]), float(ind[2])])
                        E_list.append(float(rr["one_minus_r2"]))
                except Exception:
                    continue

        if len(X_list) >= 5:
            X = np.asarray(X_list, dtype=float)
            E = np.asarray(E_list, dtype=float)
            best_idx = int(np.nanargmin(np.where(np.isfinite(E), E, np.inf)))

            labels = [
                f"G [{UNIT_G}]",
                f"Gm [{UNIT_G}]",
                f"n [{UNIT_DIMLESS}]",
            ]
            diag_pairplot_3vars_compact_3x2(
                X, E, labels,
                title=f"Pair-plot parameters — {case_name}",
                outbase=cfg.paths.figures_dir / "DIAG_PairPlot_G_Gm_n",
                best_idx=best_idx,
            )

    finally:
        clear_all_read_flags(cfg.paths.mechanical_csv)
        cleanup_repo_root(repo_root)

        try:
            for p in cfg.paths.case_dir.iterdir():
                if p.is_dir():
                    continue
                if p.suffix.lower() in (".txt", ".csv"):
                    continue
                try:
                    p.unlink()
                except Exception:
                    pass
        except Exception:
            pass


if __name__ == "__main__":
    main()
