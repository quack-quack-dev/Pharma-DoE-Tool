#!/usr/bin/env python3
"""
===============================================================================
Group 5 — Offline Pharmaceutical Powder Characterization Tool (Group_5_tool.py)
===============================================================================
Completely offline laboratory web application (127.0.0.1) for Statistics &
Design of Experiments (HHU Düsseldorf), supporting portable ZIP distribution
with external experiment media loaded dynamically from the relative media/
subdirectory. Implements:
  1. Analytical Sieve Analysis (Q3 cumulative & q3 density distributions,
     linear / logarithmic D10, D50, D90 quantile interpolation, span, relative
     span, and embedded Matplotlib plots).
  2. Bulk & Tapped Density, Carr's Compressibility Index, and Hausner Ratio.
  3. Angle of Repose (fixed-base cone method with optional pedestal subtraction
     and cm/mm SI unit conversion).
  4. Flow Rate Through an Orifice (fixed-mass 100 g, fixed-time 10 s, custom).
  5. Batch Flowability CSV Analyzer (per-run metrics & grouped Mean +/- SD).

Academic & Pharmacopeial Standards:
  - Ph. Eur. 2.9.16: Flowability
  - Ph. Eur. 2.9.34: Bulk Density and Tapped Density of Powders
  - Ph. Eur. 2.9.36: Powder Flow (Tables 2.9.36.-1 & 2.9.36.-2)
  - Ph. Eur. 2.9.38: Particle-Size Distribution Estimation by Analytical Sieving
  - Ph. Eur. 5.3: Statistical Analysis of Results of Biological Assays and Tests
  - USP <429>: Light Diffraction / Particle Size Quantiles & Relative Span
  - ISO 9276-1: Representation of Results of Particle Size Analysis — Part 1
  - ISO 3310-1: Test Sieves — Technical Requirements and Testing
  - ISO 3534-1: Statistics — Vocabulary and Symbols (Mean & Unbiased SD, ddof=1)
  - ISO 80000-1 & IEC 80000-13: Quantities, Units & Information Science Prefixes
  - SI Brochure, 9th ed. (BIPM): Metric Prefixes and Exact Length Conversions
  - IEEE 754-2019: Floating-Point Arithmetic (Finite Real Numbers)
  - IETF RFC 6335 & RFC 4180: Local User Port Range & CSV Parsing
===============================================================================
"""

from __future__ import annotations

# =============================================================================
# DEPENDENCY CHECK & IMPORTS (ZERO-CRASH STARTUP PROTOCOL)
# =============================================================================
# Standard library modules are imported first so check_dependencies() can report
# missing scientific packages cleanly before any third-party import occurs.

import base64
import csv
from http.server import BaseHTTPRequestHandler, HTTPServer
import importlib
import io
import json
import math
from pathlib import Path
import sys
import threading
from urllib.parse import quote, unquote
import webbrowser

REQUIRED_PACKAGES = ("numpy", "pandas", "matplotlib")  # Source: Group 5 Offline Runtime Specification


def check_dependencies() -> None:
    """
    Verifies required third-party scientific packages are installed offline.
    Prints an actionable pip installation command and exits cleanly if missing.
    """
    missing = []
    for pkg in REQUIRED_PACKAGES:
        try:
            importlib.import_module(pkg)
        except ImportError:
            missing.append(pkg)

    if missing:
        print("Error: Missing required Python package(s): " + ", ".join(missing))
        print("Install them once using:")
        print("    pip install " + " ".join(REQUIRED_PACKAGES))
        sys.exit(1)


# Run dependency check BEFORE importing third-party packages (Zero-Crash Protocol)
check_dependencies()

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # Offline non-interactive raster backend for PNG rendering
import matplotlib.pyplot as plt


# =============================================================================
# PHYSICAL, VALIDATION, MEDIA & SERVER CONSTANTS
# =============================================================================

# Source: Python Standard Library (pathlib) — dynamic relative path resolution
# Always resolves media/ relative to this script's runtime location (never hardcoded absolute paths).
SCRIPT_DIR: Path = Path(__file__).resolve().parent
MEDIA_DIR: Path = SCRIPT_DIR / "media"

# Source: IANA Media Types (RFC 6838) — offline MIME mapping for external media/ assets
ALLOWED_MEDIA_TYPES: dict[str, str] = {
    ".svg": "image/svg+xml; charset=utf-8",
    ".png": "image/png",
    ".jpg": "image/jpeg",
    ".jpeg": "image/jpeg",
    ".gif": "image/gif",
    ".webp": "image/webp",
    ".mp4": "video/mp4",
    ".webm": "video/webm",
    ".pdf": "application/pdf",
}

# Source: Group 5 Experiment Module Mapping (Ph. Eur. 2.9.16, 2.9.34, 2.9.36, 2.9.38)
EXPERIMENT_MEDIA_PREFIXES: dict[str, tuple[str, ...]] = {
    "sieve": ("sieve_analysis", "sieve"),
    "carr": ("carr_hausner", "carr", "density", "hausner"),
    "angle": ("angle_of_repose", "angle", "repose"),
    "flow": ("orifice_flow", "flow_rate", "orifice", "flow"),
    "batch": ("batch_flowability", "batch"),
}

# Source: IETF RFC 6335 (local unregistered user port on 127.0.0.1 loopback)
PORT = 8765

# Source: IEC 80000-13 (10 MiB = 10 * 1024^2 bytes payload cap)
MAX_REQUEST_BYTES = 10 * 1024 * 1024

# Source: Ph. Eur. 2.9.16 (Flowability): fixed test mass = 100 g
ORIFICE_FIXED_MASS_G = 100.0

# Source: Ph. Eur. 2.9.36 (Flow through an orifice): fixed-time method = 10 s
ORIFICE_FIXED_TIME_S = 10.0

# Source: SI Brochure, 9th ed. (BIPM), Table 7 (centi = 10^-2, milli = 10^-3): 1 cm = 10 mm
MM_PER_CM = 10.0

# Source: Ph. Eur. 2.9.38 / ISO 9276-1: Default assumed upper particle-size boundary
# above the coarsest sieve (1600 µm) in the HHU Düsseldorf PD01 sieve stack
DEFAULT_X_MAX_UM = 2000.0

# Source: IEEE 754-2019 / ISO 9276-1 — finite upper particle-size bound (1e12 µm = 10^6 m)
# well above physical granular materials and safely inside Matplotlib LogLocator range
MAX_SIEVE_SIZE_UM = 1e12

# Source: Ph. Eur. 2.9.38 / ISO 3310-1: A valid analytical sieve stack requires at
# least one wire-mesh sieve plus the collecting pan (minimum 2 stages)
MIN_SIEVE_STAGES = 2

# Source: Ph. Eur. 2.9.38 / ISO 9276-1: Collecting pan represents aperture x = 0.0 µm
PAN_MESH_UM = 0.0

# Source: HHU Düsseldorf Lab Protocol (PD01) / Analytical balance readability precision
# (rounding to 6 decimal places eliminates IEEE 754 binary subtraction artifacts)
RETAINED_MASS_DECIMALS = 6

# Source: SI Brochure, 9th ed. (BIPM): Supported length units for angle of repose measurements
SUPPORTED_LENGTH_UNITS = ("cm", "mm")

# Source: Ph. Eur. 2.9.38 / HHU Düsseldorf PD01 CSV schema for sieve gravimetric data
REQUIRED_SIEVE_COLUMNS = ("mesh_um", "empty_g", "full_g")

# Source: Ph. Eur. 2.9.34 & 2.9.36 / HHU Düsseldorf PD01 CSV schema for flowability batch data
REQUIRED_FLOWABILITY_COLUMNS = (
    "sample_code",
    "group_run",
    "mass_g",
    "v0_ml",
    "vt_ml",
    "total_height_mm",
    "h_base_mm",
    "diameter_mm",
)


# =============================================================================
# MODULE A: DATA VALIDATION, UNIT HELPERS, EMBEDDED SAMPLE DATA & CSV PARSING
# =============================================================================

# =============================================================================
# EMBEDDED OFFLINE SAMPLE DATASETS (RULE #1: ZERO EXTERNAL FILE DEPENDENCIES)
# =============================================================================

# Source: HHU Düsseldorf Pharmaceutical Technology Lab Dataset — sieve_PD01.csv
SAMPLE_SIEVE_PD01_CSV = (
    "mesh_um,empty_g,full_g\n"
    "1600,390.0,391.8\n"
    "1400,365.3,368.0\n"
    "1000,357.6,375.4\n"
    "800,323.9,329.3\n"
    "630,320.4,324.5\n"
    "400,336.6,351.3\n"
    "315,285.5,286.8\n"
    "125,388.7,393.2\n"
    "0,721.1,726.9\n"
)

# Source: HHU Düsseldorf Pharmaceutical Technology Lab Dataset — sieve_reference.csv
SAMPLE_SIEVE_REFERENCE_CSV = (
    "mesh_um,empty_g,full_g\n"
    "1600,389.9,391.7\n"
    "1400,365.4,370.8\n"
    "1000,357.5,373.0\n"
    "800,323.9,330.5\n"
    "630,320.6,325.3\n"
    "400,336.3,343.4\n"
    "315,285.5,287.2\n"
    "125,388.7,392.6\n"
    "0,721.1,724.1\n"
)

# Source: HHU Düsseldorf Pharmaceutical Technology Lab Dataset — flowability_complete_dataset.csv
SAMPLE_FLOWABILITY_CSV = (
    "sample_code,group_run,mass_g,v0_ml,vt_ml,total_height_mm,h_base_mm,diameter_mm\n"
    "PD01-EFF (Effervescent),G3 (n=1),100.4,143.0,121.0,54.0,28.0,100.0\n"
    "PD01-EFF (Effervescent),G3 (n=2),102.5,142.0,116.0,55.0,28.0,100.0\n"
    "PD01-EFF (Effervescent),G4 (n=1),100.0,138.0,126.0,54.0,28.0,100.0\n"
    "PD01-EFF (Effervescent),G4 (n=2),100.0,132.0,124.0,55.0,28.0,100.0\n"
    "PD01-CON-S1 (L/S=0.1),G3 (n=1),75.6,240.0,149.0,90.0,28.0,100.0\n"
    "PD01-CON-S1 (L/S=0.1),G4 (n=1),78.0,240.0,157.0,90.0,28.0,100.0\n"
    "PD01-CON-S1 (L/S=0.1),G4 (n=2),78.0,240.0,154.0,87.0,28.0,100.0\n"
    "PD01-CON-S2 (L/S=0.3),G3 (n=1),72.9,200.0,160.0,85.0,28.0,100.0\n"
    "PD01-CON-S2 (L/S=0.3),G4 (n=1),60.0,150.0,140.0,86.0,28.0,100.0\n"
    "Product 3: High-Shear F1 (PVP K30),G3 (n=1),100.5,188.0,152.0,54.0,28.0,100.0\n"
    "Product 3: High-Shear F1 (PVP K30),G4 (n=1),100.0,176.0,148.0,54.0,28.0,100.0\n"
    "Product 3: High-Shear F1 (PVP K30),G4 (n=2),100.0,182.0,155.0,54.0,28.0,100.0\n"
    "Product 3: High-Shear F2 (PVP K90),G4 (n=1),100.0,190.0,153.0,54.0,28.0,100.0\n"
)


# =============================================================================
# INTERNAL HELPERS
# =============================================================================

def _coerce_numeric_sequence(seq, field_name: str) -> list:
    """
    Validates that an input is an ordered iterable sequence (not a scalar, string,
    bytes, dict, or unordered set) and converts it to a Python list.
    """
    # Source: Defensive Input Validation (Rule #2) — reject strings, bytes, mappings, and unordered sets
    if seq is None or isinstance(seq, (str, bytes, bytearray, dict, set, frozenset)):
        raise ValueError(
            f"'{field_name}' must be an ordered sequence of numbers (e.g., list or tuple), "
            f"received {type(seq).__name__}."
        )
    try:
        return list(seq)
    except TypeError as exc:
        raise ValueError(
            f"'{field_name}' must be an ordered sequence of numbers, "
            f"received non-iterable {type(seq).__name__}."
        ) from exc


def _read_csv_table(csv_text: str, required_columns: tuple[str, ...]) -> pd.DataFrame:
    """
    Validates raw CSV text, strips UTF-8 BOM, normalizes line endings, parses all
    cells as strings without silent NA coercion, and verifies required columns.
    """
    if not isinstance(csv_text, str):
        raise ValueError(
            f"CSV input must be a text string, received {type(csv_text).__name__}."
        )

    # Source: RFC 4180 / Unicode Standard — strip optional UTF-8 BOM and normalize CRLF/CR to LF
    cleaned = csv_text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")
    if not cleaned.strip():
        raise ValueError("CSV input is empty or contains only whitespace.")

    try:
        # Source: RFC 4180 / Zero-Assumption Protocol (Rule #2) — validate exact column counts
        # and duplicate required headers via stdlib csv.reader before pandas can mangle duplicates
        parsed_rows = [
            row
            for row in csv.reader(io.StringIO(cleaned))
            if row and not (len(row) == 1 and not row[0].strip())
        ]
    except csv.Error as exc:
        raise ValueError(f"Failed to parse CSV text: {exc}") from exc

    if not parsed_rows:
        raise ValueError("CSV input is empty or contains only whitespace.")

    # Strip leading/trailing whitespace and any stray BOM on header names
    raw_columns = [str(col).lstrip("\ufeff").strip() for col in parsed_rows[0]]
    expected_cols = len(raw_columns)

    # Detect duplicate required headers BEFORE pd.read_csv so pandas cannot rename them to '.1'
    for col in required_columns:
        if raw_columns.count(col) > 1:
            raise ValueError(f"CSV contains duplicate required column '{col}'.")

    missing = [col for col in required_columns if col not in set(raw_columns)]
    if missing:
        raise ValueError(
            f"CSV is missing required column(s): {', '.join(missing)}. "
            f"Found columns: {raw_columns}."
        )

    if len(parsed_rows) < 2:
        raise ValueError("CSV contains headers but no data rows.")

    for row_num, row_cells in enumerate(parsed_rows[1:], start=2):
        if len(row_cells) != expected_cols:
            raise ValueError(
                f"CSV line {row_num} has {len(row_cells)} column(s), expected {expected_cols}."
            )

    try:
        # Source: Zero-Assumption Protocol (Rule #2) — dtype=str, keep_default_na=False & index_col=False prevent silent coercion
        df = pd.read_csv(io.StringIO(cleaned), dtype=str, keep_default_na=False, index_col=False)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeError, ValueError) as exc:
        raise ValueError(f"Failed to parse CSV text: {exc}") from exc

    df.columns = raw_columns

    if df.empty:
        raise ValueError("CSV contains headers but no data rows.")

    return df


# =============================================================================
# PUBLIC VALIDATION & PARSING FUNCTIONS
# =============================================================================

def parse_finite_float(value, field_name: str) -> float:
    """
    Strictly parses a scalar value into a finite Python float.

    Rejects None, booleans, empty strings, ambiguous mixed decimal/thousand
    separators, multiple decimal separators, non-numeric values, NaN, and +/-Inf.
    Safely accepts a single European decimal comma (e.g., '12,5' -> 12.5).
    """
    # Source: Zero-Assumption Protocol (Rule #2) — reject missing values, booleans, and raw byte buffers
    # (in Python bool is a subclass of int, so float(True) would silently yield 1.0)
    if value is None or isinstance(value, (bool, bytes, bytearray)) or type(value).__name__ == "bool_":
        raise ValueError(
            f"Invalid value for '{field_name}': expected a finite number, received {value!r}."
        )

    if isinstance(value, str):
        text = value.strip()
        if text == "":
            raise ValueError(f"Invalid value for '{field_name}': input cannot be empty.")

        # Source: ISO 80000-1, Section 7.3.2 — decimal sign is either a point or a comma,
        # never both simultaneously and never repeated within a single number
        has_comma = "," in text
        has_dot = "." in text
        if has_comma and has_dot:
            raise ValueError(
                f"Invalid value for '{field_name}' ({value!r}): ambiguous number format containing "
                "both ',' and '.'. Remove thousand separators and use a single decimal separator."
            )
        if text.count(",") > 1 or text.count(".") > 1:
            raise ValueError(
                f"Invalid value for '{field_name}' ({value!r}): multiple decimal separators are not allowed."
            )
        if "_" in text:
            raise ValueError(
                f"Invalid value for '{field_name}' ({value!r}): underscores are not valid in numeric inputs."
            )

        # Source: ISO 80000-1, Section 7.3.2 — convert unambiguous European decimal comma to decimal point
        if has_comma:
            text = text.replace(",", ".")

        try:
            num = float(text)
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError(
                f"Invalid value for '{field_name}': could not parse {value!r} as a number."
            ) from exc
    else:
        try:
            num = float(value)
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError(
                f"Invalid value for '{field_name}': could not parse {value!r} as a number."
            ) from exc

    # Source: IEEE 754-2019 — physical pharmaceutical measurements must be finite real numbers (not NaN/Inf)
    if not math.isfinite(num):
        raise ValueError(
            f"Invalid value for '{field_name}': value must be finite (received {value!r})."
        )

    # Normalize negative zero (-0.0) to canonical 0.0
    return 0.0 if num == 0.0 else num


def validate_sieve_data(
    mesh_sizes,
    empty_weights,
    full_weights,
    x_max_um=DEFAULT_X_MAX_UM,
) -> dict:
    """
    Validates gravimetric sieve analysis inputs and computes retained masses.

    Returns a dictionary containing cleaned lists of mesh sizes [µm], empty sieve
    weights [g], full sieve weights [g], retained masses [g], and upper class
    boundary x_max_um [µm].
    """
    mesh_raw = _coerce_numeric_sequence(mesh_sizes, "mesh_sizes")
    empty_raw = _coerce_numeric_sequence(empty_weights, "empty_weights")
    full_raw = _coerce_numeric_sequence(full_weights, "full_weights")

    if not (len(mesh_raw) == len(empty_raw) == len(full_raw)):
        raise ValueError(
            "Input sequences 'mesh_sizes', 'empty_weights', and 'full_weights' must have equal lengths "
            f"(received {len(mesh_raw)}, {len(empty_raw)}, {len(full_raw)})."
        )

    # Source: Ph. Eur. 2.9.38 / ISO 3310-1 — minimum 2 stages (at least 1 test sieve + collecting pan)
    if len(mesh_raw) < MIN_SIEVE_STAGES:
        raise ValueError(
            f"At least {MIN_SIEVE_STAGES} sieve stages (at least one mesh sieve and the 0 µm pan) "
            f"are required (received {len(mesh_raw)})."
        )

    x_max_clean = parse_finite_float(x_max_um, "x_max_um")
    # Source: IEEE 754-2019 / ISO 9276-1 — enforce safe finite upper aperture bound
    if x_max_clean > MAX_SIEVE_SIZE_UM:
        raise ValueError(
            f"Upper size bound x_max_um ({x_max_clean:g} µm) exceeds maximum supported "
            f"particle size limit ({MAX_SIEVE_SIZE_UM:g} µm)."
        )

    mesh_clean: list[float] = []
    empty_clean: list[float] = []
    full_clean: list[float] = []

    for idx, (m_raw, e_raw, f_raw) in enumerate(zip(mesh_raw, empty_raw, full_raw)):
        stage = idx + 1
        m_val = parse_finite_float(m_raw, f"Stage {stage} mesh_um")
        e_val = parse_finite_float(e_raw, f"Stage {stage} empty_g")
        f_val = parse_finite_float(f_raw, f"Stage {stage} full_g")

        # Source: Ph. Eur. 2.9.38 / SI physical non-negativity — apertures and masses cannot be negative
        if m_val < 0.0:
            raise ValueError(f"Stage {stage}: mesh size cannot be negative ({m_val} µm).")
        if e_val < 0.0:
            raise ValueError(f"Stage {stage}: empty sieve weight cannot be negative ({e_val} g).")
        if f_val < 0.0:
            raise ValueError(f"Stage {stage}: full sieve weight cannot be negative ({f_val} g).")

        # Source: Ph. Eur. 2.9.38 — retained mass m_i = m_full,i - m_empty,i must satisfy m_i >= 0
        if f_val < e_val:
            raise ValueError(
                f"Stage {stage} ({m_val} µm): full weight ({f_val} g) is less than empty weight "
                f"({e_val} g). Retained mass cannot be negative."
            )

        mesh_clean.append(m_val)
        empty_clean.append(e_val)
        full_clean.append(f_val)

    # Source: Ph. Eur. 2.9.38 / ISO 9276-1 — sieves must be ordered strictly decreasing from top sieve to pan
    for idx in range(len(mesh_clean) - 1):
        if mesh_clean[idx] <= mesh_clean[idx + 1]:
            raise ValueError(
                f"Mesh sizes must be strictly decreasing from coarsest top sieve to pan. "
                f"Stage {idx + 1} ({mesh_clean[idx]} µm) is not greater than "
                f"Stage {idx + 2} ({mesh_clean[idx + 1]} µm)."
            )

    # Source: Ph. Eur. 2.9.38 / ISO 9276-1 — final row must be the collecting pan (aperture = 0 µm)
    if mesh_clean[-1] != PAN_MESH_UM:
        raise ValueError(
            f"Final row must be the collecting pan with mesh size {PAN_MESH_UM} µm "
            f"(received {mesh_clean[-1]} µm)."
        )

    # Source: ISO 9276-1 — upper boundary x_max must strictly exceed top sieve aperture so dx_0 = x_max - x_0 > 0
    if x_max_clean <= mesh_clean[0]:
        raise ValueError(
            f"Upper size bound x_max_um ({x_max_clean} µm) must be strictly greater than "
            f"the largest sieve mesh size ({mesh_clean[0]} µm)."
        )

    # Source: Ph. Eur. 2.9.38 — gravimetric retained mass m_i = m_full,i - m_empty,i
    raw_retained = [f - e for f, e in zip(full_clean, empty_clean)]
    retained_masses = [round(m, RETAINED_MASS_DECIMALS) for m in raw_retained]

    # Source: Ph. Eur. 2.9.38 / ISO 9276-1 — total retained mass m_total = sum(m_i) must be > 0
    if sum(raw_retained) <= 0.0 or sum(retained_masses) <= 0.0:
        raise ValueError(
            "Total retained powder mass across all sieves is 0 g. "
            "At least one sieve or pan must retain powder (> 0 g)."
        )

    return {
        "mesh_sizes": mesh_clean,
        "empty_weights": empty_clean,
        "full_weights": full_clean,
        "retained_masses": retained_masses,
        "x_max_um": x_max_clean,
    }


def validate_carr_hausner_data(mass_g, v0_ml, vf_ml) -> tuple[float, float, float]:
    """
    Validates inputs for Bulk Density, Tapped Density, Carr's Compressibility Index,
    and Hausner Ratio.

    Returns (mass_g, v0_ml, vf_ml) as finite floats.
    """
    mass_clean = parse_finite_float(mass_g, "mass_g")
    v0_clean = parse_finite_float(v0_ml, "v0_ml")
    vf_clean = parse_finite_float(vf_ml, "vf_ml")

    # Source: Ph. Eur. 2.9.34 (Bulk Density and Tapped Density of Powders) — m > 0, V_0 > 0, V_f > 0
    if mass_clean <= 0.0:
        raise ValueError(f"Powder mass 'mass_g' must be > 0 g (received {mass_clean} g).")
    if v0_clean <= 0.0:
        raise ValueError(f"Unsettled bulk volume 'v0_ml' must be > 0 mL (received {v0_clean} mL).")
    if vf_clean <= 0.0:
        raise ValueError(f"Tapped volume 'vf_ml' must be > 0 mL (received {vf_clean} mL).")

    # Source: Ph. Eur. 2.9.34 & 2.9.36 — tapping consolidates the powder bed, requiring V_0 >= V_f
    if v0_clean < vf_clean:
        raise ValueError(
            f"Unsettled volume v0_ml ({v0_clean} mL) cannot be smaller than tapped volume "
            f"vf_ml ({vf_clean} mL). Ensure V_0 >= V_f."
        )

    return (mass_clean, v0_clean, vf_clean)


def validate_angle_of_repose_data(
    total_height,
    diameter,
    base_height=0.0,
    unit="cm",
) -> dict:
    """
    Validates inputs for the Angle of Repose (fixed-base cone method) with optional
    pedestal base-height subtraction and SI unit conversion between cm and mm.

    Returns a dictionary containing validated dimensions in both the input unit
    and normalized cm/mm values.
    """
    if not isinstance(unit, str):
        raise ValueError(
            f"Length 'unit' must be a string ('cm' or 'mm'), received {type(unit).__name__}."
        )
    unit_clean = unit.strip().lower()
    if unit_clean not in SUPPORTED_LENGTH_UNITS:
        raise ValueError(
            f"Unsupported length unit {unit!r}. Allowed units are 'cm' or 'mm'."
        )

    h_total = parse_finite_float(total_height, "total_height")
    d_val = parse_finite_float(diameter, "diameter")
    h_base = parse_finite_float(base_height, "base_height")

    # Source: Ph. Eur. 2.9.36 (Powder Flow, Angle of Repose) — physical dimensions h_total > 0, d > 0, h_base >= 0
    if h_total <= 0.0:
        raise ValueError(
            f"Total height must be > 0 {unit_clean} (received {h_total} {unit_clean})."
        )
    if d_val <= 0.0:
        raise ValueError(
            f"Cone base diameter must be > 0 {unit_clean} (received {d_val} {unit_clean})."
        )
    if h_base < 0.0:
        raise ValueError(
            f"Base/pedestal height cannot be negative (received {h_base} {unit_clean})."
        )

    # Source: HHU Düsseldorf Pharmaceutical Technology Lab Protocol (PD01) / Ph. Eur. 2.9.36 —
    # Net powder cone height h_cone = h_total - h_base must be strictly positive
    cone_height = h_total - h_base
    if cone_height <= 0.0:
        raise ValueError(
            f"Net cone height (total_height - base_height = {h_total} - {h_base} = {cone_height} "
            f"{unit_clean}) must be > 0 {unit_clean}."
        )

    # Source: SI Brochure, 9th ed. (BIPM): 1 cm = 10 mm
    if unit_clean == "cm":
        cone_height_cm = cone_height
        cone_height_mm = cone_height * MM_PER_CM  # Source: SI Brochure, 9th ed. (BIPM): 1 cm = 10 mm
        diameter_cm = d_val
        diameter_mm = d_val * MM_PER_CM  # Source: SI Brochure, 9th ed. (BIPM): 1 cm = 10 mm
    else:
        cone_height_mm = cone_height
        cone_height_cm = cone_height / MM_PER_CM  # Source: SI Brochure, 9th ed. (BIPM): 1 cm = 10 mm
        diameter_mm = d_val
        diameter_cm = d_val / MM_PER_CM  # Source: SI Brochure, 9th ed. (BIPM): 1 cm = 10 mm

    # Source: IEEE 754-2019 — verify unit-converted dimensions remain finite (guard against * 10 overflow)
    if not all(math.isfinite(v) for v in (cone_height_cm, cone_height_mm, diameter_cm, diameter_mm)):
        raise ValueError(
            "Angle of repose dimensions overflow finite floating-point limits during unit conversion."
        )

    return {
        "unit": unit_clean,
        "total_height": h_total,
        "base_height": h_base,
        "cone_height": cone_height,
        "diameter": d_val,
        "cone_height_cm": cone_height_cm,
        "cone_height_mm": cone_height_mm,
        "diameter_cm": diameter_cm,
        "diameter_mm": diameter_mm,
    }


def validate_flow_rate_data(mass_g, flow_time_s) -> tuple[float, float]:
    """
    Validates inputs for mass flow rate through an orifice (Q = m / t).

    Returns (mass_g, flow_time_s) as finite floats.
    """
    mass_clean = parse_finite_float(mass_g, "mass_g")
    time_clean = parse_finite_float(flow_time_s, "flow_time_s")

    # Source: Ph. Eur. 2.9.16 (Flowability) — sample mass m > 0 g and efflux time t > 0 s
    if mass_clean <= 0.0:
        raise ValueError(f"Powder mass 'mass_g' must be > 0 g (received {mass_clean} g).")
    if time_clean <= 0.0:
        raise ValueError(f"Flow time 'flow_time_s' must be > 0 s (received {time_clean} s).")

    return (mass_clean, time_clean)


def load_sieve_csv(csv_text: str) -> tuple[list[float], list[float], list[float]]:
    """
    Parses and validates a sieve analysis CSV string containing columns
    'mesh_um', 'empty_g', and 'full_g'.

    Returns (mesh_list, empty_list, full_list).
    """
    df = _read_csv_table(csv_text, REQUIRED_SIEVE_COLUMNS)

    mesh_list: list[float] = []
    empty_list: list[float] = []
    full_list: list[float] = []

    for line_no, (_, row) in enumerate(df.iterrows(), start=2):
        mesh_list.append(parse_finite_float(row["mesh_um"], f"CSV line {line_no} 'mesh_um'"))
        empty_list.append(parse_finite_float(row["empty_g"], f"CSV line {line_no} 'empty_g'"))
        full_list.append(parse_finite_float(row["full_g"], f"CSV line {line_no} 'full_g'"))

    # Source: Ph. Eur. 2.9.38 / ISO 9276-1 — ensure upper boundary strictly exceeds coarsest sieve aperture without ULP collapse
    x_max_auto = max(DEFAULT_X_MAX_UM, mesh_list[0] * 1.25)
    if not math.isfinite(x_max_auto):
        raise ValueError(
            f"Coarsest mesh size ({mesh_list[0]} µm) is too large to compute a finite upper boundary x_max_um."
        )
    validated = validate_sieve_data(mesh_list, empty_list, full_list, x_max_um=x_max_auto)

    return (
        validated["mesh_sizes"],
        validated["empty_weights"],
        validated["full_weights"],
    )


def load_flowability_batch_csv(csv_text: str) -> list[dict]:
    """
    Parses and validates a batch flowability CSV string containing columns
    'sample_code', 'group_run', 'mass_g', 'v0_ml', 'vt_ml', 'total_height_mm',
    'h_base_mm', and 'diameter_mm'.

    Returns a list of validated record dictionaries.
    """
    df = _read_csv_table(csv_text, REQUIRED_FLOWABILITY_COLUMNS)

    records: list[dict] = []
    for line_no, (_, row) in enumerate(df.iterrows(), start=2):
        sample_code = str(row["sample_code"]).strip()
        group_run = str(row["group_run"]).strip()

        if not sample_code:
            raise ValueError(f"CSV line {line_no}: 'sample_code' cannot be empty.")
        if not group_run:
            raise ValueError(f"CSV line {line_no}: 'group_run' cannot be empty.")

        try:
            mass_g, v0_ml, vt_ml = validate_carr_hausner_data(
                row["mass_g"],
                row["v0_ml"],
                row["vt_ml"],
            )
            # Source: HHU Düsseldorf PD01 flowability dataset records cone geometry in mm
            aor = validate_angle_of_repose_data(
                row["total_height_mm"],
                row["diameter_mm"],
                base_height=row["h_base_mm"],
                unit="mm",
            )
        except ValueError as exc:
            raise ValueError(
                f"CSV line {line_no} ({sample_code} / {group_run}): {exc}"
            ) from exc

        records.append({
            "sample_code": sample_code,
            "group_run": group_run,
            "mass_g": mass_g,
            "v0_ml": v0_ml,
            "vt_ml": vt_ml,
            "total_height_mm": aor["total_height"],
            "h_base_mm": aor["base_height"],
            "cone_height_mm": aor["cone_height_mm"],
            "diameter_mm": aor["diameter_mm"],
        })

    return records


# =============================================================================
# MODULE B: SIEVE ANALYSIS MATH, INTERPOLATION & MATPLOTLIB PLOTTING
# Sources: Ph. Eur. 2.9.38 (Particle-size distribution estimation by analytical
#          sieving); ISO 9276-1 (Representation of results of particle size
#          analysis - Part 1: Graphical representation); USP <429>
# =============================================================================

# Visualization color constants (Group 5 UI palette specification)
PLOT_COLOR = "#3b8a75"  # Primary teal color for Q3 curves and q3 histogram bars
D_COLOR = "#4b3f63"     # Muted dark purple for D10, D50, and D90 quantile markers


def _require_finite_float(value, name: str) -> float:
    """
    Validates that a scalar input is a finite real number (excluding booleans)
    and returns it as a Python float.
    """
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"{name} must be a numeric value (received {type(value).__name__}).")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite (received {value}).")
    return result


def _to_finite_float_list(values, name: str) -> list[float]:
    """
    Validates that a sequence contains only finite numeric values (no booleans,
    NaN, Inf, or non-numeric entries) and converts it to a list of floats.
    """
    if isinstance(values, (str, bytes)) or not hasattr(values, "__iter__"):
        raise ValueError(f"{name} must be a sequence of finite numeric values.")
    cleaned = []
    for idx, val in enumerate(values):
        cleaned.append(_require_finite_float(val, f"{name}[{idx}]"))
    return cleaned


def _unpack_validated_sieve(validated, fallback_x_max: float) -> tuple[list[float], list[float], list[float], list[float], float]:
    """
    Extracts cleaned sieve arrays and x_max_um from the return value of
    validate_sieve_data (supports dict or 5-tuple conventions from Module A).
    """
    if isinstance(validated, dict):
        mesh_raw = validated.get("mesh_sizes", validated.get("mesh_um", validated.get("mesh_clean")))
        empty_raw = validated.get("empty_weights", validated.get("empty_g", validated.get("empty_clean")))
        full_raw = validated.get("full_weights", validated.get("full_g", validated.get("full_clean")))
        retained_raw = validated.get("retained_masses", validated.get("retained_g", validated.get("retained_mass_g")))
        x_max_raw = validated.get("x_max_um", validated.get("x_max_clean", fallback_x_max))
        if any(item is None for item in (mesh_raw, empty_raw, full_raw, retained_raw, x_max_raw)):
            raise ValueError("Validated sieve payload is missing one or more required fields.")
    elif isinstance(validated, tuple) and len(validated) == 5:
        mesh_raw, empty_raw, full_raw, retained_raw, x_max_raw = validated
    else:
        raise TypeError(
            "validate_sieve_data must return a dict or a 5-tuple "
            "(mesh_clean, empty_clean, full_clean, retained_g, x_max_clean)."
        )

    mesh_clean = _to_finite_float_list(mesh_raw, "mesh_clean")
    empty_clean = _to_finite_float_list(empty_raw, "empty_clean")
    full_clean = _to_finite_float_list(full_raw, "full_clean")
    retained_g = _to_finite_float_list(retained_raw, "retained_g")
    x_max_clean = _require_finite_float(x_max_raw, "x_max_um")

    if not (len(mesh_clean) == len(empty_clean) == len(full_clean) == len(retained_g) >= 2):
        raise ValueError("Validated sieve arrays must have equal length and at least 2 entries.")

    return mesh_clean, empty_clean, full_clean, retained_g, x_max_clean


def interpolate_quantile(
    q_target: float,
    q3_asc: list[float],
    x_asc: list[float],
    mode: str = "linear",
) -> tuple[float, bool]:
    """
    Interpolates a particle size quantile x(q_target) from an ascending
    cumulative undersize distribution (q3_asc [%], x_asc [µm]).

    Returns:
        (interpolated_size_um, fell_back_to_linear)
    """
    if not isinstance(mode, str) or mode not in {"linear", "log"}:
        raise ValueError(f"Interpolation mode must be 'linear' or 'log' (received {mode!r}).")

    q_val = _require_finite_float(q_target, "q_target")
    # Source: ISO 9276-1; Ph. Eur. 2.9.38 (interior cumulative undersize quantile range 0% < Q3 < 100%)
    if not (0.0 < q_val < 100.0):
        raise ValueError(f"Target quantile q_target must satisfy 0 < q_target < 100 (received {q_val}).")

    q_vals = _to_finite_float_list(q3_asc, "q3_asc")
    x_vals = _to_finite_float_list(x_asc, "x_asc")

    if len(q_vals) != len(x_vals) or len(q_vals) < 2:
        raise ValueError("q3_asc and x_asc must have equal length and contain at least 2 points.")

    for i in range(len(q_vals) - 1):
        # Source: ISO 9276-1 (cumulative undersize Q3 is monotonically non-decreasing with particle size)
        if q_vals[i] > q_vals[i + 1]:
            raise ValueError(
                f"q3_asc must be monotonically non-decreasing (found {q_vals[i]} followed by {q_vals[i + 1]})."
            )
        # Source: ISO 9276-1 (sieve aperture grid x_asc must be strictly increasing)
        if x_vals[i] >= x_vals[i + 1]:
            raise ValueError(
                f"x_asc must be strictly increasing (found {x_vals[i]} followed by {x_vals[i + 1]})."
            )

    if not (q_vals[0] <= q_val <= q_vals[-1]):
        raise ValueError(
            f"Target quantile {q_val}% lies outside cumulative distribution range "
            f"[{q_vals[0]}%, {q_vals[-1]}%]."
        )

    # Locate bracketing interval with q_l <= q_target <= q_u and strictly positive step q_u > q_l
    # (skips flat zero-mass plateaus safely to prevent division by zero).
    for i in range(len(q_vals) - 1):
        q_l, q_u = q_vals[i], q_vals[i + 1]
        x_l, x_u = x_vals[i], x_vals[i + 1]

        if q_l <= q_val <= q_u and q_u > q_l:
            # Exact match at a distribution knot or interior plateau boundary returns the knot directly.
            if q_val == q_l:
                return float(x_l), False
            if q_val == q_u:
                return float(x_u), False

            # Normalized interpolation fraction within the bracketing cumulative interval
            # Source: ISO 9276-1; Ph. Eur. 2.9.38
            fraction = (q_val - q_l) / (q_u - q_l)

            if mode == "log":
                if x_l > 0.0 and x_u > 0.0:
                    # Source: ISO 9276-1; Ph. Eur. 2.9.38 (logarithmic interpolation on Q3 vs log10(x))
                    log_x = math.log10(x_l) + fraction * (math.log10(x_u) - math.log10(x_l))
                    return float(10.0 ** log_x), False
                # Explicit fallback to linear interpolation when lower aperture <= 0 µm (e.g. pan)
                # Source: ISO 9276-1; Ph. Eur. 2.9.38 (linear interpolation formula)
                linear_val = x_l + fraction * (x_u - x_l)
                return float(linear_val), True

            # Source: ISO 9276-1; Ph. Eur. 2.9.38 (linear interpolation formula)
            linear_val = x_l + fraction * (x_u - x_l)
            return float(linear_val), False

    raise ValueError(
        f"Cannot interpolate q_target={q_val}%: no strictly increasing Q3 interval brackets the target."
    )


def analyze_sieve_data(
    mesh_sizes,
    empty_weights,
    full_weights,
    x_max_um: float = 2000.0,
    interpolation_mode: str = "linear",
) -> tuple[pd.DataFrame, dict]:
    """
    Performs analytical sieve distribution calculations, quantile interpolation
    (D10, D50, D90), and span computation in accordance with Ph. Eur. 2.9.38
    and ISO 9276-1.
    """
    if not isinstance(interpolation_mode, str) or interpolation_mode not in {"linear", "log"}:
        raise ValueError(
            f"interpolation_mode must be 'linear' or 'log' (received {interpolation_mode!r})."
        )

    validated = validate_sieve_data(
        mesh_sizes, empty_weights, full_weights, x_max_um=x_max_um
    )
    mesh_clean, empty_clean, full_clean, retained_g, x_max_clean = _unpack_validated_sieve(
        validated, fallback_x_max=x_max_um
    )

    # Source: Ph. Eur. 2.9.38; ISO 9276-1 (total recovered powder mass across all sieves and pan)
    total_mass = float(sum(retained_g))
    if not math.isfinite(total_mass) or total_mass <= 0.0:
        raise ValueError("Total retained powder mass must be strictly positive (> 0 g).")

    # Source: ISO 9276-1; Ph. Eur. 2.9.38 (mass fraction Delta Q3_i [%] = (m_i / m_total) * 100)
    delta_q3_pct = [(m / total_mass) * 100.0 for m in retained_g]

    # Source: ISO 9276-1; Ph. Eur. 2.9.38 (cumulative undersize Q3(x_i) [%] = sum of finer fractions)
    # Clamped to [0.0, 100.0] to prevent IEEE-754 floating-point overshoot beyond 100.0%.
    q3_undersize_pct = [
        min(100.0, max(0.0, float(sum(delta_q3_pct[i + 1:]))))
        for i in range(len(mesh_clean))
    ]

    # Source: ISO 9276-1 (upper aperture limit of size class i is adjacent coarser sieve or assumed x_max)
    upper_bounds = [x_max_clean] + mesh_clean[:-1]
    # Source: ISO 9276-1 (particle size class width Delta x_i = x_upper,i - x_lower,i [µm])
    class_widths = [hi - lo for lo, hi in zip(mesh_clean, upper_bounds)]
    # Source: ISO 9276-1 (arithmetic class midpoint x_m,i = x_lower,i + (x_upper,i - x_lower,i) / 2 [µm]; overflow-safe)
    class_midpoints = [lo + (hi - lo) / 2.0 for lo, hi in zip(mesh_clean, upper_bounds)]
    # Source: ISO 9276-1 (mass density distribution q3(x_i) = Delta Q3_i / Delta x_i [%/µm])
    q3_density = [
        (dq / dx) if dx > 0.0 else 0.0
        for dq, dx in zip(delta_q3_pct, class_widths)
    ]

    n_stages = len(mesh_clean)
    sieve_labels = [f"Sieve {i + 1}" for i in range(n_stages - 1)] + ["Pan"]

    # Store unrounded finite floats in df_sieve so downstream plotting preserves exact small-mesh apertures
    df_sieve = pd.DataFrame({
        "Sieve": sieve_labels,
        "Mesh size [µm]": mesh_clean,
        "Empty sieve [g]": empty_clean,
        "Sieve + residue [g]": full_clean,
        "Retained mass [g]": retained_g,
        "dQ3 [%]": delta_q3_pct,
        "Q3 (undersize) [%]": q3_undersize_pct,
        "Class width dx [µm]": class_widths,
        "Midpoint [µm]": class_midpoints,
        "q3 [%/µm]": q3_density,
    })

    # Ascending (x, Q3) grid from pan (0 µm, 0%) to assumed maximum particle size (x_max_um, 100%)
    # Source: ISO 9276-1; Ph. Eur. 2.9.38
    x_asc = list(reversed(mesh_clean)) + [x_max_clean]
    q3_asc = list(reversed(q3_undersize_pct)) + [100.0]

    # Source: ISO 9276-1; Ph. Eur. 2.9.38; USP <429> (characteristic quantiles at Q3 = 10%, 50%, 90%)
    d10, fb_d10 = interpolate_quantile(10.0, q3_asc, x_asc, mode=interpolation_mode)
    d50, fb_d50 = interpolate_quantile(50.0, q3_asc, x_asc, mode=interpolation_mode)
    d90, fb_d90 = interpolate_quantile(90.0, q3_asc, x_asc, mode=interpolation_mode)

    fallback_quantiles = [
        label
        for label, fell_back in (("D10", fb_d10), ("D50", fb_d50), ("D90", fb_d90))
        if fell_back
    ]
    if fallback_quantiles:
        interpolation_note = (
            "Logarithmic interpolation fell back to linear interpolation for "
            f"{', '.join(fallback_quantiles)} because the lower bracketing aperture is <= 0 µm."
        )
    else:
        interpolation_note = None

    # Source: ISO 9276-1 (interdecile distribution width = D90 - D10 [µm])
    span_um = d90 - d10
    if not math.isfinite(d50) or d50 <= 0.0:
        raise ValueError(
            f"Median particle size D50 must be strictly positive (> 0 µm) to compute relative span (got {d50})."
        )
    # Source: ISO 9276-1; USP <429> (dimensionless relative span = (D90 - D10) / D50)
    relative_span = span_um / d50

    metrics = {
        "total_mass_g": round(total_mass, 4),
        "x_max_um": round(x_max_clean, 4),
        "interpolation_mode": interpolation_mode,
        "interpolation_note": interpolation_note,
        "D10_um": round(d10, 4),
        "D50_um": round(d50, 4),
        "D90_um": round(d90, 4),
        "span_um": round(span_um, 4),
        "relative_span": round(relative_span, 4),
    }

    return df_sieve, metrics


def new_plot(title: str, x_label: str, y_label: str, x_limits: tuple[float, float], y_top: float):
    """
    Creates a formatted Matplotlib figure and axes configured with a logarithmic
    particle size abscissa in accordance with ISO 9276-1 graphical conventions.
    """
    if not hasattr(x_limits, "__len__") or len(x_limits) != 2:
        raise ValueError("x_limits must be a 2-element sequence (x_min, x_max).")

    x_min = _require_finite_float(x_limits[0], "x_limits[0]")
    x_max = _require_finite_float(x_limits[1], "x_limits[1]")
    y_upper = _require_finite_float(y_top, "y_top")

    if not (math.isfinite(x_min) and math.isfinite(x_max)) or x_min <= 0.0 or x_max <= x_min:
        raise ValueError(
            f"Log-scale x_limits must satisfy 0 < x_min < x_max (received ({x_min}, {x_max}))."
        )
    if x_max > 1e100 or x_min < 1e-100:
        raise ValueError(
            f"Log-scale x_limits ({x_min}, {x_max}) exceed safe Matplotlib LogLocator range [1e-100, 1e100]."
        )
    if y_upper <= 0.0:
        raise ValueError(f"y_top must be strictly positive (received {y_upper}).")

    # Visualization layout parameter: 8.6 x 4.6 inch aspect ratio for embedded browser panels
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    # Source: ISO 9276-1; Ph. Eur. 2.9.38 (logarithmic particle size abscissa)
    ax.set_xscale("log")
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.0, y_upper)
    ax.set_xlabel(x_label, fontsize=11, fontweight="bold")
    ax.set_ylabel(y_label, fontsize=11, fontweight="bold")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.set_axisbelow(True)
    ax.set_title(title, fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig, ax


def overlap_area(box_a, box_b) -> float:
    """
    Computes the pixel overlap area between two Matplotlib Bbox objects.
    Returns 0.0 when the boxes do not intersect.
    """
    # Standard 2D axis-aligned rectangle intersection area
    width = min(float(box_a.x1), float(box_b.x1)) - max(float(box_a.x0), float(box_b.x0))
    height = min(float(box_a.y1), float(box_b.y1)) - max(float(box_a.y0), float(box_b.y0))
    return float(max(width, 0.0) * max(height, 0.0))


def mark_quantiles(ax, metrics: dict, quantile_lines: bool = False) -> None:
    """
    Overlays D10, D50, and D90 vertical dashed markers (and optional horizontal
    dotted lines at 10%, 50%, 90%) onto an axes object, placing text callouts
    at candidate positions that minimize overlap with plot lines, bars, other
    labels, and axes boundaries.
    """
    if not isinstance(metrics, dict):
        raise ValueError("metrics must be a dictionary containing D10_um, D50_um, and D90_um.")

    # Source: ISO 9276-1; Ph. Eur. 2.9.38; USP <429> (characteristic quantiles at 10%, 50%, 90%)
    raw_quantiles = [
        (10.0, "D10", metrics.get("D10_um")),
        (50.0, "D50", metrics.get("D50_um")),
        (90.0, "D90", metrics.get("D90_um")),
    ]

    valid_quantiles = []
    for q_pct, name, raw_val in raw_quantiles:
        if quantile_lines:
            ax.axhline(q_pct, color="gray", linestyle=":", alpha=0.6)
        if raw_val is None or isinstance(raw_val, bool):
            continue
        try:
            val = float(raw_val)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(val) or val <= 0.0:
            # Skip non-positive or non-finite values safely on logarithmic abscissa
            continue
        ax.axvline(val, color=D_COLOR, linestyle="--", linewidth=1.6)
        valid_quantiles.append((q_pct, name, val))

    if not valid_quantiles:
        return

    # Force a canvas draw before querying renderer and pixel bounding boxes
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    axes_box = ax.get_window_extent(renderer)

    # Sample plotted line segments in display pixel coordinates for collision scoring
    sampled_segments = []
    for line in ax.lines:
        xydata = np.asarray(line.get_xydata(), dtype=float)
        if xydata.ndim != 2 or xydata.shape[0] < 2:
            continue
        pixels = line.get_transform().transform(xydata)
        for start, end in zip(pixels[:-1], pixels[1:]):
            if np.all(np.isfinite(start)) and np.all(np.isfinite(end)):
                # Visualization collision sampling: 60 points per line segment
                sampled_segments.append(np.linspace(start, end, 60))

    if sampled_segments:
        line_points = np.vstack(sampled_segments)
    else:
        line_points = np.empty((0, 2), dtype=float)

    bar_boxes = [bar.get_window_extent(renderer) for bar in ax.patches]
    label_boxes = []

    y_bottom, y_top = ax.get_ylim()
    y_span = y_top - y_bottom
    if not math.isfinite(y_span) or y_span <= 0.0:
        return

    for q_pct, name, val in valid_quantiles:
        # Visualization layout fractions along vertical quantile guide from top to bottom
        heights = [
            y_bottom + y_span * frac
            for frac in (0.90, 0.75, 0.60, 0.45, 0.30, 0.15)
        ]
        if quantile_lines and y_bottom < q_pct < y_top:
            heights.insert(0, q_pct)

        label = ax.annotate(
            f"{name} = {val:.1f} µm",
            xy=(val, heights[0]),
            xytext=(0, 0),
            textcoords="offset points",
            fontsize=9,
            fontweight="bold",
            color=D_COLOR,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.5),
        )

        def move_label(y_pos: float, dx_pt: float, ha: str, dy_pt: float, va: str):
            label.xy = (val, y_pos)
            label.set(position=(dx_pt, dy_pt), ha=ha, va=va)
            # Visualization padding: 3 px safety margin around text bounding box
            return label.get_window_extent(renderer).padded(3.0)

        def penalty(spot) -> float:
            box = move_label(*spot)
            score = 0.0
            if line_points.size > 0:
                score += float(box.count_contains(line_points))
            box_area = float(box.width * box.height)
            if box_area > 0.0 and bar_boxes:
                # Visualization overlap weight: 50x fractional overlap with histogram bars
                score += 50.0 * sum(overlap_area(box, bar_box) for bar_box in bar_boxes) / box_area
            if label_boxes and box.count_overlaps(label_boxes) > 0:
                # Visualization overlap weight: strong penalty (1000) for overlapping another label
                score += 1000.0
            if not (axes_box.contains(box.x0, box.y0) and axes_box.contains(box.x1, box.y1)):
                # Visualization boundary weight: strong penalty (1000) for clipping outside axes
                score += 1000.0
            return score

        # Candidate anchor offsets in points (right/left and above/below each height)
        spots = [
            (y_pos, dx_pt, ha, dy_pt, va)
            for y_pos in heights
            for dx_pt, ha in ((6.0, "left"), (-6.0, "right"))
            for dy_pt, va in ((6.0, "bottom"), (-6.0, "top"))
        ]
        best_spot = min(spots, key=penalty)
        label_boxes.append(move_label(*best_spot))


def plot_q3_cumulative(df_sieve: pd.DataFrame, metrics: dict):
    """
    Plots the cumulative undersize particle size distribution Q3 [%] against
    logarithmic sieve aperture x [µm] according to Ph. Eur. 2.9.38 and ISO 9276-1.
    """
    required_cols = {"Mesh size [µm]", "Q3 (undersize) [%]"}
    if not isinstance(df_sieve, pd.DataFrame) or not required_cols.issubset(df_sieve.columns):
        raise ValueError(f"df_sieve must be a DataFrame containing columns {sorted(required_cols)}.")

    # Exclude pan (0 µm) on logarithmic abscissa and sort ascending by aperture
    # Source: ISO 9276-1 (log10(0) is undefined; plot measured positive sieve apertures ascending)
    df_valid = df_sieve[df_sieve["Mesh size [µm]"] > 0.0].sort_values("Mesh size [µm]")
    if df_valid.empty:
        raise ValueError("No positive sieve mesh sizes available to plot cumulative distribution Q3.")

    x_sizes = df_valid["Mesh size [µm]"].to_numpy(dtype=float)
    q3_cum_vals = df_valid["Q3 (undersize) [%]"].to_numpy(dtype=float)

    x_min = float(x_sizes.min())
    x_max = float(x_sizes.max())
    x_max_um = _require_finite_float(metrics.get("x_max_um", x_max), "metrics['x_max_um']")

    # Visualization axis margins: default 50–2500 µm window expanded with 0.6x / 1.25x padding
    x_limits = (min(50.0, x_min * 0.6), max(2500.0, max(x_max, x_max_um) * 1.25))
    # Visualization y-axis ceiling: 105% provides headroom above the 100% cumulative plateau
    y_top = 105.0

    fig, ax = new_plot(
        title=r"Cumulative Size Distribution $Q_3$ (Ph. Eur. 2.9.38)",
        x_label=r"Particle Size Aperture ($x$) [µm] (log scale)",
        y_label=r"Cumulative Distribution ($Q_3$) [%]",
        x_limits=x_limits,
        y_top=y_top,
    )

    ax.plot(
        x_sizes,
        q3_cum_vals,
        color=PLOT_COLOR,
        marker="o",
        linewidth=2.4,
        markersize=6.5,
    )

    # Dashed extension to assumed upper particle size boundary (x_max_um, 100%)
    # Source: ISO 9276-1 (closure of cumulative undersize curve at Q3(x_max) = 100%)
    if x_max_um > x_max:
        ax.plot(
            [x_max, x_max_um],
            [float(q3_cum_vals[-1]), 100.0],
            color=PLOT_COLOR,
            linestyle="--",
            linewidth=1.8,
            marker="o",
            markersize=5.0,
            alpha=0.75,
        )

    mark_quantiles(ax, metrics, quantile_lines=True)
    return fig


def plot_q3_density(df_sieve: pd.DataFrame, metrics: dict):
    """
    Plots the mass density particle size distribution q3 [%/µm] as a
    class-width histogram against logarithmic particle size x [µm] according
    to ISO 9276-1.
    """
    required_cols = {"Mesh size [µm]", "Class width dx [µm]", "q3 [%/µm]"}
    if not isinstance(df_sieve, pd.DataFrame) or not required_cols.issubset(df_sieve.columns):
        raise ValueError(f"df_sieve must be a DataFrame containing columns {sorted(required_cols)}.")
    if df_sieve.empty:
        raise ValueError("df_sieve is empty; cannot plot density distribution q3.")

    lower = df_sieve["Mesh size [µm]"].to_numpy(dtype=float)
    # Source: ISO 9276-1 (upper interval boundary = lower aperture + class width dx)
    upper = lower + df_sieve["Class width dx [µm]"].to_numpy(dtype=float)
    q3_vals = df_sieve["q3 [%/µm]"].to_numpy(dtype=float)

    positive_lower = lower[lower > 0.0]
    if positive_lower.size == 0:
        raise ValueError("At least one positive mesh size is required for logarithmic density plotting.")

    # Visualization axis margins: pan class lower edge clipped at x_left on log abscissa
    x_left = min(50.0, float(positive_lower.min()) * 0.6)
    lower_plot = np.maximum(lower, x_left)
    x_right = max(2500.0, float(upper.max()) * 1.25)

    q3_max = float(q3_vals.max()) if q3_vals.size > 0 else 0.0
    # Visualization y-axis headroom: 1.25x peak density (fallback to 1.0 if all densities are 0)
    y_top = q3_max * 1.25 if (math.isfinite(q3_max) and q3_max > 0.0) else 1.0

    fig, ax = new_plot(
        title=r"Density Size Distribution $q_3$ (ISO 9276-1)",
        x_label=r"Particle Size ($x$) [µm] (log scale)",
        y_label=r"Density Distribution ($q_3$) [%/µm]",
        x_limits=(x_left, x_right),
        y_top=y_top,
    )

    ax.bar(
        lower_plot,
        q3_vals,
        width=np.maximum(upper - lower_plot, 0.0),
        align="edge",
        color=PLOT_COLOR,
        edgecolor="#1f2933",
        alpha=0.88,
    )

    mark_quantiles(ax, metrics, quantile_lines=False)
    return fig


def figure_to_base64_png(fig, dpi: int = 150) -> str:
    """
    Encodes a Matplotlib figure as a base64 ASCII PNG string in memory and
    always closes the figure to release resources.
    """
    dpi_val = _require_finite_float(dpi, "dpi")
    if dpi_val <= 0.0:
        plt.close(fig)
        raise ValueError(f"dpi must be strictly positive (received {dpi}).")

    buffer = io.BytesIO()
    try:
        # Visualization export resolution: default 150 dpi raster PNG
        fig.savefig(buffer, format="png", dpi=int(round(dpi_val)))
        return base64.b64encode(buffer.getvalue()).decode("ascii")
    finally:
        plt.close(fig)


# =============================================================================
# MODULE C: POWDER DENSITY & FLOWABILITY MATH
# Sources: Ph. Eur. 2.9.34, Ph. Eur. 2.9.36 (Tables 2.9.36.-1 & 2.9.36.-2),
#          Ph. Eur. 2.9.16, ISO 3534-1 / Ph. Eur. 5.3, SI Brochure (BIPM)
# =============================================================================

# =============================================================================
# CLASSIFICATION FUNCTIONS
# =============================================================================

def classify_carr_index(carr_index: float) -> str:
    """
    Maps Carr's Compressibility Index (%) to the Ph. Eur. 2.9.36 Table 2.9.36.-2
    flowability classification scale.
    """
    ci = parse_finite_float(carr_index, "Carr's Index [%]")
    # Compressibility Index CI = 100 * (V_0 - V_f) / V_0 is bounded to [0, 100) for V_0 >= V_f > 0
    # Source: Ph. Eur. 2.9.36 (Powder flow)
    if ci < 0.0 or ci >= 100.0:
        raise ValueError(
            f"Carr's Index [%] must be in [0, 100) (received {ci})."
        )

    if ci <= 10.0:    # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (<= 10%: Excellent)
        return "Excellent"
    elif ci <= 15.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (11-15%: Good)
        return "Good"
    elif ci <= 20.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (16-20%: Fair)
        return "Fair"
    elif ci <= 25.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (21-25%: Passable)
        return "Passable"
    elif ci <= 31.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (26-31%: Poor)
        return "Poor"
    elif ci <= 37.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (32-37%: Very poor)
        return "Very poor"
    else:             # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (> 37%: Very, very poor)
        return "Very, very poor"


def classify_hausner_ratio(hausner_ratio: float) -> str:
    """
    Maps the Hausner Ratio (-) to the Ph. Eur. 2.9.36 Table 2.9.36.-2
    flowability classification scale.
    """
    hr = parse_finite_float(hausner_ratio, "Hausner Ratio [-]")
    # Hausner Ratio HR = V_0 / V_f >= 1.0 for V_0 >= V_f > 0
    # Source: Ph. Eur. 2.9.36 (Powder flow)
    if hr < 1.0:
        raise ValueError(
            f"Hausner Ratio [-] must be >= 1.0 (received {hr})."
        )

    if hr <= 1.11:    # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (1.00-1.11: Excellent)
        return "Excellent"
    elif hr <= 1.18:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (1.12-1.18: Good)
        return "Good"
    elif hr <= 1.25:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (1.19-1.25: Fair)
        return "Fair"
    elif hr <= 1.34:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (1.26-1.34: Passable)
        return "Passable"
    elif hr <= 1.45:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (1.35-1.45: Poor)
        return "Poor"
    elif hr <= 1.59:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (1.46-1.59: Very poor)
        return "Very poor"
    else:             # Source: Ph. Eur. 2.9.36, Table 2.9.36.-2 (> 1.60: Very, very poor)
        return "Very, very poor"


def classify_angle_of_repose(alpha_deg: float) -> str:
    """
    Maps the Angle of Repose (degrees) to the Ph. Eur. 2.9.36 Table 2.9.36.-1
    flowability classification scale.
    """
    alpha = parse_finite_float(alpha_deg, "Angle of Repose [deg]")
    # Angle of repose alpha = arctan(2*h / d) is strictly within (0, 90) degrees for h > 0, d > 0
    # Source: Ph. Eur. 2.9.36 (Powder flow)
    if alpha <= 0.0 or alpha >= 90.0:
        raise ValueError(
            f"Angle of Repose [deg] must be strictly between 0 and 90 degrees (received {alpha})."
        )

    if alpha < 25.0:    # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (lower bound is 25 deg; < 25 deg is below table range)
        return "Excellent (< 25°; below Table 2.9.36.-1 range)"
    elif alpha <= 30.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (25-30 deg: Excellent)
        return "Excellent"
    elif alpha <= 35.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (31-35 deg: Good)
        return "Good"
    elif alpha <= 40.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (36-40 deg: Fair - aid not needed)
        return "Fair (aid not needed)"
    elif alpha <= 45.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (41-45 deg: Passable - may hang up)
        return "Passable (may hang up)"
    elif alpha <= 55.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (46-55 deg: Poor - must agitate, vibrate)
        return "Poor (must agitate, vibrate)"
    elif alpha <= 65.0:  # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (56-65 deg: Very poor)
        return "Very poor"
    else:                # Source: Ph. Eur. 2.9.36, Table 2.9.36.-1 (> 66 deg: Very, very poor)
        return "Very, very poor"


# =============================================================================
# SINGLE-RUN ANALYSIS FUNCTIONS
# =============================================================================

def analyze_carr_hausner(mass_g, v0_ml, vf_ml) -> dict:
    """
    Calculates Bulk Density, Tapped Density, Carr's Compressibility Index,
    Hausner Ratio, and Ph. Eur. 2.9.36 classifications.
    """
    validated = validate_carr_hausner_data(mass_g, v0_ml, vf_ml)
    if isinstance(validated, dict):
        m = validated["mass_g"]
        v0 = validated["v0_ml"]
        vf = validated["vf_ml"]
    elif isinstance(validated, (tuple, list)) and len(validated) == 3:
        m, v0, vf = validated
    else:
        m = parse_finite_float(mass_g, "Powder mass [g]")
        v0 = parse_finite_float(v0_ml, "Poured volume V0 [mL]")
        vf = parse_finite_float(vf_ml, "Tapped volume Vf [mL]")

    # Bulk (poured) density: rho_bulk = m / V_0 [g/mL]
    # Source: Ph. Eur. 2.9.34 (Bulk density and tapped density of powders)
    rho_bulk = m / v0

    # Tapped density: rho_tapped = m / V_f [g/mL]
    # Source: Ph. Eur. 2.9.34 (Bulk density and tapped density of powders)
    rho_tapped = m / vf

    # Compressibility Index (Carr's Index): CI = 100 * (V_0 - V_f) / V_0 [%]
    # Source: Ph. Eur. 2.9.36 (Powder flow)
    carr_index = ((v0 - vf) / v0) * 100.0

    # Hausner Ratio: HR = V_0 / V_f [-]
    # Source: Ph. Eur. 2.9.36 (Powder flow)
    hausner_ratio = v0 / vf

    for label, val in (
        ("Bulk Density [g/mL]", rho_bulk),
        ("Tapped Density [g/mL]", rho_tapped),
        ("Carr's Index [%]", carr_index),
        ("Hausner Ratio [-]", hausner_ratio),
    ):
        if not math.isfinite(val):
            raise ValueError(f"Numerical overflow while computing {label}.")

    ci_rating = classify_carr_index(carr_index)
    hr_rating = classify_hausner_ratio(hausner_ratio)

    return {
        "bulk_density_g_ml": round(rho_bulk, 4),
        "tapped_density_g_ml": round(rho_tapped, 4),
        "carr_index_pct": round(carr_index, 4),
        "hausner_ratio": round(hausner_ratio, 4),
        "ci_rating": ci_rating,
        "hr_rating": hr_rating,
    }


def analyze_angle_of_repose(total_height, diameter, base_height=0.0, unit="cm") -> dict:
    """
    Calculates the Angle of Repose from cone total height, optional base height,
    and cone base diameter, and classifies flowability per Ph. Eur. 2.9.36.
    """
    valid = validate_angle_of_repose_data(
        total_height,
        diameter,
        base_height=base_height,
        unit=unit,
    )

    cone_height = valid["cone_height"]
    diam = valid["diameter"]
    unit_norm = valid["unit"]

    # Angle of repose: tan(alpha) = h / r = (2 * h) / d  =>  alpha = atan2(h, d / 2) [deg] (overflow-safe)
    # Source: Ph. Eur. 2.9.36 (Powder flow: Angle of repose)
    alpha_deg = math.degrees(math.atan2(cone_height, diam / 2.0))
    aor_rating = classify_angle_of_repose(alpha_deg)

    if "cone_height_cm" in valid and "cone_height_mm" in valid:
        cone_height_cm = valid["cone_height_cm"]
        cone_height_mm = valid["cone_height_mm"]
    elif unit_norm == "cm":
        cone_height_cm = cone_height
        cone_height_mm = cone_height * MM_PER_CM  # Source: SI Brochure (BIPM): 1 cm = 10 mm
    elif unit_norm == "mm":
        cone_height_mm = cone_height
        cone_height_cm = cone_height / MM_PER_CM  # Source: SI Brochure (BIPM): 1 mm = 0.1 cm
    else:
        raise ValueError(f"Unsupported length unit '{unit_norm}'. Expected 'cm' or 'mm'.")

    return {
        "unit": unit_norm,
        "total_height": round(valid["total_height"], 4),
        "base_height": round(valid["base_height"], 4),
        "cone_height": round(cone_height, 4),
        "cone_height_cm": round(cone_height_cm, 4),
        "cone_height_mm": round(cone_height_mm, 4),
        "diameter": round(diam, 4),
        "angle_of_repose_deg": round(alpha_deg, 4),
        "aor_rating": aor_rating,
    }


def analyze_flow_rate(method: str, measured_value=None, mass_g=None, flow_time_s=None) -> dict:
    """
    Calculates orifice mass flow rate [g/s] and normalized flow time per 100 g [s/100 g]
    without premature input rounding.
    """
    if not isinstance(method, str):
        raise ValueError(
            "Flow rate method must be a string ('fixed_mass', 'fixed_time', or 'custom')."
        )

    method_norm = method.strip().lower()

    if method_norm == "fixed_mass":
        if measured_value is None:
            raise ValueError(
                "Method 'fixed_mass' requires 'measured_value' (measured flow time in seconds)."
            )
        m = ORIFICE_FIXED_MASS_G  # Source: Ph. Eur. 2.9.16 (Flowability): fixed test mass = 100 g
        t = measured_value
    elif method_norm == "fixed_time":
        if measured_value is None:
            raise ValueError(
                "Method 'fixed_time' requires 'measured_value' (measured powder mass in grams)."
            )
        m = measured_value
        t = ORIFICE_FIXED_TIME_S  # Source: Ph. Eur. 2.9.36 (Flow through an orifice): fixed-time method = 10 s
    elif method_norm == "custom":
        if mass_g is None or flow_time_s is None:
            raise ValueError(
                "Method 'custom' requires both 'mass_g' and 'flow_time_s'."
            )
        m = mass_g
        t = flow_time_s
    else:
        raise ValueError(
            f"Unknown flow rate method '{method}'. Allowed methods: 'fixed_mass', 'fixed_time', 'custom'."
        )

    validated = validate_flow_rate_data(m, t)
    if isinstance(validated, dict):
        m_val = validated["mass_g"]
        t_val = validated["flow_time_s"]
    elif isinstance(validated, (tuple, list)) and len(validated) == 2:
        m_val, t_val = validated
    else:
        m_val = parse_finite_float(m, "Powder mass [g]")
        t_val = parse_finite_float(t, "Flow time [s]")

    # Mass flow rate through an orifice: Q = m / t [g/s]
    # Source: Ph. Eur. 2.9.36 (Flow through an orifice)
    flow_rate_g_s = m_val / t_val

    # Normalized flow time per 100 g of powder: t_100g = (100 * t) / m [s / 100 g]
    # Source: Ph. Eur. 2.9.16 (Flowability: time required for 100 g sample to flow through funnel)
    time_per_100g_s = (ORIFICE_FIXED_MASS_G * t_val) / m_val

    if not (math.isfinite(flow_rate_g_s) and math.isfinite(time_per_100g_s)):
        raise ValueError("Numerical overflow while computing flow rate metrics.")

    return {
        "method": method_norm,
        "mass_g": round(m_val, 4),
        "flow_time_s": round(t_val, 4),
        "flow_rate_g_s": round(flow_rate_g_s, 4),
        "time_per_100g_s": round(time_per_100g_s, 4),
    }


# =============================================================================
# BATCH FLOWABILITY ANALYSIS
# =============================================================================

def _compute_mean_sd(values):
    """
    Computes arithmetic mean and unbiased sample standard deviation (ddof=1)
    using scale-safe normalization to prevent IEEE 754 overflow on extreme values.
    Returns (raw_mean, rounded_mean, rounded_sd_or_none) so n < 2 never yields NaN.
    """
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        raise ValueError("Cannot compute summary statistics on empty or non-finite values.")

    # Scale-safe normalization prevents overflow in sum(x_i) and sum((x_i - x_bar)^2) near float64 limits
    # Source: IEEE 754-2019 / ISO 3534-1
    scale = float(np.max(np.abs(arr)))
    if scale == 0.0:
        raw_mean = 0.0
        raw_sd = 0.0 if arr.size >= 2 else None
    else:
        scaled = arr / scale
        # Arithmetic sample mean: x_bar = (1 / n) * sum(x_i)
        # Source: ISO 3534-1 (Statistics - Vocabulary and symbols)
        raw_mean = float(np.mean(scaled)) * scale
        # Unbiased sample standard deviation (Bessel's correction, ddof=1):
        # s = sqrt(sum((x_i - x_bar)^2) / (n - 1)) for n >= 2; undefined (None) for n < 2
        # Source: ISO 3534-1 / Ph. Eur. 5.3 (Statistical analysis of results of biological assays and tests)
        if arr.size >= 2:
            raw_sd = float(np.std(scaled, ddof=1)) * scale
        else:
            raw_sd = None

    if not math.isfinite(raw_mean) or (raw_sd is not None and not math.isfinite(raw_sd)):
        raise ValueError("Numerical overflow while computing batch summary statistics.")

    rounded_mean = round(raw_mean, 4)
    rounded_sd = round(raw_sd, 4) if raw_sd is not None else None
    return raw_mean, rounded_mean, rounded_sd


def analyze_flowability_batch(csv_text: str) -> dict:
    """
    Runs batch Carr/Hausner and Angle of Repose calculations on a validated
    flowability CSV dataset and produces per-sample_code summary statistics
    (mean and unbiased sample SD with ddof=1 computed from unrounded row metrics,
    preserving insertion order).
    """
    raw_rows = load_flowability_batch_csv(csv_text)
    if hasattr(raw_rows, "to_dict"):
        raw_rows = raw_rows.to_dict(orient="records")

    if not raw_rows:
        raise ValueError("Flowability batch CSV contains no data rows.")

    processed_rows = []
    grouped_raw = {}

    for idx, row in enumerate(raw_rows, start=1):
        sample_code = str(row["sample_code"]).strip()
        group_run = str(row.get("group_run", "")).strip()
        vf_raw = row["vt_ml"] if "vt_ml" in row else row["vf_ml"]

        try:
            ch = analyze_carr_hausner(
                mass_g=row["mass_g"],
                v0_ml=row["v0_ml"],
                vf_ml=vf_raw,
            )
            aor = analyze_angle_of_repose(
                total_height=row["total_height_mm"],
                diameter=row["diameter_mm"],
                base_height=row["h_base_mm"],
                unit="mm",
            )
            m_val = parse_finite_float(row["mass_g"], "mass_g")
            v0_val = parse_finite_float(row["v0_ml"], "v0_ml")
            vt_val = parse_finite_float(vf_raw, "vt_ml")
            h_total_val = parse_finite_float(row["total_height_mm"], "total_height_mm")
            h_base_val = parse_finite_float(row["h_base_mm"], "h_base_mm")
            diam_val = parse_finite_float(row["diameter_mm"], "diameter_mm")

            # Unrounded per-row metrics for unbiased group Mean & SD aggregation
            # Source: Ph. Eur. 2.9.34 (rho_bulk = m / V_0, rho_tapped = m / V_f)
            raw_bd = m_val / v0_val
            raw_td = m_val / vt_val
            # Source: Ph. Eur. 2.9.36 (CI = 100 * (V_0 - V_f) / V_0, HR = V_0 / V_f)
            raw_ci = ((v0_val - vt_val) / v0_val) * 100.0
            raw_hr = v0_val / vt_val
            # Source: Ph. Eur. 2.9.36 (alpha = arctan(2 * h_cone / d) = atan2(h_cone, d / 2))
            raw_aor = math.degrees(math.atan2(h_total_val - h_base_val, diam_val / 2.0))
        except ValueError as exc:
            raise ValueError(
                f"Batch row {idx} (sample '{sample_code}'): {exc}"
            ) from exc

        row_record = {
            "sample_code": sample_code,
            "group_run": group_run,
            "mass_g": round(m_val, 4),
            "v0_ml": round(v0_val, 4),
            "vt_ml": round(vt_val, 4),
            "total_height_mm": aor["total_height"],
            "h_base_mm": aor["base_height"],
            "cone_height_mm": aor["cone_height_mm"],
            "diameter_mm": aor["diameter"],
            "bulk_density_g_ml": ch["bulk_density_g_ml"],
            "tapped_density_g_ml": ch["tapped_density_g_ml"],
            "carr_index_pct": ch["carr_index_pct"],
            "hausner_ratio": ch["hausner_ratio"],
            "ci_rating": ch["ci_rating"],
            "hr_rating": ch["hr_rating"],
            "angle_of_repose_deg": aor["angle_of_repose_deg"],
            "aor_rating": aor["aor_rating"],
        }
        processed_rows.append(row_record)

        if sample_code not in grouped_raw:
            grouped_raw[sample_code] = []
        grouped_raw[sample_code].append({
            "bd": raw_bd,
            "td": raw_td,
            "ci": raw_ci,
            "hr": raw_hr,
            "aor": raw_aor,
        })

    summary_rows = []
    for sample_code, items in grouped_raw.items():
        n = len(items)
        _, bd_mean, bd_sd = _compute_mean_sd([r["bd"] for r in items])
        _, td_mean, td_sd = _compute_mean_sd([r["td"] for r in items])
        ci_raw_mean, ci_mean, ci_sd = _compute_mean_sd([r["ci"] for r in items])
        hr_raw_mean, hr_mean, hr_sd = _compute_mean_sd([r["hr"] for r in items])
        aor_raw_mean, aor_mean, aor_sd = _compute_mean_sd([r["aor"] for r in items])

        summary_rows.append({
            "sample_code": sample_code,
            "n": n,
            "bulk_density_mean_g_ml": bd_mean,
            "bulk_density_sd_g_ml": bd_sd,
            "tapped_density_mean_g_ml": td_mean,
            "tapped_density_sd_g_ml": td_sd,
            "carr_index_mean_pct": ci_mean,
            "carr_index_sd_pct": ci_sd,
            "hausner_ratio_mean": hr_mean,
            "hausner_ratio_sd": hr_sd,
            "angle_of_repose_mean_deg": aor_mean,
            "angle_of_repose_sd_deg": aor_sd,
            "ci_rating": classify_carr_index(ci_raw_mean),
            "hr_rating": classify_hausner_ratio(hr_raw_mean),
            "aor_rating": classify_angle_of_repose(aor_raw_mean),
        })

    return {
        "row_count": len(processed_rows),
        "group_count": len(summary_rows),
        "rows": processed_rows,
        "summary": summary_rows,
    }


# =============================================================================
# MODULE D: SELF-CONTAINED EMBEDDED WEB UI & LOCAL OFFLINE HTTP SERVER
# Serves an offline laboratory web dashboard over 127.0.0.1 with zero external
# files, zero external URLs/CDNs, and strict anti-crash input validation.
# =============================================================================

# =============================================================================
# EMBEDDED WEB ASSETS (100% OFFLINE — ZERO EXTERNAL FILES OR CDNs)
# =============================================================================

HTML_INDEX: str = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Group 5 — Powder Characterization Tool</title>
  <link rel="stylesheet" href="/style.css">
</head>
<body>
  <header class="app-header">
    <div class="header-inner">
      <div class="header-brand">
        <h1>Group 5 &mdash; Powder Characterization Tool</h1>
        <p class="subtitle">Statistics &amp; DoE &middot; HHU D&uuml;sseldorf &middot; Ph. Eur. 2.9.16 / 2.9.34 / 2.9.36 / 2.9.38 &amp; ISO 9276-1</p>
      </div>
      <div class="header-actions">
        <span class="status-pill" aria-label="Application status">
          <span class="status-dot" aria-hidden="true"></span>
          Offline Local Mode
        </span>
        <button type="button" id="btn-quit" class="btn btn-danger-outline">Quit Server</button>
      </div>
    </div>
  </header>

  <main class="container">
    <nav class="tab-bar" role="tablist" aria-label="Powder characterization modules">
      <button type="button" class="tab-btn active" role="tab" id="tab-btn-sieve" aria-selected="true" aria-controls="panel-sieve" data-tab="sieve">
        Sieve Analysis
      </button>
      <button type="button" class="tab-btn" role="tab" id="tab-btn-carr" aria-selected="false" aria-controls="panel-carr" data-tab="carr">
        Carr &amp; Hausner
      </button>
      <button type="button" class="tab-btn" role="tab" id="tab-btn-angle" aria-selected="false" aria-controls="panel-angle" data-tab="angle">
        Angle of Repose
      </button>
      <button type="button" class="tab-btn" role="tab" id="tab-btn-flow" aria-selected="false" aria-controls="panel-flow" data-tab="flow">
        Orifice Flow Rate
      </button>
      <button type="button" class="tab-btn" role="tab" id="tab-btn-batch" aria-selected="false" aria-controls="panel-batch" data-tab="batch">
        Batch Flowability CSV
      </button>
    </nav>

    <!-- =================================================================== -->
    <!-- TAB 1: SIEVE ANALYSIS (Ph. Eur. 2.9.38 / ISO 9276-1)                -->
    <!-- =================================================================== -->
    <section id="panel-sieve" class="tab-panel active" role="tabpanel" aria-labelledby="tab-btn-sieve">
      <div class="card">
        <div class="card-header">
          <div>
            <h2>Particle Size Distribution by Sieving</h2>
            <p class="card-desc">Enter sieve weights from largest mesh aperture (top sieve) down to the pan (0 &micro;m), or import a CSV file.</p>
          </div>
        </div>

        <div id="error-sieve" class="alert alert-error hidden" role="alert"></div>

        <div class="toolbar">
          <label class="btn btn-secondary file-label">
            Upload Sieve CSV
            <input type="file" id="sieve-csv-input" accept=".csv,text/csv" hidden>
          </label>
          <button type="button" id="btn-load-pd01" class="btn btn-secondary">Load PD01 Sample</button>
          <button type="button" id="btn-load-ref" class="btn btn-secondary">Load Reference Sample</button>
          <button type="button" id="btn-add-sieve-row" class="btn btn-accent-outline">+ Add Sieve Row</button>
        </div>

        <div class="form-grid form-grid-2">
          <div class="form-field">
            <label for="sieve-xmax">Upper Size Limit x<sub>max</sub> [&micro;m]</label>
            <input type="text" id="sieve-xmax" inputmode="decimal" value="2000" placeholder="e.g. 2000">
            <span class="field-hint">Assumed upper bound above top sieve (must exceed largest mesh size).</span>
          </div>
          <div class="form-field">
            <label for="sieve-interp">Quantile Interpolation Mode</label>
            <select id="sieve-interp">
              <option value="linear" selected>Linear (Ph. Eur. 2.9.38)</option>
              <option value="log">Logarithmic (ISO 9276-1)</option>
            </select>
            <span class="field-hint">Interpolation applied on cumulative undersize curve Q<sub>3</sub>(x).</span>
          </div>
        </div>

        <div class="table-scroll">
          <table class="data-table editable-table" id="sieve-input-table">
            <thead>
              <tr>
                <th>Stage</th>
                <th>Mesh Size [&micro;m]</th>
                <th>Empty Sieve [g]</th>
                <th>Sieve + Residue [g]</th>
                <th class="col-action">Action</th>
              </tr>
            </thead>
            <tbody id="sieve-tbody"></tbody>
          </table>
        </div>

        <div class="action-bar">
          <button type="button" id="btn-calc-sieve" class="btn btn-primary">Calculate Sieve Analysis</button>
        </div>
      </div>

      <div id="results-sieve" class="results-block hidden">
        <div id="sieve-interp-note" class="alert alert-info hidden" role="status"></div>

        <div class="kpi-grid kpi-grid-6">
          <div class="kpi-card">
            <span class="kpi-label">D10 [&micro;m]</span>
            <span class="kpi-value" id="sieve-ans-d10">&mdash;</span>
          </div>
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">D50 (Median) [&micro;m]</span>
            <span class="kpi-value" id="sieve-ans-d50">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">D90 [&micro;m]</span>
            <span class="kpi-value" id="sieve-ans-d90">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">Span (D90 &minus; D10) [&micro;m]</span>
            <span class="kpi-value" id="sieve-ans-span">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">Relative Span [&minus;]</span>
            <span class="kpi-value" id="sieve-ans-relspan">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">Total Retained Mass [g]</span>
            <span class="kpi-value" id="sieve-ans-mass">&mdash;</span>
          </div>
        </div>

        <div class="charts-grid">
          <div class="card chart-card">
            <h3>Cumulative Undersize Distribution Q<sub>3</sub></h3>
            <img id="sieve-chart-q3" class="chart-img" alt="Cumulative undersize particle size distribution Q3 plot">
          </div>
          <div class="card chart-card">
            <h3>Density Size Distribution q<sub>3</sub></h3>
            <img id="sieve-chart-density" class="chart-img" alt="Density particle size distribution q3 histogram">
          </div>
        </div>

        <div class="card">
          <h3 id="sieve-summary-title">Sieve Analysis Summary Table</h3>
          <div class="table-scroll">
            <table class="data-table" id="sieve-summary-table">
              <thead id="sieve-summary-head"></thead>
              <tbody id="sieve-summary-body"></tbody>
            </table>
          </div>
        </div>
      </div>

      <details class="card media-card" open>
        <summary class="media-summary">Experiment Guide &mdash; Analytical Sieving</summary>
        <div class="media-body" id="media-slot-sieve">
          <p class="media-empty">Loading experiment media&hellip;</p>
        </div>
      </details>
    </section>

    <!-- =================================================================== -->
    <!-- TAB 2: CARR'S INDEX & HAUSNER RATIO (Ph. Eur. 2.9.34 / 2.9.36)      -->
    <!-- =================================================================== -->
    <section id="panel-carr" class="tab-panel" role="tabpanel" aria-labelledby="tab-btn-carr">
      <div class="card">
        <div class="card-header">
          <div>
            <h2>Bulk / Tapped Density, Carr&rsquo;s Index &amp; Hausner Ratio</h2>
            <p class="card-desc">Evaluates powder compressibility and flowability according to Ph. Eur. 2.9.34 and Ph. Eur. 2.9.36 Table 2.9.36.-2.</p>
          </div>
        </div>

        <div id="error-carr" class="alert alert-error hidden" role="alert"></div>

        <div class="form-grid form-grid-3">
          <div class="form-field">
            <label for="carr-mass">Powder Mass m [g]</label>
            <input type="text" id="carr-mass" inputmode="decimal" placeholder="e.g. 100.4">
          </div>
          <div class="form-field">
            <label for="carr-v0">Poured Volume V<sub>0</sub> [mL]</label>
            <input type="text" id="carr-v0" inputmode="decimal" placeholder="e.g. 143.0">
          </div>
          <div class="form-field">
            <label for="carr-vf">Tapped Volume V<sub>f</sub> [mL]</label>
            <input type="text" id="carr-vf" inputmode="decimal" placeholder="e.g. 121.0">
          </div>
        </div>

        <div class="action-bar">
          <button type="button" id="btn-calc-carr" class="btn btn-primary">Calculate</button>
          <button type="button" id="btn-example-carr" class="btn btn-secondary">Load Example</button>
        </div>
      </div>

      <div id="results-carr" class="results-block hidden">
        <div class="kpi-grid kpi-grid-3">
          <div class="kpi-card">
            <span class="kpi-label">Bulk Density &rho;<sub>bulk</sub> [g/mL]</span>
            <span class="kpi-value" id="carr-ans-bulk">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">Tapped Density &rho;<sub>tapped</sub> [g/mL]</span>
            <span class="kpi-value" id="carr-ans-tapped">&mdash;</span>
          </div>
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">Carr&rsquo;s Index CI [%]</span>
            <span class="kpi-value" id="carr-ans-ci">&mdash;</span>
          </div>
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">Hausner Ratio HR [&minus;]</span>
            <span class="kpi-value" id="carr-ans-hr">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">CI Rating (Ph. Eur. 2.9.36)</span>
            <span class="kpi-badge-wrap" id="carr-ans-ci-rating">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">HR Rating (Ph. Eur. 2.9.36)</span>
            <span class="kpi-badge-wrap" id="carr-ans-hr-rating">&mdash;</span>
          </div>
        </div>
      </div>

      <div class="card ref-card">
        <h3>Reference: Ph. Eur. 2.9.36 &mdash; Table 2.9.36.-2 (Scale of Flowability)</h3>
        <div class="table-scroll">
          <table class="data-table ref-table">
            <thead>
              <tr>
                <th>Flow Character</th>
                <th>Compressibility Index (CI) [%]</th>
                <th>Hausner Ratio (HR) [&minus;]</th>
              </tr>
            </thead>
            <tbody>
              <tr><td><span class="badge badge-good">Excellent</span></td><td>&le; 10 (1&ndash;10)</td><td>1.00&ndash;1.11</td></tr>
              <tr><td><span class="badge badge-good">Good</span></td><td>11&ndash;15</td><td>1.12&ndash;1.18</td></tr>
              <tr><td><span class="badge badge-warn">Fair</span></td><td>16&ndash;20</td><td>1.19&ndash;1.25</td></tr>
              <tr><td><span class="badge badge-warn">Passable</span></td><td>21&ndash;25</td><td>1.26&ndash;1.34</td></tr>
              <tr><td><span class="badge badge-poor">Poor</span></td><td>26&ndash;31</td><td>1.35&ndash;1.45</td></tr>
              <tr><td><span class="badge badge-poor">Very poor</span></td><td>32&ndash;37</td><td>1.46&ndash;1.59</td></tr>
              <tr><td><span class="badge badge-poor">Very, very poor</span></td><td>&gt; 38 (&gt; 37)</td><td>&gt; 1.60 (&gt; 1.59)</td></tr>
            </tbody>
          </table>
        </div>
      </div>

      <details class="card media-card" open>
        <summary class="media-summary">Experiment Guide &mdash; Bulk &amp; Tapped Density</summary>
        <div class="media-body" id="media-slot-carr">
          <p class="media-empty">Loading experiment media&hellip;</p>
        </div>
      </details>
    </section>

    <!-- =================================================================== -->
    <!-- TAB 3: ANGLE OF REPOSE (Ph. Eur. 2.9.36 Table 2.9.36.-1)            -->
    <!-- =================================================================== -->
    <section id="panel-angle" class="tab-panel" role="tabpanel" aria-labelledby="tab-btn-angle">
      <div class="card">
        <div class="card-header">
          <div>
            <h2>Angle of Repose</h2>
            <p class="card-desc">Computes effective cone height h<sub>cone</sub> = h<sub>total</sub> &minus; h<sub>base</sub> and &alpha; = arctan(2 &middot; h<sub>cone</sub> / d) per Ph. Eur. 2.9.36.</p>
          </div>
        </div>

        <div id="error-angle" class="alert alert-error hidden" role="alert"></div>

        <div class="form-grid form-grid-4">
          <div class="form-field">
            <label for="angle-unit">Length Unit</label>
            <select id="angle-unit">
              <option value="cm" selected>Centimeters (cm)</option>
              <option value="mm">Millimeters (mm)</option>
            </select>
          </div>
          <div class="form-field">
            <label for="angle-htotal">Total Height h<sub>total</sub> [<span class="angle-unit-lbl">cm</span>]</label>
            <input type="text" id="angle-htotal" inputmode="decimal" placeholder="e.g. 5.4">
          </div>
          <div class="form-field">
            <label for="angle-hbase">Pedestal Base Height h<sub>base</sub> [<span class="angle-unit-lbl">cm</span>]</label>
            <input type="text" id="angle-hbase" inputmode="decimal" value="0" placeholder="0">
          </div>
          <div class="form-field">
            <label for="angle-diameter">Base Diameter d [<span class="angle-unit-lbl">cm</span>]</label>
            <input type="text" id="angle-diameter" inputmode="decimal" placeholder="e.g. 10.0">
          </div>
        </div>

        <div class="action-bar">
          <button type="button" id="btn-calc-angle" class="btn btn-primary">Calculate</button>
          <button type="button" id="btn-example-angle" class="btn btn-secondary">Load Example</button>
        </div>
      </div>

      <div id="results-angle" class="results-block hidden">
        <div class="kpi-grid kpi-grid-3">
          <div class="kpi-card">
            <span class="kpi-label">Effective Cone Height h<sub>cone</sub></span>
            <span class="kpi-value" id="angle-ans-hcone">&mdash;</span>
          </div>
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">Angle of Repose &alpha; [&deg;]</span>
            <span class="kpi-value" id="angle-ans-deg">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">Flow Property (Table 2.9.36.-1)</span>
            <span class="kpi-badge-wrap" id="angle-ans-rating">&mdash;</span>
          </div>
        </div>
      </div>

      <div class="card ref-card">
        <h3>Reference: Ph. Eur. 2.9.36 &mdash; Table 2.9.36.-1 (Flow Properties &amp; Angle of Repose)</h3>
        <div class="table-scroll">
          <table class="data-table ref-table">
            <thead>
              <tr>
                <th>Flow Property</th>
                <th>Angle of Repose &alpha; [&deg;]</th>
              </tr>
            </thead>
            <tbody>
              <tr><td><span class="badge badge-good">Excellent</span></td><td>25&ndash;30 (&le; 30)</td></tr>
              <tr><td><span class="badge badge-good">Good</span></td><td>31&ndash;35</td></tr>
              <tr><td><span class="badge badge-warn">Fair (aid not needed)</span></td><td>36&ndash;40</td></tr>
              <tr><td><span class="badge badge-warn">Passable (may hang up)</span></td><td>41&ndash;45</td></tr>
              <tr><td><span class="badge badge-poor">Poor (must agitate, vibrate)</span></td><td>46&ndash;55</td></tr>
              <tr><td><span class="badge badge-poor">Very poor</span></td><td>56&ndash;65</td></tr>
              <tr><td><span class="badge badge-poor">Very, very poor</span></td><td>&gt; 66 (&gt; 65)</td></tr>
            </tbody>
          </table>
        </div>
      </div>

      <details class="card media-card" open>
        <summary class="media-summary">Experiment Guide &mdash; Angle of Repose</summary>
        <div class="media-body" id="media-slot-angle">
          <p class="media-empty">Loading experiment media&hellip;</p>
        </div>
      </details>
    </section>

    <!-- =================================================================== -->
    <!-- TAB 4: ORIFICE FLOW RATE (Ph. Eur. 2.9.16 / 2.9.36)                 -->
    <!-- =================================================================== -->
    <section id="panel-flow" class="tab-panel" role="tabpanel" aria-labelledby="tab-btn-flow">
      <div class="card">
        <div class="card-header">
          <div>
            <h2>Flow Rate Through an Orifice</h2>
            <p class="card-desc">Computes mass flow rate (g/s) and standardized flow time per 100 g (s/100g) according to Ph. Eur. 2.9.16 and 2.9.36.</p>
          </div>
        </div>

        <div id="error-flow" class="alert alert-error hidden" role="alert"></div>

        <div class="form-grid form-grid-3">
          <div class="form-field">
            <label for="flow-method">Measurement Method</label>
            <select id="flow-method">
              <option value="fixed_mass" selected>Fixed Mass (100 g) &mdash; measure time t</option>
              <option value="fixed_time">Fixed Time (10 s) &mdash; measure mass m</option>
              <option value="custom">Custom (enter both m &amp; t)</option>
            </select>
          </div>
          <div class="form-field hidden" id="flow-mass-group">
            <label for="flow-mass">Powder Mass m [g]</label>
            <input type="text" id="flow-mass" inputmode="decimal" placeholder="e.g. 85.0">
            <span class="field-hint" id="flow-mass-hint">Mass flowing through the orifice.</span>
          </div>
          <div class="form-field" id="flow-time-group">
            <label for="flow-time">Flow Time t [s]</label>
            <input type="text" id="flow-time" inputmode="decimal" placeholder="e.g. 12.4">
            <span class="field-hint" id="flow-time-hint">Fixed test mass = 100.0 g (Ph. Eur. 2.9.16).</span>
          </div>
        </div>

        <div class="action-bar">
          <button type="button" id="btn-calc-flow" class="btn btn-primary">Calculate</button>
          <button type="button" id="btn-example-flow" class="btn btn-secondary">Load Example</button>
        </div>
      </div>

      <div id="results-flow" class="results-block hidden">
        <div class="kpi-grid kpi-grid-4">
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">Flow Rate [g/s]</span>
            <span class="kpi-value" id="flow-ans-rate">&mdash;</span>
          </div>
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">Flow Time per 100 g [s/100g]</span>
            <span class="kpi-value" id="flow-ans-t100">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">Mass Used [g]</span>
            <span class="kpi-value" id="flow-ans-mass">&mdash;</span>
          </div>
          <div class="kpi-card">
            <span class="kpi-label">Time Used [s]</span>
            <span class="kpi-value" id="flow-ans-time">&mdash;</span>
          </div>
        </div>
      </div>

      <details class="card media-card" open>
        <summary class="media-summary">Experiment Guide &mdash; Flow Through an Orifice</summary>
        <div class="media-body" id="media-slot-flow">
          <p class="media-empty">Loading experiment media&hellip;</p>
        </div>
      </details>
    </section>

    <!-- =================================================================== -->
    <!-- TAB 5: BATCH FLOWABILITY CSV                                        -->
    <!-- =================================================================== -->
    <section id="panel-batch" class="tab-panel" role="tabpanel" aria-labelledby="tab-btn-batch">
      <div class="card">
        <div class="card-header">
          <div>
            <h2>Batch Flowability CSV Analyzer</h2>
            <p class="card-desc">Upload a multi-run flowability CSV or load the embedded HHU course dataset for automated per-run and grouped Mean &plusmn; SD analysis.</p>
          </div>
        </div>

        <div id="error-batch" class="alert alert-error hidden" role="alert"></div>

        <div class="toolbar">
          <label class="btn btn-secondary file-label">
            Upload Batch CSV
            <input type="file" id="batch-csv-input" accept=".csv,text/csv" hidden>
          </label>
          <button type="button" id="btn-load-course-csv" class="btn btn-secondary">Load Course Dataset</button>
          <button type="button" id="btn-analyze-batch" class="btn btn-primary">Analyze Batch CSV</button>
        </div>

        <div class="form-field">
          <label for="batch-csv-text">CSV Content Preview / Input</label>
          <textarea id="batch-csv-text" rows="7" spellcheck="false" placeholder="sample_code,group_run,mass_g,v0_ml,vt_ml,total_height_mm,h_base_mm,diameter_mm"></textarea>
          <span class="field-hint">Upload a .csv file, click &ldquo;Load Course Dataset&rdquo;, or paste CSV rows directly.</span>
        </div>
      </div>

      <div id="results-batch" class="results-block hidden">
        <div class="kpi-grid kpi-grid-2">
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">Total Runs Analyzed</span>
            <span class="kpi-value" id="batch-ans-runs">&mdash;</span>
          </div>
          <div class="kpi-card kpi-highlight">
            <span class="kpi-label">Unique Formulations</span>
            <span class="kpi-value" id="batch-ans-formulations">&mdash;</span>
          </div>
        </div>

        <div class="card">
          <h3>Grouped Formulation Summary (Mean &plusmn; SD)</h3>
          <div class="table-scroll">
            <table class="data-table" id="batch-summary-table">
              <thead id="batch-summary-head"></thead>
              <tbody id="batch-summary-body"></tbody>
            </table>
          </div>
        </div>

        <div class="card">
          <h3>Individual Run Results</h3>
          <div class="table-scroll">
            <table class="data-table" id="batch-runs-table">
              <thead id="batch-runs-head"></thead>
              <tbody id="batch-runs-body"></tbody>
            </table>
          </div>
        </div>
      </div>

      <details class="card media-card">
        <summary class="media-summary">Lab Media Gallery</summary>
        <div class="media-body" id="media-slot-batch">
          <p class="media-empty">Loading media gallery&hellip;</p>
        </div>
      </details>
    </section>
  </main>

  <footer class="app-footer">
    <p>Group 5 &mdash; Powder Characterization Tool &middot; Offline Local Execution (127.0.0.1) &middot; Ph. Eur. 2.9.16, 2.9.34, 2.9.36, 2.9.38 &amp; ISO 9276-1</p>
  </footer>

  <div id="quit-overlay" class="quit-overlay hidden" role="dialog" aria-modal="true" aria-labelledby="quit-title">
    <div class="quit-modal">
      <div class="quit-icon" aria-hidden="true">&#10003;</div>
      <h2 id="quit-title">Local Server Stopped</h2>
      <p>The offline Powder Characterization Tool server has shut down cleanly.</p>
      <p class="quit-sub">You may now close this browser tab safely.</p>
    </div>
  </div>

  <script src="/app.js"></script>
</body>
</html>
"""

CSS_STYLES: str = """/* Group 5 Powder Characterization Tool — Embedded Offline Stylesheet */
:root {
  --bg: #f4f7f6;
  --surface: #ffffff;
  --surface-alt: #f9fbfa;
  --border: #d8e2df;
  --border-strong: #bccbc6;
  --text: #1c2825;
  --text-muted: #53645f;
  --primary: #2c6e5d;
  --primary-hover: #215548;
  --primary-soft: #e7f2ef;
  --accent: #4b3f63;
  --accent-hover: #3a304d;
  --accent-soft: #efecf5;
  --danger: #a62626;
  --danger-bg: #fdecec;
  --danger-border: #f3bcbc;
  --info-bg: #eef4fa;
  --info-text: #1d4e7a;
  --info-border: #bdd6ee;
  --good-bg: #e3f5ec;
  --good-text: #15623e;
  --good-border: #9fd9bc;
  --warn-bg: #fef4dc;
  --warn-text: #8a5600;
  --warn-border: #f0cf80;
  --poor-bg: #fde8e8;
  --poor-text: #9b1c1c;
  --poor-border: #f2b4b4;
  --neutral-bg: #edf1f0;
  --neutral-text: #475853;
  --neutral-border: #cfd9d6;
  --radius: 10px;
  --shadow: 0 2px 10px rgba(28, 40, 37, 0.06);
}

*, *::before, *::after {
  box-sizing: border-box;
}

body {
  margin: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  background-color: var(--bg);
  color: var(--text);
  line-height: 1.5;
  -webkit-font-smoothing: antialiased;
}

.hidden {
  display: none !important;
}

/* Header */
.app-header {
  background: var(--surface);
  border-bottom: 1px solid var(--border);
  box-shadow: 0 1px 4px rgba(28, 40, 37, 0.04);
  position: sticky;
  top: 0;
  z-index: 20;
}

.header-inner {
  max-width: 1200px;
  margin: 0 auto;
  padding: 16px 24px;
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
  flex-wrap: wrap;
}

.header-brand h1 {
  margin: 0;
  font-size: 1.35rem;
  font-weight: 700;
  color: var(--primary);
  letter-spacing: -0.01em;
}

.subtitle {
  margin: 3px 0 0;
  font-size: 0.85rem;
  color: var(--text-muted);
}

.header-actions {
  display: flex;
  align-items: center;
  gap: 12px;
  flex-wrap: wrap;
}

.status-pill {
  display: inline-flex;
  align-items: center;
  gap: 7px;
  padding: 5px 12px;
  border-radius: 999px;
  background: var(--good-bg);
  color: var(--good-text);
  border: 1px solid var(--good-border);
  font-size: 0.8rem;
  font-weight: 600;
}

.status-dot {
  width: 8px;
  height: 8px;
  border-radius: 50%;
  background: #1d8f59;
  display: inline-block;
}

/* Main Container */
.container {
  max-width: 1200px;
  margin: 22px auto 40px;
  padding: 0 24px;
}

/* Tabs */
.tab-bar {
  display: flex;
  gap: 6px;
  flex-wrap: wrap;
  margin-bottom: 20px;
  border-bottom: 2px solid var(--border);
  padding-bottom: 0;
}

.tab-btn {
  appearance: none;
  background: transparent;
  border: none;
  border-bottom: 3px solid transparent;
  margin-bottom: -2px;
  padding: 11px 18px;
  font-size: 0.94rem;
  font-weight: 600;
  color: var(--text-muted);
  cursor: pointer;
  border-radius: 8px 8px 0 0;
  transition: background-color 0.15s ease, color 0.15s ease, border-color 0.15s ease;
}

.tab-btn:hover {
  background: var(--primary-soft);
  color: var(--primary);
}

.tab-btn:focus-visible {
  outline: 2px solid var(--accent);
  outline-offset: -2px;
}

.tab-btn.active {
  color: var(--primary);
  border-bottom-color: var(--primary);
  background: var(--surface);
}

.tab-panel {
  display: none;
}

.tab-panel.active {
  display: block;
}

/* Cards */
.card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 22px 24px;
  margin-bottom: 20px;
  box-shadow: var(--shadow);
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
  margin-bottom: 16px;
}

.card h2 {
  margin: 0;
  font-size: 1.15rem;
  font-weight: 700;
  color: var(--text);
}

.card h3 {
  margin: 0 0 12px;
  font-size: 1rem;
  font-weight: 700;
  color: var(--accent);
}

.card-desc {
  margin: 4px 0 0;
  font-size: 0.88rem;
  color: var(--text-muted);
}

/* Alerts */
.alert {
  padding: 12px 15px;
  border-radius: 8px;
  margin-bottom: 16px;
  font-size: 0.9rem;
  font-weight: 500;
  border: 1px solid transparent;
}

.alert-error {
  background: var(--danger-bg);
  color: var(--danger);
  border-color: var(--danger-border);
}

.alert-info {
  background: var(--info-bg);
  color: var(--info-text);
  border-color: var(--info-border);
}

/* Buttons & Toolbars */
.toolbar,
.action-bar {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  align-items: center;
}

.toolbar {
  margin-bottom: 18px;
}

.action-bar {
  margin-top: 18px;
}

.btn {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  gap: 6px;
  padding: 9px 16px;
  border-radius: 8px;
  font-size: 0.9rem;
  font-weight: 600;
  cursor: pointer;
  border: 1px solid transparent;
  transition: background-color 0.15s ease, border-color 0.15s ease, color 0.15s ease;
  font-family: inherit;
  text-decoration: none;
}

.btn:disabled {
  opacity: 0.6;
  cursor: not-allowed;
}

.btn-primary {
  background: var(--primary);
  color: #ffffff;
}

.btn-primary:hover:not(:disabled) {
  background: var(--primary-hover);
}

.btn-secondary {
  background: var(--surface-alt);
  color: var(--text);
  border-color: var(--border-strong);
}

.btn-secondary:hover:not(:disabled) {
  background: var(--primary-soft);
  border-color: var(--primary);
  color: var(--primary);
}

.btn-accent-outline {
  background: var(--accent-soft);
  color: var(--accent);
  border-color: #c9c0db;
}

.btn-accent-outline:hover:not(:disabled) {
  background: var(--accent);
  color: #ffffff;
}

.btn-danger-outline {
  background: #ffffff;
  color: var(--danger);
  border-color: var(--danger-border);
  padding: 7px 14px;
  font-size: 0.85rem;
}

.btn-danger-outline:hover:not(:disabled) {
  background: var(--danger);
  color: #ffffff;
  border-color: var(--danger);
}

.btn-row-remove {
  background: transparent;
  color: var(--danger);
  border: 1px solid var(--danger-border);
  border-radius: 6px;
  padding: 4px 10px;
  font-size: 0.8rem;
  font-weight: 600;
  cursor: pointer;
}

.btn-row-remove:hover {
  background: var(--danger-bg);
}

/* Forms */
.form-grid {
  display: grid;
  gap: 16px;
  margin-bottom: 16px;
}

.form-grid-2 {
  grid-template-columns: repeat(2, minmax(0, 1fr));
}

.form-grid-3 {
  grid-template-columns: repeat(3, minmax(0, 1fr));
}

.form-grid-4 {
  grid-template-columns: repeat(4, minmax(0, 1fr));
}

.form-field {
  display: flex;
  flex-direction: column;
  gap: 5px;
}

.form-field label {
  font-size: 0.86rem;
  font-weight: 600;
  color: var(--text);
}

.form-field input,
.form-field select,
.form-field textarea {
  width: 100%;
  padding: 9px 11px;
  border: 1px solid var(--border-strong);
  border-radius: 7px;
  font-size: 0.92rem;
  font-family: inherit;
  color: var(--text);
  background: #ffffff;
}

.form-field textarea {
  font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
  font-size: 0.83rem;
  resize: vertical;
}

.form-field input:focus,
.form-field select:focus,
.form-field textarea:focus {
  outline: none;
  border-color: var(--primary);
  box-shadow: 0 0 0 3px rgba(44, 110, 93, 0.14);
}

.field-hint {
  font-size: 0.78rem;
  color: var(--text-muted);
}

/* Tables */
.table-scroll {
  width: 100%;
  overflow-x: auto;
  border: 1px solid var(--border);
  border-radius: 8px;
}

.data-table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.87rem;
  background: var(--surface);
}

.data-table th,
.data-table td {
  padding: 9px 12px;
  border-bottom: 1px solid var(--border);
  text-align: left;
  white-space: nowrap;
}

.data-table th {
  background: var(--surface-alt);
  color: var(--text);
  font-weight: 700;
  font-size: 0.82rem;
  text-transform: none;
}

.data-table tbody tr:last-child td {
  border-bottom: none;
}

.data-table tbody tr:hover {
  background: #f5faf8;
}

.editable-table input {
  width: 100%;
  min-width: 95px;
  padding: 6px 8px;
  border: 1px solid var(--border-strong);
  border-radius: 6px;
  font-size: 0.88rem;
  font-family: inherit;
}

.editable-table input:focus {
  outline: none;
  border-color: var(--primary);
  box-shadow: 0 0 0 2px rgba(44, 110, 93, 0.15);
}

.col-action {
  width: 90px;
  text-align: center;
}

.stage-tag {
  font-weight: 600;
  color: var(--accent);
}

/* KPI Grid */
.results-block {
  margin-top: 8px;
}

.kpi-grid {
  display: grid;
  gap: 14px;
  margin-bottom: 20px;
}

.kpi-grid-2 {
  grid-template-columns: repeat(2, minmax(0, 1fr));
}

.kpi-grid-3 {
  grid-template-columns: repeat(3, minmax(0, 1fr));
}

.kpi-grid-4 {
  grid-template-columns: repeat(4, minmax(0, 1fr));
}

.kpi-grid-6 {
  grid-template-columns: repeat(6, minmax(0, 1fr));
}

.kpi-card {
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius);
  padding: 15px 16px;
  box-shadow: var(--shadow);
  display: flex;
  flex-direction: column;
  gap: 6px;
}

.kpi-card.kpi-highlight {
  border-top: 3px solid var(--primary);
}

.kpi-label {
  font-size: 0.78rem;
  font-weight: 600;
  color: var(--text-muted);
}

.kpi-value {
  font-size: 1.35rem;
  font-weight: 700;
  color: var(--text);
  letter-spacing: -0.01em;
}

.kpi-badge-wrap {
  margin-top: 2px;
}

/* Flowability Rating Badges */
.badge {
  display: inline-block;
  padding: 4px 10px;
  border-radius: 999px;
  font-size: 0.8rem;
  font-weight: 700;
  border: 1px solid transparent;
  line-height: 1.3;
}

.badge-good {
  background: var(--good-bg);
  color: var(--good-text);
  border-color: var(--good-border);
}

.badge-warn {
  background: var(--warn-bg);
  color: var(--warn-text);
  border-color: var(--warn-border);
}

.badge-poor {
  background: var(--poor-bg);
  color: var(--poor-text);
  border-color: var(--poor-border);
}

.badge-neutral {
  background: var(--neutral-bg);
  color: var(--neutral-text);
  border-color: var(--neutral-border);
}

/* Charts */
.charts-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 20px;
  margin-bottom: 20px;
}

.chart-card {
  margin-bottom: 0;
}

.chart-img {
  width: 100%;
  height: auto;
  display: block;
  border-radius: 6px;
  border: 1px solid var(--border);
  background: #ffffff;
}

/* Footer */
.app-footer {
  max-width: 1200px;
  margin: 0 auto 30px;
  padding: 0 24px;
  text-align: center;
  font-size: 0.8rem;
  color: var(--text-muted);
}

/* Quit Overlay */
.quit-overlay {
  position: fixed;
  inset: 0;
  background: rgba(28, 40, 37, 0.72);
  display: flex;
  align-items: center;
  justify-content: center;
  padding: 20px;
  z-index: 100;
}

.quit-modal {
  background: var(--surface);
  border-radius: 14px;
  padding: 32px 36px;
  max-width: 440px;
  width: 100%;
  text-align: center;
  box-shadow: 0 18px 40px rgba(0, 0, 0, 0.22);
  border: 1px solid var(--border);
}

.quit-icon {
  width: 52px;
  height: 52px;
  margin: 0 auto 14px;
  border-radius: 50%;
  background: var(--primary-soft);
  color: var(--primary);
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 1.5rem;
  font-weight: 700;
}

.quit-modal h2 {
  margin: 0 0 8px;
  color: var(--primary);
}

.quit-modal p {
  margin: 6px 0;
  color: var(--text);
}

.quit-sub {
  font-size: 0.88rem;
  color: var(--text-muted) !important;
}

/* Responsive Breakpoints */
@media (max-width: 1024px) {
  .kpi-grid-6 {
    grid-template-columns: repeat(3, minmax(0, 1fr));
  }
  .charts-grid {
    grid-template-columns: 1fr;
  }
}

@media (max-width: 768px) {
  .container,
  .header-inner {
    padding-left: 14px;
    padding-right: 14px;
  }
  .form-grid-2,
  .form-grid-3,
  .form-grid-4,
  .kpi-grid-2,
  .kpi-grid-3,
  .kpi-grid-4,
  .kpi-grid-6 {
    grid-template-columns: 1fr;
  }
  .tab-btn {
    padding: 9px 12px;
    font-size: 0.86rem;
  }
}

.media-card {
  margin-top: 16px;
}

.media-summary {
  cursor: pointer;
  font-weight: 700;
  font-size: 0.98rem;
  color: var(--primary-dark, #1f5346);
  list-style: none;
  display: flex;
  align-items: center;
  gap: 8px;
}

.media-summary::-webkit-details-marker {
  display: none;
}

.media-summary::before {
  content: "\\25B8";
  font-size: 0.85rem;
  color: var(--primary, #2c6e5d);
  transition: transform 0.15s ease;
}

.media-card[open] > .media-summary::before {
  transform: rotate(90deg);
}

.media-body {
  margin-top: 14px;
}

.media-grid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
  gap: 16px;
}

.media-figure {
  margin: 0;
  padding: 12px;
  border: 1px solid var(--border, #d9e2df);
  border-radius: 10px;
  background: var(--surface-alt, #f9fbfa);
  display: flex;
  flex-direction: column;
  gap: 8px;
}

.media-figure img,
.media-figure video {
  width: 100%;
  height: auto;
  max-height: 360px;
  object-fit: contain;
  border-radius: 6px;
  background: #ffffff;
}

.media-empty {
  margin: 0;
  font-size: 0.88rem;
  color: var(--text-muted, #526066);
}
"""

JS_APP: str = """(function () {
  "use strict";

  function $(id) {
    return document.getElementById(id);
  }

  function escapeHtml(value) {
    return String(value == null ? "" : value)
      .replace(/&/g, "&amp;")
      .replace(/</g, "&lt;")
      .replace(/>/g, "&gt;")
      .replace(/"/g, "&quot;")
      .replace(/'/g, "&#39;");
  }

  function showError(tabKey, message) {
    const el = $("error-" + tabKey);
    if (!el) return;
    el.textContent = message;
    el.classList.remove("hidden");
  }

  function clearError(tabKey) {
    const el = $("error-" + tabKey);
    if (!el) return;
    el.textContent = "";
    el.classList.add("hidden");
  }

  function ratingBadgeClass(ratingText) {
    const s = String(ratingText || "").toLowerCase();
    if (!s || s === "n/a" || s === "-") return "badge-neutral";
    if (s.includes("poor")) return "badge-poor";
    if (s.includes("excellent") || s.includes("good")) return "badge-good";
    if (s.includes("fair") || s.includes("passable")) return "badge-warn";
    return "badge-neutral";
  }

  function renderBadgeHtml(text) {
    const cls = ratingBadgeClass(text);
    return '<span class="badge ' + cls + '">' + escapeHtml(text) + "</span>";
  }

  function isRatingHeader(header) {
    const h = String(header || "").toLowerCase();
    return h.includes("rating") || h.includes("flowability") || h.includes("property") || h.includes("class");
  }

  function renderDataTable(headId, bodyId, headers, rows, badgeRatingCols) {
    const thead = $(headId);
    const tbody = $(bodyId);
    if (!thead || !tbody) return;

    const headCells = headers.map(function (h) {
      return "<th>" + escapeHtml(h) + "</th>";
    }).join("");
    thead.innerHTML = "<tr>" + headCells + "</tr>";

    const bodyRows = rows.map(function (row) {
      const cells = row.map(function (cell, idx) {
        const text = cell == null ? "" : String(cell);
        if (badgeRatingCols && isRatingHeader(headers[idx]) && text && text !== "N/A") {
          return "<td>" + renderBadgeHtml(text) + "</td>";
        }
        return "<td>" + escapeHtml(text) + "</td>";
      }).join("");
      return "<tr>" + cells + "</tr>";
    }).join("");
    tbody.innerHTML = bodyRows;
  }

  async function postJson(endpoint, payload) {
    const response = await fetch(endpoint, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload || {})
    });

    let data;
    try {
      data = await response.json();
    } catch (_err) {
      throw new Error("Server returned an unreadable response (HTTP " + response.status + ").");
    }

    if (!response.ok && (!data || !data.error)) {
      throw new Error("Request failed with HTTP " + response.status + ".");
    }
    if (!data || data.ok === false) {
      throw new Error((data && data.error) ? data.error : "Calculation failed.");
    }
    return data;
  }

  /* -----------------------------------------------------------------------
   * TAB SWITCHING
   * --------------------------------------------------------------------- */
  function initTabs() {
    const buttons = document.querySelectorAll(".tab-btn");
    const panels = document.querySelectorAll(".tab-panel");

    buttons.forEach(function (btn) {
      btn.addEventListener("click", function () {
        const target = btn.getAttribute("data-tab");
        buttons.forEach(function (b) {
          const active = b === btn;
          b.classList.toggle("active", active);
          b.setAttribute("aria-selected", active ? "true" : "false");
        });
        panels.forEach(function (p) {
          p.classList.toggle("active", p.id === "panel-" + target);
        });
      });
    });
  }

  /* -----------------------------------------------------------------------
   * TAB 1: SIEVE ANALYSIS
   * --------------------------------------------------------------------- */
  function updateSieveStageLabels() {
    const tbody = $("sieve-tbody");
    const rows = tbody.querySelectorAll("tr");
    rows.forEach(function (tr, idx) {
      const meshInput = tr.querySelector('input[data-col="mesh"]');
      const stageEl = tr.querySelector(".stage-tag");
      const meshVal = meshInput ? meshInput.value.trim().replace(",", ".") : "";
      const isPan = (idx === rows.length - 1 && meshVal !== "" && (meshVal === "0" || Number(meshVal) === 0));
      if (stageEl) {
        stageEl.textContent = isPan ? "Pan" : ("Sieve " + (idx + 1));
      }
    });
  }

  function createSieveRow(mesh, emptyW, fullW) {
    const tr = document.createElement("tr");
    tr.innerHTML =
      '<td><span class="stage-tag"></span></td>' +
      '<td><input type="text" inputmode="decimal" data-col="mesh" value="' + escapeHtml(mesh) + '" placeholder="e.g. 500"></td>' +
      '<td><input type="text" inputmode="decimal" data-col="empty" value="' + escapeHtml(emptyW) + '" placeholder="e.g. 320.0"></td>' +
      '<td><input type="text" inputmode="decimal" data-col="full" value="' + escapeHtml(fullW) + '" placeholder="e.g. 325.0"></td>' +
      '<td class="col-action"><button type="button" class="btn-row-remove" title="Remove row">Remove</button></td>';

    const meshInput = tr.querySelector('input[data-col="mesh"]');
    meshInput.addEventListener("input", updateSieveStageLabels);

    const removeBtn = tr.querySelector(".btn-row-remove");
    removeBtn.addEventListener("click", function () {
      const tbody = $("sieve-tbody");
      if (tbody.querySelectorAll("tr").length <= 2) {
        showError("sieve", "At least 2 rows (sieve + pan) are required.");
        return;
      }
      clearError("sieve");
      tr.remove();
      updateSieveStageLabels();
    });

    return tr;
  }

  function setSieveTableRows(rows) {
    const tbody = $("sieve-tbody");
    tbody.innerHTML = "";
    rows.forEach(function (row) {
      tbody.appendChild(createSieveRow(row[0], row[1], row[2]));
    });
    updateSieveStageLabels();
  }

  function addSieveRow() {
    clearError("sieve");
    const tbody = $("sieve-tbody");
    const newRow = createSieveRow("", "", "");
    const existing = tbody.querySelectorAll("tr");
    if (existing.length > 0) {
      const last = existing[existing.length - 1];
      const lastMesh = last.querySelector('input[data-col="mesh"]');
      if (lastMesh && lastMesh.value.trim().replace(",", ".") === "0") {
        tbody.insertBefore(newRow, last);
        updateSieveStageLabels();
        return;
      }
    }
    tbody.appendChild(newRow);
    updateSieveStageLabels();
  }

  function collectSieveRows() {
    const tbody = $("sieve-tbody");
    const rows = [];
    tbody.querySelectorAll("tr").forEach(function (tr) {
      const mesh = tr.querySelector('input[data-col="mesh"]').value.trim();
      const emptyW = tr.querySelector('input[data-col="empty"]').value.trim();
      const fullW = tr.querySelector('input[data-col="full"]').value.trim();
      rows.push([mesh, emptyW, fullW]);
    });
    return rows;
  }

  async function loadSieveSample(sampleKey) {
    clearError("sieve");
    try {
      const data = await postJson("/api/sample", { sample: sampleKey });
      if (!Array.isArray(data.rows) || data.rows.length === 0) {
        throw new Error("Sample returned no sieve rows.");
      }
      setSieveTableRows(data.rows);
    } catch (err) {
      showError("sieve", err.message || "Could not load sample data.");
    }
  }

  function handleSieveCsvUpload(event) {
    const file = event.target.files && event.target.files[0];
    event.target.value = "";
    if (!file) return;
    clearError("sieve");

    const reader = new FileReader();
    reader.onload = async function (e) {
      try {
        const text = String(e.target.result || "");
        const data = await postJson("/api/read_sieve_csv", { text: text });
        setSieveTableRows(data.rows);
      } catch (err) {
        showError("sieve", err.message || "Could not parse sieve CSV file.");
      }
    };
    reader.onerror = function () {
      showError("sieve", "Failed to read the selected CSV file.");
    };
    reader.readAsText(file);
  }

  async function calculateSieve() {
    clearError("sieve");
    const btn = $("btn-calc-sieve");
    btn.disabled = true;
    try {
      const payload = {
        rows: collectSieveRows(),
        x_max_um: $("sieve-xmax").value.trim(),
        interp_mode: $("sieve-interp").value
      };
      const data = await postJson("/api/sieve", payload);
      const ans = data.answers || {};

      $("sieve-ans-d10").textContent = ans.d10 != null ? ans.d10 : "—";
      $("sieve-ans-d50").textContent = ans.d50 != null ? ans.d50 : "—";
      $("sieve-ans-d90").textContent = ans.d90 != null ? ans.d90 : "—";
      $("sieve-ans-span").textContent = ans.span != null ? ans.span : "—";
      $("sieve-ans-relspan").textContent = ans.relative_span != null ? ans.relative_span : "—";
      $("sieve-ans-mass").textContent = ans.mass != null ? ans.mass : "—";

      const noteEl = $("sieve-interp-note");
      if (data.interpolation_note) {
        noteEl.textContent = data.interpolation_note;
        noteEl.classList.remove("hidden");
      } else {
        noteEl.textContent = "";
        noteEl.classList.add("hidden");
      }

      if (Array.isArray(data.charts) && data.charts.length >= 2) {
        $("sieve-chart-q3").src = "data:image/png;base64," + data.charts[0];
        $("sieve-chart-density").src = "data:image/png;base64," + data.charts[1];
      }

      $("sieve-summary-title").textContent = data.table_title || "Sieve Analysis Summary Table";
      renderDataTable(
        "sieve-summary-head",
        "sieve-summary-body",
        data.table_header || [],
        data.table_rows || [],
        false
      );

      $("results-sieve").classList.remove("hidden");
    } catch (err) {
      showError("sieve", err.message || "Sieve calculation failed.");
    } finally {
      btn.disabled = false;
    }
  }

  /* -----------------------------------------------------------------------
   * TAB 2: CARR'S INDEX & HAUSNER RATIO
   * --------------------------------------------------------------------- */
  function loadCarrExample() {
    clearError("carr");
    $("carr-mass").value = "100.4";
    $("carr-v0").value = "143.0";
    $("carr-vf").value = "121.0";
  }

  async function calculateCarr() {
    clearError("carr");
    const btn = $("btn-calc-carr");
    btn.disabled = true;
    try {
      const data = await postJson("/api/carr_hausner", {
        mass: $("carr-mass").value.trim(),
        v0: $("carr-v0").value.trim(),
        vf: $("carr-vf").value.trim()
      });
      const ans = data.answers || {};
      $("carr-ans-bulk").textContent = ans.bulk != null ? ans.bulk : "—";
      $("carr-ans-tapped").textContent = ans.tapped != null ? ans.tapped : "—";
      $("carr-ans-ci").textContent = ans.carr != null ? ans.carr : "—";
      $("carr-ans-hr").textContent = ans.hausner != null ? ans.hausner : "—";
      $("carr-ans-ci-rating").innerHTML = renderBadgeHtml(ans.ci_rating || ans.property || "N/A");
      $("carr-ans-hr-rating").innerHTML = renderBadgeHtml(ans.hr_rating || ans.property || "N/A");
      $("results-carr").classList.remove("hidden");
    } catch (err) {
      showError("carr", err.message || "Carr & Hausner calculation failed.");
    } finally {
      btn.disabled = false;
    }
  }

  /* -----------------------------------------------------------------------
   * TAB 3: ANGLE OF REPOSE
   * --------------------------------------------------------------------- */
  function updateAngleUnitLabels() {
    const unit = $("angle-unit").value;
    document.querySelectorAll(".angle-unit-lbl").forEach(function (span) {
      span.textContent = unit;
    });
  }

  function loadAngleExample() {
    clearError("angle");
    $("angle-unit").value = "mm";
    updateAngleUnitLabels();
    $("angle-htotal").value = "54.0";
    $("angle-hbase").value = "28.0";
    $("angle-diameter").value = "100.0";
  }

  async function calculateAngle() {
    clearError("angle");
    const btn = $("btn-calc-angle");
    btn.disabled = true;
    try {
      const unit = $("angle-unit").value;
      const data = await postJson("/api/angle", {
        unit: unit,
        h_total: $("angle-htotal").value.trim(),
        h_base: $("angle-hbase").value.trim(),
        diameter: $("angle-diameter").value.trim()
      });
      const ans = data.answers || {};
      const outUnit = ans.unit || unit;
      $("angle-ans-hcone").textContent = (ans.h_cone != null ? ans.h_cone : "—") + " " + outUnit;
      $("angle-ans-deg").textContent = ans.angle != null ? ans.angle : "—";
      $("angle-ans-rating").innerHTML = renderBadgeHtml(ans.property || "N/A");
      $("results-angle").classList.remove("hidden");
    } catch (err) {
      showError("angle", err.message || "Angle of Repose calculation failed.");
    } finally {
      btn.disabled = false;
    }
  }

  /* -----------------------------------------------------------------------
   * TAB 4: ORIFICE FLOW RATE
   * --------------------------------------------------------------------- */
  function syncFlowMethodInputs() {
    clearError("flow");
    const method = $("flow-method").value;
    const massGroup = $("flow-mass-group");
    const timeGroup = $("flow-time-group");
    const massHint = $("flow-mass-hint");
    const timeHint = $("flow-time-hint");

    if (method === "fixed_mass") {
      massGroup.classList.add("hidden");
      timeGroup.classList.remove("hidden");
      timeHint.textContent = "Fixed test mass = 100.0 g (Ph. Eur. 2.9.16).";
    } else if (method === "fixed_time") {
      massGroup.classList.remove("hidden");
      timeGroup.classList.add("hidden");
      massHint.textContent = "Fixed test time = 10.0 s (Ph. Eur. 2.9.36).";
    } else {
      massGroup.classList.remove("hidden");
      timeGroup.classList.remove("hidden");
      massHint.textContent = "Measured powder mass m [g].";
      timeHint.textContent = "Measured discharge time t [s].";
    }
  }

  function loadFlowExample() {
    clearError("flow");
    const method = $("flow-method").value;
    if (method === "fixed_mass") {
      $("flow-time").value = "12.5";
    } else if (method === "fixed_time") {
      $("flow-mass").value = "80.0";
    } else {
      $("flow-mass").value = "75.0";
      $("flow-time").value = "9.4";
    }
  }

  async function calculateFlowRate() {
    clearError("flow");
    const btn = $("btn-calc-flow");
    btn.disabled = true;
    try {
      const method = $("flow-method").value;
      const payload = { method: method };
      if (method === "fixed_mass") {
        payload.time = $("flow-time").value.trim();
      } else if (method === "fixed_time") {
        payload.mass = $("flow-mass").value.trim();
      } else {
        payload.mass = $("flow-mass").value.trim();
        payload.time = $("flow-time").value.trim();
      }

      const data = await postJson("/api/flow_rate", payload);
      const ans = data.answers || {};
      $("flow-ans-rate").textContent = ans.rate != null ? ans.rate : "—";
      $("flow-ans-t100").textContent = ans.time_per_100g != null ? ans.time_per_100g : "—";
      $("flow-ans-mass").textContent = ans.mass_used != null ? ans.mass_used : "—";
      $("flow-ans-time").textContent = ans.time_used != null ? ans.time_used : "—";
      $("results-flow").classList.remove("hidden");
    } catch (err) {
      showError("flow", err.message || "Orifice flow rate calculation failed.");
    } finally {
      btn.disabled = false;
    }
  }

  /* -----------------------------------------------------------------------
   * TAB 5: BATCH FLOWABILITY CSV
   * --------------------------------------------------------------------- */
  function handleBatchCsvUpload(event) {
    const file = event.target.files && event.target.files[0];
    event.target.value = "";
    if (!file) return;
    clearError("batch");

    const reader = new FileReader();
    reader.onload = function (e) {
      $("batch-csv-text").value = String(e.target.result || "");
    };
    reader.onerror = function () {
      showError("batch", "Failed to read the selected CSV file.");
    };
    reader.readAsText(file);
  }

  async function loadCourseBatchDataset() {
    clearError("batch");
    try {
      const data = await postJson("/api/sample", { sample: "flowability_batch" });
      $("batch-csv-text").value = data.csv_text || "";
    } catch (err) {
      showError("batch", err.message || "Could not load course flowability dataset.");
    }
  }

  async function analyzeBatchCsv() {
    clearError("batch");
    const btn = $("btn-analyze-batch");
    btn.disabled = true;
    try {
      const text = $("batch-csv-text").value;
      const data = await postJson("/api/flowability_batch", { text: text });
      const kpis = data.kpis || {};
      $("batch-ans-runs").textContent = kpis.total_runs != null ? String(kpis.total_runs) : "—";
      $("batch-ans-formulations").textContent = kpis.formulations != null ? String(kpis.formulations) : "—";

      renderDataTable(
        "batch-summary-head",
        "batch-summary-body",
        data.summary_header || [],
        data.summary_rows || [],
        true
      );
      renderDataTable(
        "batch-runs-head",
        "batch-runs-body",
        data.runs_header || [],
        data.runs_rows || [],
        true
      );
      $("results-batch").classList.remove("hidden");
    } catch (err) {
      showError("batch", err.message || "Batch CSV analysis failed.");
    } finally {
      btn.disabled = false;
    }
  }

  /* -----------------------------------------------------------------------
   * QUIT SERVER
   * --------------------------------------------------------------------- */
  async function quitServer() {
    const btn = $("btn-quit");
    btn.disabled = true;
    try {
      await fetch("/api/quit", { method: "POST" });
    } catch (_err) {
      // Server may close the socket immediately while shutting down.
    }
    $("quit-overlay").classList.remove("hidden");
  }

  /* -----------------------------------------------------------------------
   * EXPERIMENT MEDIA LOADER (EXTERNAL media/ SUBDIRECTORY)
   * --------------------------------------------------------------------- */
  function renderMediaSlot(slotId, items, emptyMessage) {
    const slot = $(slotId);
    if (!slot) return;
    slot.innerHTML = "";

    if (!Array.isArray(items) || items.length === 0) {
      const p = document.createElement("p");
      p.className = "media-empty";
      p.textContent = emptyMessage;
      slot.appendChild(p);
      return;
    }

    const grid = document.createElement("div");
    grid.className = "media-grid";

    items.forEach(function (item) {
      const fig = document.createElement("figure");
      fig.className = "media-figure";

      if (item.kind === "video") {
        const video = document.createElement("video");
        video.controls = true;
        video.preload = "metadata";
        video.src = item.url;
        video.addEventListener("error", function () {
          video.remove();
          const msg = document.createElement("p");
          msg.className = "media-empty";
          msg.textContent = "Unable to play video file: media/" + item.filename;
          fig.insertBefore(msg, fig.firstChild);
        });
        fig.appendChild(video);
      } else if (item.kind === "pdf") {
        const link = document.createElement("a");
        link.href = item.url;
        link.target = "_blank";
        link.rel = "noopener noreferrer";
        link.className = "btn btn-secondary";
        link.textContent = "Open PDF Document: " + item.label;
        fig.appendChild(link);
      } else {
        const img = document.createElement("img");
        img.src = item.url;
        img.alt = item.label || item.filename;
        img.loading = "lazy";
        img.addEventListener("error", function () {
          img.remove();
          const msg = document.createElement("p");
          msg.className = "media-empty";
          msg.textContent = "Unable to load image file: media/" + item.filename;
          fig.insertBefore(msg, fig.firstChild);
        });
        fig.appendChild(img);
      }

      grid.appendChild(fig);
    });

    slot.appendChild(grid);
  }

  async function loadExperimentMedia() {
    const fallbackNote = "No media files found in media/ for this experiment. Add images (.svg, .png, .jpg, .gif, .webp) or videos (.mp4, .webm) to the media/ folder next to Group_5_tool.py.";
    try {
      const response = await fetch("/api/media");
      if (!response.ok) {
        throw new Error("HTTP " + response.status);
      }
      const data = await response.json();
      if (!data || data.ok === false) {
        throw new Error((data && data.error) ? data.error : "Unable to inspect media/ folder.");
      }
      const exps = data.experiments || {};
      const missingDirNote = data.media_dir_exists
        ? fallbackNote
        : "The optional media/ folder was not found next to Group_5_tool.py. Place experiment diagrams or photos in a media/ subfolder to display them here.";

      renderMediaSlot("media-slot-sieve", exps.sieve || [], missingDirNote);
      renderMediaSlot("media-slot-carr", exps.carr || [], missingDirNote);
      renderMediaSlot("media-slot-angle", exps.angle || [], missingDirNote);
      renderMediaSlot("media-slot-flow", exps.flow || [], missingDirNote);
      renderMediaSlot("media-slot-batch", data.all_files || [], missingDirNote);
    } catch (_err) {
      ["media-slot-sieve", "media-slot-carr", "media-slot-angle", "media-slot-flow", "media-slot-batch"].forEach(function (id) {
        renderMediaSlot(id, [], "External media preview is currently unavailable.");
      });
    }
  }

  /* -----------------------------------------------------------------------
   * BOOTSTRAP
   * --------------------------------------------------------------------- */
  document.addEventListener("DOMContentLoaded", function () {
    initTabs();

    // Tab 1: Sieve Analysis
    $("sieve-csv-input").addEventListener("change", handleSieveCsvUpload);
    $("btn-load-pd01").addEventListener("click", function () { loadSieveSample("sieve_pd01"); });
    $("btn-load-ref").addEventListener("click", function () { loadSieveSample("sieve_reference"); });
    $("btn-add-sieve-row").addEventListener("click", addSieveRow);
    $("btn-calc-sieve").addEventListener("click", calculateSieve);

    // Tab 2: Carr & Hausner
    $("btn-calc-carr").addEventListener("click", calculateCarr);
    $("btn-example-carr").addEventListener("click", loadCarrExample);

    // Tab 3: Angle of Repose
    $("angle-unit").addEventListener("change", updateAngleUnitLabels);
    $("btn-calc-angle").addEventListener("click", calculateAngle);
    $("btn-example-angle").addEventListener("click", loadAngleExample);

    // Tab 4: Orifice Flow Rate
    $("flow-method").addEventListener("change", syncFlowMethodInputs);
    $("btn-calc-flow").addEventListener("click", calculateFlowRate);
    $("btn-example-flow").addEventListener("click", loadFlowExample);
    syncFlowMethodInputs();

    // Tab 5: Batch Flowability CSV
    $("batch-csv-input").addEventListener("change", handleBatchCsvUpload);
    $("btn-load-course-csv").addEventListener("click", loadCourseBatchDataset);
    $("btn-analyze-batch").addEventListener("click", analyzeBatchCsv);

    // Quit button
    $("btn-quit").addEventListener("click", quitServer);

    // Load external experiment diagrams/photos from media/ subdirectory
    loadExperimentMedia();

    // Auto-load PD01 sieve sample for immediate one-click evaluation
    loadSieveSample("sieve_pd01");
  });
})();
"""

EMBEDDED_WEB_ASSETS: dict[str, tuple[bytes, str]] = {
    "/": (HTML_INDEX.encode("utf-8"), "text/html; charset=utf-8"),
    "/style.css": (CSS_STYLES.encode("utf-8"), "text/css; charset=utf-8"),
    "/app.js": (JS_APP.encode("utf-8"), "application/javascript; charset=utf-8"),
}


# =============================================================================
# RESPONSE FORMATTING & INPUT EXTRACTION HELPERS
# =============================================================================

def _require_key(inputs: dict, keys: tuple[str, ...], field_label: str):
    """
    Extracts the first present key from `keys` in `inputs`, raising a clear
    ValueError if none of the keys exist.
    """
    for key in keys:
        if key in inputs:
            return inputs[key]
    raise ValueError(f"Missing required input field: {field_label}.")


def _is_missing_or_nan(val) -> bool:
    """
    Returns True if a value is None, NaN, or +/-Inf.
    """
    if val is None:
        return True
    if isinstance(val, (float, np.floating)):
        return not math.isfinite(float(val))
    if isinstance(val, str):
        return val.strip().lower() in ("nan", "inf", "-inf", "none")
    try:
        fval = float(val)
        return not math.isfinite(fval)
    except (TypeError, ValueError):
        return False


def _format_num(val, decimals: int = 2, fallback: str = "N/A") -> str:
    """
    Formats a finite numeric value with fixed decimals, or returns `fallback`.
    """
    if _is_missing_or_nan(val):
        return fallback
    fval = float(val)
    if fval == 0.0:
        fval = 0.0
    return f"{fval:.{decimals}f}"


def _format_input_float(val: float) -> str:
    """
    Formats a validated float for editable table inputs with 12 significant
    digits so CSV round-trips do not lose measurement precision.
    """
    fval = float(val)
    if fval == 0.0:
        fval = 0.0
    return f"{fval:.12g}"


def _format_cell(val, col_name: str = "") -> str:
    """
    Converts any table cell value to a safe JSON string with no NaN/Inf values.
    Replaces missing SD entries with 'N/A (n=1)'.
    """
    col_lower = str(col_name).lower()
    if _is_missing_or_nan(val):
        if "sd" in col_lower or "std" in col_lower or "±" in col_lower:
            return "N/A (n=1)"
        return "N/A"
    if isinstance(val, (float, np.floating)):
        fval = float(val)
        decimals = 6 if "q3 [%/µm]" in col_lower else 4
        formatted = f"{fval:.{decimals}f}".rstrip("0").rstrip(".")
        return "0" if formatted in ("", "-0") else formatted
    return str(val)


def _format_mean_sd_pair(mean_val, sd_val, decimals: int = 2) -> str:
    """
    Formats paired Mean and SD values as 'Mean ± SD' or 'Mean ± N/A (n=1)' when
    SD is null/NaN (n = 1 sample).
    """
    if _is_missing_or_nan(mean_val):
        return "N/A"
    mean_str = _format_num(mean_val, decimals=decimals)
    if _is_missing_or_nan(sd_val):
        return f"{mean_str} ± N/A (n=1)"
    return f"{mean_str} ± {_format_num(sd_val, decimals=decimals)}"


def _table_from_dataframe(df: pd.DataFrame) -> tuple[list[str], list[list[str]]]:
    """
    Converts a pandas DataFrame into JSON-safe (headers, rows) string arrays.
    """
    if df is None or df.empty:
        return [], []
    headers = [str(col) for col in df.columns]
    records = df.to_dict(orient="records")
    rows = [[_format_cell(rec.get(h), h) for h in headers] for rec in records]
    return headers, rows


# =============================================================================
# API ROUTE ADAPTERS
# =============================================================================

def run_sieve(inputs: dict) -> dict:
    """
    API adapter for `/api/sieve`.
    Validates sieve rows, x_max_um, and interpolation mode, runs sieve analysis,
    renders Q3/q3 PNG plots as base64, and returns KPIs + summary table.
    """
    if "rows" not in inputs or not isinstance(inputs["rows"], list):
        raise ValueError("Missing or invalid 'rows': expected a list of [mesh_um, empty_g, full_g] rows.")
    if len(inputs["rows"]) < MIN_SIEVE_STAGES:
        raise ValueError(f"At least {MIN_SIEVE_STAGES} sieve rows (including pan at 0 µm) are required.")

    mesh_sizes: list[float] = []
    empty_weights: list[float] = []
    full_weights: list[float] = []

    for idx, row in enumerate(inputs["rows"]):
        if not isinstance(row, (list, tuple)) or len(row) != 3:
            raise ValueError(f"Row {idx + 1}: expected exactly 3 values [mesh_um, empty_g, full_g].")
        mesh_sizes.append(parse_finite_float(row[0], f"Row {idx + 1}, Mesh Size [µm]"))
        empty_weights.append(parse_finite_float(row[1], f"Row {idx + 1}, Empty Sieve [g]"))
        full_weights.append(parse_finite_float(row[2], f"Row {idx + 1}, Sieve + Residue [g]"))

    raw_xmax = inputs.get("x_max_um", inputs.get("x_max", DEFAULT_X_MAX_UM))
    if raw_xmax is None or (isinstance(raw_xmax, str) and not raw_xmax.strip()):
        x_max_um = DEFAULT_X_MAX_UM
    else:
        x_max_um = parse_finite_float(raw_xmax, "Upper Size Limit x_max [µm]")
    if x_max_um <= 0.0:
        raise ValueError(f"Upper Size Limit x_max must be > 0 µm (received {x_max_um}).")

    raw_mode = inputs.get("interp_mode", inputs.get("interpolation_mode", inputs.get("mode", "linear")))
    if not isinstance(raw_mode, str):
        raise ValueError("Interpolation mode must be a string ('linear' or 'log').")
    interp_mode = raw_mode.strip().lower()
    if interp_mode not in ("linear", "log"):
        raise ValueError(f"Invalid interpolation mode '{raw_mode}'. Expected 'linear' or 'log'.")

    df_sieve, metrics = analyze_sieve_data(
        mesh_sizes,
        empty_weights,
        full_weights,
        x_max_um=x_max_um,
        interpolation_mode=interp_mode,
    )

    fig_q3 = None
    fig_density = None
    try:
        fig_q3 = plot_q3_cumulative(df_sieve, metrics)
        chart_q3 = figure_to_base64_png(fig_q3)
        fig_q3 = None
        fig_density = plot_q3_density(df_sieve, metrics)
        chart_density = figure_to_base64_png(fig_density)
        fig_density = None
    finally:
        if fig_q3 is not None:
            plt.close(fig_q3)
        if fig_density is not None:
            plt.close(fig_density)

    table_header, table_rows = _table_from_dataframe(df_sieve)
    interp_note = metrics.get("interpolation_note")

    return {
        "ok": True,
        "answers": {
            "d10": _format_num(metrics["D10_um"], 1),
            "d50": _format_num(metrics["D50_um"], 1),
            "d90": _format_num(metrics["D90_um"], 1),
            "span": _format_num(metrics["span_um"], 1),
            "relative_span": _format_num(metrics["relative_span"], 3),
            "mass": _format_num(metrics["total_mass_g"], 2),
        },
        "interpolation_mode": metrics["interpolation_mode"],
        "interpolation_note": "" if interp_note is None else str(interp_note),
        "charts": [chart_q3, chart_density],
        "table_title": "Sieve Analysis Summary Table (Ph. Eur. 2.9.38 / ISO 9276-1)",
        "table_header": table_header,
        "table_rows": table_rows,
    }


def run_read_sieve_csv(inputs: dict) -> dict:
    """
    API adapter for `/api/read_sieve_csv`.
    Parses uploaded sieve CSV text via `load_sieve_csv` and returns table rows.
    """
    text = _require_key(inputs, ("text", "csv_text"), "text")
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Uploaded CSV text must be a non-empty string.")

    mesh, empty, full = load_sieve_csv(text)
    rows = [[_format_input_float(m), _format_input_float(e), _format_input_float(f)] for m, e, f in zip(mesh, empty, full)]
    return {"ok": True, "rows": rows}


def run_load_sample(inputs: dict) -> dict:
    """
    API adapter for `/api/sample`.
    Returns embedded offline sample datasets for Sieve Analysis or Batch Flowability.
    """
    raw_sample = _require_key(inputs, ("sample", "name"), "sample")
    if not isinstance(raw_sample, str) or not raw_sample.strip():
        raise ValueError("Sample identifier must be a non-empty string.")

    key = raw_sample.strip().lower()
    if key in ("sieve_pd01", "pd01"):
        mesh, empty, full = load_sieve_csv(SAMPLE_SIEVE_PD01_CSV)
        rows = [[_format_input_float(m), _format_input_float(e), _format_input_float(f)] for m, e, f in zip(mesh, empty, full)]
        return {
            "ok": True,
            "sample": "sieve_pd01",
            "rows": rows,
            "csv_text": SAMPLE_SIEVE_PD01_CSV,
        }

    if key in ("sieve_reference", "reference"):
        mesh, empty, full = load_sieve_csv(SAMPLE_SIEVE_REFERENCE_CSV)
        rows = [[_format_input_float(m), _format_input_float(e), _format_input_float(f)] for m, e, f in zip(mesh, empty, full)]
        return {
            "ok": True,
            "sample": "sieve_reference",
            "rows": rows,
            "csv_text": SAMPLE_SIEVE_REFERENCE_CSV,
        }

    if key in ("flowability_batch", "flowability", "course_dataset"):
        return {
            "ok": True,
            "sample": "flowability_batch",
            "csv_text": SAMPLE_FLOWABILITY_CSV,
        }

    raise ValueError(
        f"Unknown sample '{raw_sample}'. Expected 'sieve_pd01', 'sieve_reference', or 'flowability_batch'."
    )


def run_carr_hausner(inputs: dict) -> dict:
    """
    API adapter for `/api/carr_hausner`.
    Validates powder mass, poured volume V0, and tapped volume Vf, then calls
    `analyze_carr_hausner` (Ph. Eur. 2.9.34 & 2.9.36 Table 2.9.36.-2).
    """
    mass_g = parse_finite_float(_require_key(inputs, ("mass", "mass_g"), "mass"), "Powder mass m [g]")
    v0_ml = parse_finite_float(_require_key(inputs, ("v0", "v0_ml"), "v0"), "Poured volume V0 [mL]")
    vf_ml = parse_finite_float(_require_key(inputs, ("vf", "vf_ml", "vt_ml"), "vf"), "Tapped volume Vf [mL]")

    results = analyze_carr_hausner(mass_g, v0_ml, vf_ml)

    return {
        "ok": True,
        "answers": {
            "bulk": _format_num(results["bulk_density_g_ml"], 4),
            "tapped": _format_num(results["tapped_density_g_ml"], 4),
            "carr": _format_num(results["carr_index_pct"], 2),
            "hausner": _format_num(results["hausner_ratio"], 3),
            "ci_rating": str(results["ci_rating"]),
            "hr_rating": str(results["hr_rating"]),
            "property": str(results["ci_rating"]),
        },
    }


def run_angle(inputs: dict) -> dict:
    """
    API adapter for `/api/angle`.
    Validates unit ('cm' or 'mm'), h_total, h_base (default 0), and diameter,
    then calls `analyze_angle_of_repose` (Ph. Eur. 2.9.36 Table 2.9.36.-1).
    """
    raw_unit = inputs.get("unit", "cm")
    if not isinstance(raw_unit, str):
        raise ValueError("Unit must be a string ('cm' or 'mm').")
    unit = raw_unit.strip().lower()
    if unit not in SUPPORTED_LENGTH_UNITS:
        raise ValueError(f"Unsupported unit '{raw_unit}'. Expected 'cm' or 'mm'.")

    h_total = parse_finite_float(
        _require_key(inputs, ("h_total", "total_height", "height"), "h_total"),
        f"Total height h_total [{unit}]",
    )
    raw_hbase = inputs.get("h_base", inputs.get("base_height", 0.0))
    if raw_hbase is None or (isinstance(raw_hbase, str) and not raw_hbase.strip()):
        h_base = 0.0
    else:
        h_base = parse_finite_float(raw_hbase, f"Pedestal base height h_base [{unit}]")

    diameter = parse_finite_float(
        _require_key(inputs, ("diameter", "d"), "diameter"),
        f"Base diameter d [{unit}]",
    )

    results = analyze_angle_of_repose(
        total_height=h_total,
        diameter=diameter,
        base_height=h_base,
        unit=unit,
    )

    return {
        "ok": True,
        "answers": {
            "h_cone": _format_num(results["cone_height"], 2),
            "unit": results["unit"],
            "angle": _format_num(results["angle_of_repose_deg"], 2),
            "property": str(results["aor_rating"]),
        },
    }


def run_flow_rate(inputs: dict) -> dict:
    """
    API adapter for `/api/flow_rate`.
    Supports 'fixed_mass' (100 g), 'fixed_time' (10 s), and 'custom' (m & t)
    per Ph. Eur. 2.9.16 and Ph. Eur. 2.9.36.
    """
    raw_method = _require_key(inputs, ("method",), "method")
    if not isinstance(raw_method, str):
        raise ValueError("Method must be a string ('fixed_mass', 'fixed_time', or 'custom').")
    method = raw_method.strip().lower()
    if method not in ("fixed_mass", "fixed_time", "custom"):
        raise ValueError(
            f"Unsupported flow rate method '{raw_method}'. Expected 'fixed_mass', 'fixed_time', or 'custom'."
        )

    if method == "fixed_mass":
        raw_t = _require_key(inputs, ("time", "flow_time_s", "measured_value", "value"), "time")
        flow_time_s = parse_finite_float(raw_t, "Flow time t [s]")
        result = analyze_flow_rate("fixed_mass", measured_value=flow_time_s)
    elif method == "fixed_time":
        raw_m = _require_key(inputs, ("mass", "mass_g", "measured_value", "value"), "mass")
        mass_g = parse_finite_float(raw_m, "Powder mass m [g]")
        result = analyze_flow_rate("fixed_time", measured_value=mass_g)
    else:
        mass_g = parse_finite_float(_require_key(inputs, ("mass", "mass_g"), "mass"), "Powder mass m [g]")
        flow_time_s = parse_finite_float(
            _require_key(inputs, ("time", "flow_time_s"), "time"),
            "Flow time t [s]",
        )
        result = analyze_flow_rate("custom", mass_g=mass_g, flow_time_s=flow_time_s)

    return {
        "ok": True,
        "answers": {
            "rate": _format_num(result["flow_rate_g_s"], 2),
            "time_per_100g": _format_num(result["time_per_100g_s"], 2),
            "mass_used": _format_num(result["mass_g"], 2),
            "time_used": _format_num(result["flow_time_s"], 2),
        },
    }


def run_flowability_batch(inputs: dict) -> dict:
    """
    API adapter for `/api/flowability_batch`.
    Validates batch CSV string, invokes `analyze_flowability_batch`, and returns
    KPIs, grouped Mean ± SD summary table (with 'N/A (n=1)' when SD is null),
    and individual run results table.
    """
    text = _require_key(inputs, ("text", "csv_text"), "text")
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Batch CSV text must be a non-empty string.")

    raw_result = analyze_flowability_batch(text)
    summary_records = raw_result["summary"]
    run_records = raw_result["rows"]

    summary_header = [
        "Sample Code",
        "n",
        "Bulk Density [g/mL] (Mean ± SD)",
        "Tapped Density [g/mL] (Mean ± SD)",
        "Carr's Index [%] (Mean ± SD)",
        "Hausner Ratio [-] (Mean ± SD)",
        "CI Rating",
        "HR Rating",
        "Angle of Repose [°] (Mean ± SD)",
        "AoR Rating",
    ]
    summary_rows = [
        [
            str(rec["sample_code"]),
            str(rec["n"]),
            _format_mean_sd_pair(rec["bulk_density_mean_g_ml"], rec["bulk_density_sd_g_ml"], decimals=4),
            _format_mean_sd_pair(rec["tapped_density_mean_g_ml"], rec["tapped_density_sd_g_ml"], decimals=4),
            _format_mean_sd_pair(rec["carr_index_mean_pct"], rec["carr_index_sd_pct"], decimals=2),
            _format_mean_sd_pair(rec["hausner_ratio_mean"], rec["hausner_ratio_sd"], decimals=3),
            str(rec["ci_rating"]),
            str(rec["hr_rating"]),
            _format_mean_sd_pair(rec["angle_of_repose_mean_deg"], rec["angle_of_repose_sd_deg"], decimals=2),
            str(rec["aor_rating"]),
        ]
        for rec in summary_records
    ]

    runs_header = [
        "Sample Code",
        "Group / Run",
        "Mass [g]",
        "V0 [mL]",
        "Vt [mL]",
        "h_cone [mm]",
        "d [mm]",
        "Bulk Density [g/mL]",
        "Tapped Density [g/mL]",
        "Carr's Index [%]",
        "Hausner Ratio [-]",
        "CI Rating",
        "HR Rating",
        "Angle of Repose [°]",
        "AoR Rating",
    ]
    runs_rows = [
        [
            str(rec["sample_code"]),
            str(rec["group_run"]),
            _format_num(rec["mass_g"], 2),
            _format_num(rec["v0_ml"], 1),
            _format_num(rec["vt_ml"], 1),
            _format_num(rec["cone_height_mm"], 1),
            _format_num(rec["diameter_mm"], 1),
            _format_num(rec["bulk_density_g_ml"], 4),
            _format_num(rec["tapped_density_g_ml"], 4),
            _format_num(rec["carr_index_pct"], 2),
            _format_num(rec["hausner_ratio"], 3),
            str(rec["ci_rating"]),
            str(rec["hr_rating"]),
            _format_num(rec["angle_of_repose_deg"], 2),
            str(rec["aor_rating"]),
        ]
        for rec in run_records
    ]

    return {
        "ok": True,
        "kpis": {
            "total_runs": int(raw_result["row_count"]),
            "formulations": int(raw_result["group_count"]),
        },
        "summary_header": summary_header,
        "summary_rows": summary_rows,
        "runs_header": runs_header,
        "runs_rows": runs_rows,
    }


# Route dispatch table mapping all 7 calculation & data API endpoints
CALCULATIONS = {
    "/api/sieve": run_sieve,
    "/api/read_sieve_csv": run_read_sieve_csv,
    "/api/sample": run_load_sample,
    "/api/carr_hausner": run_carr_hausner,
    "/api/angle": run_angle,
    "/api/flow_rate": run_flow_rate,
    "/api/flowability_batch": run_flowability_batch,
}


# =============================================================================
# EXTERNAL MEDIA DIRECTORY RESOLVER (DYNAMIC RELATIVE PATHS ONLY)
# =============================================================================

def _media_kind_for_suffix(suffix: str) -> str:
    """Maps a validated lowercase file extension to a UI rendering category."""
    if suffix in (".mp4", ".webm"):
        return "video"
    if suffix == ".pdf":
        return "pdf"
    return "image"


def _humanize_media_label(rel_posix: str) -> str:
    """Creates a clean display label from a relative media path."""
    stem = Path(rel_posix).stem.replace("_", " ").replace("-", " ").strip()
    return " ".join(word.capitalize() for word in stem.split()) if stem else rel_posix


def resolve_media_file(rel_path: str) -> tuple[Path, str] | None:
    """
    Resolves a requested relative path strictly inside MEDIA_DIR
    (Path(__file__).resolve().parent / "media") without using hardcoded
    absolute paths. Blocks directory traversal, symlinks escaping media/,
    and disallowed file extensions.
    """
    if not isinstance(rel_path, str) or not rel_path:
        return None

    decoded = unquote(rel_path)
    if "\x00" in decoded or "\\" in decoded or decoded.startswith("/"):
        return None

    rel_candidate = Path(decoded)
    if rel_candidate.is_absolute() or ".." in rel_candidate.parts:
        return None

    try:
        media_root = MEDIA_DIR.resolve()
        if not media_root.is_dir():
            return None
        candidate = (media_root / rel_candidate).resolve()
        # Strict containment check prevents any path traversal or symlink escape
        candidate.relative_to(media_root)
    except (OSError, ValueError):
        return None

    if not candidate.is_file():
        return None

    suffix = candidate.suffix.lower()
    content_type = ALLOWED_MEDIA_TYPES.get(suffix)
    if content_type is None:
        return None

    return candidate, content_type


def list_available_media() -> dict:
    """
    Scans MEDIA_DIR dynamically relative to the running script and categorizes
    discovered offline media files by experiment tab. Never crashes if media/
    is absent or unreadable.
    """
    experiments: dict[str, list[dict]] = {key: [] for key in EXPERIMENT_MEDIA_PREFIXES}
    all_files: list[dict] = []

    try:
        media_root = MEDIA_DIR.resolve()
        if not media_root.is_dir():
            return {
                "ok": True,
                "media_dir_exists": False,
                "experiments": experiments,
                "all_files": all_files,
            }
    except OSError:
        return {
            "ok": True,
            "media_dir_exists": False,
            "experiments": experiments,
            "all_files": all_files,
        }

    try:
        discovered = sorted(media_root.rglob("*"), key=lambda p: p.as_posix().lower())
    except OSError:
        discovered = []

    for entry in discovered:
        try:
            if not entry.is_file():
                continue
            resolved = entry.resolve()
            rel_path = resolved.relative_to(media_root)
        except (OSError, ValueError):
            continue

        # Skip hidden files or dotfiles (e.g., .DS_Store)
        if any(part.startswith(".") for part in rel_path.parts):
            continue

        suffix = resolved.suffix.lower()
        if suffix not in ALLOWED_MEDIA_TYPES:
            continue

        rel_posix = rel_path.as_posix()
        item = {
            "filename": rel_posix,
            "label": _humanize_media_label(rel_posix),
            "kind": _media_kind_for_suffix(suffix),
            "url": "/media/" + quote(rel_posix, safe="/"),
        }
        all_files.append(item)

        stem_lower = resolved.stem.lower()
        parts_lower = [p.lower() for p in rel_path.parts[:-1]]
        for exp_key, prefixes in EXPERIMENT_MEDIA_PREFIXES.items():
            if any(
                stem_lower.startswith(pref) or pref in parts_lower
                for pref in prefixes
            ):
                experiments[exp_key].append(item)

    return {
        "ok": True,
        "media_dir_exists": True,
        "experiments": experiments,
        "all_files": all_files,
    }


# =============================================================================
# LOCAL OFFLINE HTTP SERVER
# =============================================================================

class RequestHandler(BaseHTTPRequestHandler):
    """
    HTTP request handler serving embedded offline UI assets, dynamic relative
    media/ experiment files, and JSON calculation endpoints over 127.0.0.1
    with strict anti-crash payload validation.
    """

    def send_answer(self, status: int, content_type: str, data: bytes) -> None:
        try:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.end_headers()
            self.wfile.write(data)
        except (BrokenPipeError, ConnectionResetError):
            return

    def send_json(self, status: int, payload: dict) -> None:
        try:
            body = json.dumps(payload, allow_nan=False).encode("utf-8")
        except (TypeError, ValueError) as exc:
            fallback = {"ok": False, "error": f"Response serialization failed: {exc}"}
            body = json.dumps(fallback, allow_nan=False).encode("utf-8")
            status = 500
        self.send_answer(status, "application/json; charset=utf-8", body)

    def do_GET(self) -> None:
        path = self.path.split("?", 1)[0]

        if path == "/api/media":
            self.send_json(200, list_available_media())
            return

        if path.startswith("/media/"):
            resolved = resolve_media_file(path[len("/media/"):])
            if resolved is None:
                self.send_answer(404, "text/plain; charset=utf-8", b"Media file not found")
                return
            file_path, content_type = resolved
            try:
                data = file_path.read_bytes()
            except OSError:
                self.send_answer(404, "text/plain; charset=utf-8", b"Unable to read media file")
                return
            self.send_answer(200, content_type, data)
            return

        asset = EMBEDDED_WEB_ASSETS.get(path)
        if asset is None:
            self.send_answer(404, "text/plain; charset=utf-8", b"Page not found")
            return
        data, content_type = asset
        self.send_answer(200, content_type, data)

    def do_POST(self) -> None:
        path = self.path.split("?", 1)[0]

        if path == "/api/quit":
            self.send_json(200, {"ok": True})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return

        if path not in CALCULATIONS:
            self.send_json(404, {"ok": False, "error": f"Unknown API endpoint: {path}"})
            return

        raw_length = self.headers.get("Content-Length", "0")
        try:
            length = int(raw_length)
        except ValueError:
            self.send_json(400, {"ok": False, "error": "Invalid Content-Length header."})
            return

        if length < 0:
            self.send_json(400, {"ok": False, "error": "Content-Length cannot be negative."})
            return
        if length > MAX_REQUEST_BYTES:
            self.send_json(
                413,
                {"ok": False, "error": f"Request payload exceeds maximum limit ({MAX_REQUEST_BYTES} bytes)."},
            )
            return

        raw_body = self.rfile.read(length)
        try:
            text_body = raw_body.decode("utf-8")
        except UnicodeDecodeError:
            self.send_json(400, {"ok": False, "error": "Request body must be valid UTF-8 text."})
            return

        try:
            inputs = json.loads(text_body) if text_body.strip() else {}
        except json.JSONDecodeError as exc:
            self.send_json(400, {"ok": False, "error": f"Invalid JSON request body: {exc.msg}"})
            return

        if not isinstance(inputs, dict):
            self.send_json(400, {"ok": False, "error": "JSON request body must be a JSON object."})
            return

        try:
            answer = CALCULATIONS[path](inputs)
        except (ValueError, KeyError, TypeError, OverflowError) as exc:
            self.send_json(200, {"ok": False, "error": str(exc)})
            return
        except Exception as exc:
            self.send_json(500, {"ok": False, "error": f"Unexpected server error: {exc}"})
            return

        self.send_json(200, answer)

    def log_message(self, format: str, *args) -> None:
        """Suppresses default per-request HTTP stderr logging."""
        return


# =============================================================================
# ENTRYPOINT
# =============================================================================

def main() -> None:
    """
    Verifies dependencies, starts the local offline HTTP server on 127.0.0.1,
    opens the dashboard in the default browser, and shuts down cleanly on exit.
    """
    check_dependencies()

    try:
        server = HTTPServer(("127.0.0.1", PORT), RequestHandler)
    except OSError:
        # Fallback to OS-assigned ephemeral port if PORT is already occupied
        server = HTTPServer(("127.0.0.1", 0), RequestHandler)

    address = f"http://127.0.0.1:{server.server_address[1]}/"
    print("Group 5 — Powder Characterization Tool (Offline Local Mode)")
    print(f"Open in your browser: {address}")
    print("To stop: click 'Quit Server' on the web page or press Ctrl+C here.")

    try:
        webbrowser.open(address)
    except Exception:
        pass

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()

    print("Group 5 Powder Characterization Tool stopped cleanly.")


if __name__ == "__main__":
    main()
