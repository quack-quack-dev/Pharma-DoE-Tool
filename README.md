# 💊 Group 5: Powder Characterization Tool

This repository contains our automated evaluation tool for powder characterization. The application—including the mathematical models, user interface, and sample datasets—runs completely offline from `Group_5_tool.py`, with explanatory experiment diagrams and photos loaded dynamically from the `media/` subfolder so the whole project can be submitted as a portable ZIP folder.

## 🚀 Features & How It Works

The tool is designed to process lab data without needing an internet connection. It is broken down into a few core analytical modules:

- **Sieve Analysis (Ph. Eur. 2.9.38):** Calculates particle size distribution. It computes D10, D50, D90, and the Span using linear and logarithmic interpolation, and plots the Q3 cumulative curve and q3 density histogram.
- **Density & Compressibility (Ph. Eur. 2.9.34 & 2.9.36):** Automatically calculates Bulk/Tapped Density, Carr's Index, and Hausner Ratio, classifying the flow character according to the pharmacopeia tables.
- **Angle of Repose (Ph. Eur. 2.9.36):** Calculates the angle factoring in the pedestal base height and grades the flowability.
- **Orifice Flow Rate (Ph. Eur. 2.9.16):** Processes Fixed Mass, Fixed Time, and Custom modes to standardize the flow rate to grams per second.
- **Batch CSV Analyzer:** Upload a dataset (like our embedded `flowability_complete_dataset.csv`) to compute per-run metrics and group sample statistics (Mean ± SD) automatically.
- **Experiment Media Guides (`media/`):** Each experiment tab automatically displays schematics, lab photos, or videos (`.svg`, `.png`, `.jpg`, `.gif`, `.webp`, `.mp4`, `.webm`, `.pdf`) placed in the `media/` folder next to `Group_5_tool.py`.

All mathematical operations, constants, and unit conversions include inline academic citations for easy reference during the code walkthrough.

## 💻 How to Run the Tool

Because all paths are resolved dynamically relative to the script, getting it running after extracting the ZIP on a fresh laptop is very straightforward.

1. **Install the required libraries:**
   Make sure you have Python installed, then open your terminal and run:

   ```bash
   pip install numpy pandas matplotlib
   ```

2. **Run the script:**
   Extract the ZIP folder, navigate to the folder containing `Group_5_tool.py` and `media/`, and run:

   ```bash
   python Group_5_tool.py
   ```

   *(Depending on your system setup, you might need to use `python3` instead of `python`.)*

## 📝 Notes for the Group

- **Keep `media/` next to `Group_5_tool.py`:** All paths to `media/` are resolved dynamically with `pathlib.Path(__file__).resolve().parent / "media"` (no hardcoded absolute paths), so it works regardless of where the ZIP is extracted or which directory the script is launched from.
- **Adding or replacing media:** Drop images or videos into `media/` using prefixes like `sieve_analysis`, `carr_hausner`, `angle_of_repose`, or `orifice_flow` to show them automatically under the matching experiment tab (and all files in `media/` appear in the Batch tab's Lab Media Gallery).
- **Error Handling:** The tool is equipped with strict input validation and safe media fallbacks. If bad data is entered—or even if `media/` is empty or missing—it will return a clean message rather than crashing the program.
