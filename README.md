# 💊 Group 5: Powder Characterization Tool

This repository contains our automated evaluation tool for powder characterization. To comply with the strict project requirements, the entire application—including the mathematical models, the user interface, and the sample datasets—is bundled into a single, offline-capable Python script.

## 🚀 Features & How It Works

The tool is designed to process lab data without needing an internet connection. It is broken down into a few core analytical modules:

- **Sieve Analysis (Ph. Eur. 2.9.38):** Calculates particle size distribution. It computes D10, D50, D90, and the Span using linear and logarithmic interpolation, and plots the Q3 cumulative curve and q3 density histogram.
- **Density & Compressibility (Ph. Eur. 2.9.34 & 2.9.36):** Automatically calculates Bulk/Tapped Density, Carr's Index, and Hausner Ratio, classifying the flow character according to the pharmacopeia tables.
- **Angle of Repose (Ph. Eur. 2.9.36):** Calculates the angle factoring in the pedestal base height and grades the flowability.
- **Orifice Flow Rate (Ph. Eur. 2.9.16):** Processes Fixed Mass, Fixed Time, and Custom modes to standardize the flow rate to grams per second.
- **Batch CSV Analyzer:** Upload a dataset (like our embedded `flowability_complete_dataset.csv`) to compute per-run metrics and group sample statistics (Mean ± SD) automatically.

All mathematical operations, constants, and unit conversions include inline academic citations for easy reference during the code walkthrough.

## 💻 How to Run the Tool

Because everything is packaged into one file, getting it running on a fresh laptop is very straightforward.

1. **Install the required libraries:**
   Make sure you have Python installed, then open your terminal and run:

   ```bash
   pip install numpy pandas matplotlib
   ```

2. **Run the script:**
   Navigate to the folder containing the script and run:

   ```bash
   python Group_5_tool.py
   ```

   *(Depending on your system setup, you might need to use `python3` instead of `python`.)*

## 📝 Notes for the Group

- **Do not split the file:** The professor explicitly requires a single `Group_5_tool.py` script. All HTML, CSS, JavaScript, and data loaders are intentionally embedded inside this Python file to guarantee it runs flawlessly offline on the testing machine.
- **Error Handling:** The tool is equipped with strict input validation. If bad data is entered, it will return a clean error message rather than crashing the program.
