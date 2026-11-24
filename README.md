# Group 5: Powder Characterization Tool

This repository contains our automated evaluation tool for powder characterization. The application—including the mathematical models, user interface, and sample datasets—runs completely offline from `Group_5_tool.py`, with explanatory experiment diagrams and photos loaded dynamically from the `media/` subfolder so the whole project can be submitted as a portable ZIP folder.

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
