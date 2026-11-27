# HANDOFF.md — Group 5 Powder Characterization Tool

## Current Repository & Branch State (`quack-quack-dev/Pharma-DoE-Tool`)
*(Note: GitHub account was renamed from `quack-quack-pimper` to `quack-quack-dev`; ensure `quack-quack-dev/Pharma-DoE-Tool` is attached in Antigravity Project Resources for push access.)*

- **`main` (`599fd50` / `6ba862f`):** Stable baseline with 5 analytical tabs (`Sieve Analysis`, `Carr & Hausner`, `Angle of Repose`, `Orifice Flow Rate`, `Batch Flowability CSV`), dynamic relative `media/` loading (`GET /api/media`, `GET /media/*`), 4 offline SVG schematics, cleaned-up media headers (no redundant `media/` labels or figcaptions), fixed CSS chevron escape (`\\25B8`), and concise `README.md`. Passed 311/311 blind QA checks.
- **`LEO-prototype-merge`:** Branched from `main` (`6ba862f`). Cleanly merges Leonardo's prototype features directly into `Group_5_tool.py` (144 `# Source:` citations preserved) and adds full standalone Windows (`.exe`) + macOS PyInstaller support:
  - Tab 6: **Sieve Comparison** (`POST /api/compare`, 1–10 CSV files, combined vs. individual views, demo loader)
  - **Switchable histogram type** (`grouped` `dQ3 [%]` bars vs. `panels` `q3 [%/µm]`, auto-fallback when sieve stacks differ)
  - **Metadata footers** on exported PNGs (monospace footer + PNG `tEXt` tags + `metadata.txt` with SHA-256 data hash)
  - **Carr / Hausner / Angle comparisons** from optional CSV columns (`powder_mass_g`, `v0_ml`, `vf_ml`, `aor_height_mm`, `aor_diameter_mm`, `aor_base_height_mm`) + Ph. Eur. 2.9.36 Good-flow shaded bar chart
  - **ZIP download buttons** (`POST /api/sieve_report` and `POST /api/compare`, including German-Excel `comparison_table.csv` with UTF-8 BOM, `;` delimiter, decimal comma)
  - **Retained vs. lost powder balance & limit checks** (progress bar, Ph. Eur. 2.9.38 `<= 5%` loss check, optional min/max acceptance limits)
  - **Auto-calculating CSV upload** (auto-fills optional `initial_mass_g` & formulation details in Tab 1 and auto-calculates)
  - **PyInstaller Standalone Executable Support (`Group_5_tool.py`):**
    - Dynamic frozen/unfrozen `media/` resolution (`_get_executable_base_dir`, `_get_bundled_media_dir`, `get_media_roots`): checks optional external `media/` beside `.exe` / macOS binary / `.app` bundle first, falling back to bundled `sys._MEIPASS / "media"` so single-file executables work standalone with zero external folders.
    - `LocalHTTPServer(ThreadingHTTPServer)` with `_CALC_LOCK` around Matplotlib rendering so idle browser preconnect sockets never stall requests or shutdown.
    - `_safe_print` guard for detached `sys.stdout` (`pythonw`) and legacy Windows console codepages (`cp437`/`cp850`), `--no-browser` CLI flag for automated smoke tests, and `multiprocessing.freeze_support()`.
  - **Cross-Platform Build Workflow (`.github/workflows/build-executables.yml`):**
    - Builds `Group_5_tool.exe` (`windows-latest`) and `Group_5_tool` (`macos-latest`, preserved with `+x` permissions via `ditto -c -k`), runs isolated offline HTTP + Matplotlib + `_MEIPASS` smoke tests against each compiled binary, and uploads downloadable Actions artifacts (`Group_5_tool-Windows-exe` and `Group_5_tool-macOS`).

## Key Files
- `Group_5_tool.py` — Main offline application (embedded UI + math + server + PyInstaller-aware media resolver)
- `media/` — `sieve_analysis.svg`, `carr_hausner.svg`, `angle_of_repose.svg`, `orifice_flow.svg`
- `sample_data/` — `example_PD01.csv`, `sieve_PD01.csv`, `sieve_reference.csv`, `flowability_complete_dataset.csv`, `dummy_batch_{A,B,C}.csv`
- `AGENTS.md` — Permanent agent rules (Master Law)
- `HANDOFF.md` — Session state & executable build documentation
- `.github/workflows/build-executables.yml` — Windows & macOS PyInstaller CI build + smoke-test workflow

## Next Steps
1. Re-attach the renamed repository `quack-quack-dev/Pharma-DoE-Tool` in Antigravity Project Resources (replacing stale `quack-quack-pimper/Pharma-DoE-Tool`), then push `LEO-prototype-merge` (if GitHub App blocks `.github/workflows/build-executables.yml`, add that workflow file via GitHub UI or local git push).
2. Confirm both `windows-latest` (`Group_5_tool-Windows-exe`) and `macos-latest` (`Group_5_tool-macOS`) workflow artifacts pass and download cleanly from GitHub Actions.
