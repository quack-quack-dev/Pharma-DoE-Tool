# AGENTS.md — Group 5 Powder Characterization Tool (HHU Düsseldorf)

Permanent rules governing every AI agent working on `quack-quack-pimper/Pharma-DoE-Tool`.

## Core Project Rules (Master Law)

1. **Strict Offline Execution:** The tool must run 100% offline on a fresh grading laptop (Windows or macOS). Never add external API calls, internet downloads at runtime, external CDNs, or external web fonts.
2. **Dynamic Relative Paths Only (Never Hardcode Absolute Paths):** Always resolve bundled/external assets (`media/`) relative to the script or PyInstaller bundle (`sys._MEIPASS` / `Path(sys.executable).resolve().parent` / `Path(__file__).resolve().parent`). Never write absolute paths (`/home/...`, `/Users/...`, `C:\...`).
3. **Zero-Crash Protocol:** Never silently assume valid input or crash on bad data. Validate all numeric inputs (`parse_finite_float` rejecting `None`, `bool`, `bytes`, `NaN`, `±Inf`), CSV schemas (trailing commas, duplicate headers, blank cells), HTTP payloads (`Content-Length`, malformed UTF-8/JSON, `allow_nan=False`), and missing/empty `media/` folders with clear human-readable errors or fallbacks.
4. **Explicit Scientific Traceability:** Every mathematical equation, constant, unit conversion, and pharmacopeial threshold must include a brief inline `# Source:` comment (e.g., `# Source: Ph. Eur. 2.9.36, Table 2.9.36.-2`, `# Source: Ph. Eur. 2.9.38`, `# Source: ISO 9276-1`). Preserve all existing `# Source:` comments.
5. **Sanitizer → Blind QA Protocol:** Before finalizing major changes or merging to `main`, verify logic via a comment- and docstring-stripped copy (no comments/docstrings) using automated adversarial tests, keeping the fully commented `Group_5_tool.py` as the submission file.
6. **UI & Packaging Standards:**
   - Keep all HTML, CSS, and JS embedded inside `Group_5_tool.py` (no external `web/` folder); only explanatory experiment schematics/photos/videos live in `media/`.
   - Always call `matplotlib.use("Agg")` before importing `matplotlib.pyplot`.
   - Always run `py_compile` and `node --check` on `JS_APP` after editing `Group_5_tool.py` (watch Python triple-quoted string escapes like `\\n` and `\\25B8`).
7. **Git Branch Discipline:** Never overwrite `main` without explicit user approval. Develop and test new work on feature branches (currently `LEO-prototype-merge` for merged prototype features and standalone executables).
