# 🏛️ MASTER README & AGENT CONSTITUTION

> **Purpose of this document:** This file serves as the universal "Master Law" governing all AI agents (Architect, Executor, Integrator, Sanitizer, and QA) working on the Statistics and DoE semester project (Group 5 — Powder Characterization, HHU Düsseldorf). All prompts, workflows, and code generation MUST strictly obey these core principles.

---

## 🔒 CORE IMMUTABLE RULES

### 1️⃣ RULE #1: STRICT ENVIRONMENT CONSTRAINTS

* **The Principle:** The final deliverable must be exactly one Python script named `Group_5_tool.py`.
* **Execution:** It must run completely offline on a fresh testing laptop. Agents are strictly forbidden from writing code that requires external API calls, internet downloads during runtime, or multi-file dependencies.

### 2️⃣ RULE #2: THE ZERO-ASSUMPTION & ANTI-CRASH PROTOCOL

* **The Principle:** The tool must never crash.
* **Execution:** Agents must not silently fill in blanks or assume perfect user input. The code must include strict error handling: bad or invalid inputs must return a clear, human-readable error message and safely halt the specific operation without crashing the program.

### 3️⃣ RULE #3: EXPLICIT GROUNDING & TRACEABILITY

* **The Principle:** Every mathematical equation, constant, unit conversion, and threshold must be explicitly traced to its academic source to prepare for the live code walkthrough.
* **Execution:** When writing code, Executor agents MUST include brief inline comments strictly citing the source (e.g., `# Source: Ph. Eur. 2.9.36, Table 2.9.36.-2`). Agents must never invent variables or hallucinate formulas.

### 4️⃣ RULE #4: THE SANITIZER PROTOCOL (COMMENT MANAGEMENT)

* **The Principle:** QA Agents must evaluate pure logic without being influenced or "polluted" by the Executor's comments.
* **Execution:**
  1. **Executors** write the fully commented code (per Rule #3). This version is preserved and saved as the final submission file.
  2. A **Sanitizer Agent** (or script) intercepts this code and strips out every line starting with `#` and all docstrings (`"""` / `'''`).
  3. The **QA Agent** receives ONLY the bare, sanitized code to evaluate its mathematical logic and error handling blindly.

### 5️⃣ RULE #5: MODULARITY & CLEAN HANDOFFS

* **The Principle:** Agents exist in an automated assembly line. Output must be perfectly primed for the next agent.
* **Execution:** Agents must output ONLY the requested code, plan, or review. Absolutely no conversational filler, greetings, or meta-commentary (e.g., "Here is the updated script!").

### 6️⃣ RULE #6: TASK-SPECIFIC FORMATTING OBEDIENCE

* **The Principle:** Formatting rules are dynamic and dictated by the Architect agent per task, not globally.
* **Execution:** Agents must perfectly mirror the output structure requested in their specific sub-prompt. If instructed to output raw Python logic, they output only Python logic. If asked for a grading matrix by the QA agent, they output only the matrix.

---

## 🤖 AGENT ASSEMBLY LINE ROLES

1. **The Architect Agent** — Breaks project requirements into modular Developer Briefs enforcing the Master Law.
2. **The Executor Agents** — Implement individual offline Python modules with strict inline academic citations.
3. **The Integrator Agent** — Stitches modules into a single crash-proof, self-contained `Group_5_tool.py` (with embedded HTML/CSS/JS UI).
4. **The Sanitizer Agent** — Strips all comments and docstrings into a bare copy for blind QA review.
5. **The QA Agent** — Performs blind adversarial testing (dependency, vulnerability, and math logic checks) on the sanitized copy only.

---

## 🚀 RUNNING THE TOOL

```bash
pip install numpy pandas matplotlib
python3 Group_5_tool.py
```

The tool launches a self-contained local web interface (`http://127.0.0.1:<port>/`) with **zero external files** and **zero internet access required**.

Sample laboratory datasets are provided in [`sample_data/`](sample_data/):
- `example_PD01.csv`
- `sieve_PD01.csv`
- `sieve_reference.csv`
- `flowability_complete_dataset.csv`
