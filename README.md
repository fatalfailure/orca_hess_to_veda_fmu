# ORCA Hessian → VEDA .fmu Converter

This repository provides small, dependency-free tools to use **ORCA** vibrational results in **VEDA**.

- `orca_hess_to_veda_fmu.py` / `.pyw`  
  Converts an ORCA `*.hess` (optionally `*.out`) into a **VEDA-readable** `*.fmu`
  (a minimal, Gaussian fchk-like text file) and optionally a helper `*.mpo`.

- `build_mpo.py` / `.pyw`  
  Standalone **.mpo generator/repair tool** (useful when you already have `.fmu`/`.xyz`,
  or when connectivity needs tuning for metal complexes).

Japanese documentation: see [`README_JA.md`](README_JA.md).

## Requirements
- Python 3.8+ (recommended)
- Standard library only (GUI uses `tkinter`)

## Quick start (GUI)
1. Double-click `orca_hess_to_veda_fmu.pyw` (Windows), or run:
   ```bash
   python orca_hess_to_veda_fmu.py
   ```
2. Select `molecule.hess` (required)
3. (Optional) Select `molecule.out`
4. Click **Convert**
5. Output files are written next to the input:
   - `molecule.fmu`
   - `molecule.mpo` (optional)

## Quick start (VEDA)
VEDA typically prompts for `.fmu` / `.fmt`.

- Provide the generated `molecule.fmu`.
- If VEDA also expects a log-like file (`.fmt`), you can optionally provide ORCA output as `.fmt`
  (e.g., copy/rename `molecule.out` to `molecule.fmt`) depending on your workflow.

## Files
- `docs/MANUAL_EN.md` — detailed English manual
- `docs/MANUAL_JA.md` — 詳細マニュアル（日本語）

## Notes / limitations
- The `.fmu` produced here is **not** a Gaussian `.log`. It is a **minimal fchk-like** text container
  that includes atoms, coordinates, and Cartesian force constants (Hessian).
- `.mpo` connectivity is distance-based and heuristic. For challenging systems (especially metals),
  use `build_mpo.py` to regenerate/repair connectivity.

## License
MIT License (see `LICENSE`).
