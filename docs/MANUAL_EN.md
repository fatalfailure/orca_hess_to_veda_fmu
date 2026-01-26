# Manual (English) — ORCA .hess → VEDA .fmu

## 1. What this tool does
`orca_hess_to_veda_fmu.py` converts an ORCA Hessian file (`*.hess`) into a **VEDA-readable**
`*.fmu` file. The generated `.fmu` is a **minimal Gaussian fchk-like** text container that
includes the data VEDA needs for vibrational analyses:
- number of atoms
- atomic numbers
- Cartesian coordinates
- Cartesian force constants (Hessian; stored in a packed lower-triangle format)

Optionally, the tool generates `*.mpo` (bond connectivity helper).

This tool does **not** generate a Gaussian `.log`/`.out`.

## 2. Inputs
- Required: `molecule.hess`
- Optional: `molecule.out` (used to extract charge and multiplicity when available)

## 3. Outputs
- `molecule.fmu` (always)
- `molecule.mpo` (optional)

Both are written next to the input file.

## 4. Usage
### GUI
Run:
- Windows: double-click `orca_hess_to_veda_fmu.pyw`
- Or:
  ```bash
  python orca_hess_to_veda_fmu.py
  ```

### VEDA
VEDA prompts for `.fmu`/`.fmt`.
- Select the generated `.fmu`.
- If your workflow requires `.fmt`, place a suitable log-like file next to it.
  For ORCA, you can often use `molecule.out` copied/renamed to `molecule.fmt`.

## 5. Connectivity (.mpo)
The built-in `.mpo` generator is distance-based (covalent radii heuristics). It works well
for many organic systems, but may require tuning for metal complexes.

If connectivity is wrong:
- Use `build_mpo.py` to regenerate/repair `.mpo`.
- Adjust thresholds in the GUI (for `build_mpo.py`) if needed.

## 6. Troubleshooting
- **VEDA refuses the file**: check file paths, file permissions, and that the `.fmu` is not empty.
- **Unexpected frequency scaling**: indicates a likely mismatch in Hessian units or assumptions.
  Ensure `.hess` corresponds to the intended optimized structure.
