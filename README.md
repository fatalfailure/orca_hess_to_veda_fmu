# ORCA Hessian → VEDA .fmu Converter

Current release: **v1.1.0**

## Author

- Kaoru Yamamoto (Okayama University of Science)

This repository provides small, dependency-light tools to use **ORCA** vibrational results in **VEDA**.

- `orca_hess_to_veda_fmu.py` / `.pyw`  
  Converts an ORCA `*.out` plus the corresponding `*.hess` into a **VEDA-readable** `*.fmu`
  (a minimal, Gaussian fchk-like text file) and, by default, a helper `*.mpo` connectivity file.

- `adjust_mpo.py` / `.pyw`  
  Standalone **.mpo generator/repair tool**. Use this when you already have `.fmu`/`.xyz`,
  or when connectivity needs manual tuning for metal complexes.

Japanese documentation: see [`README_JA.md`](README_JA.md).

## Requirements

- Python 3.8+ (recommended)
- Standard library only for the main converter GUI (`tkinter`)
- Optional: `tkinterdnd2` for in-window drag-and-drop support
- Optional: `numpy` for the Hessian/frequency sanity check. Conversion itself can still run without it.

## Quick start (GUI)

1. Double-click `orca_hess_to_veda_fmu.pyw` on Windows, or run:
   ```bash
   python orca_hess_to_veda_fmu.py
   ```
2. Select the ORCA `molecule.out` file.
3. Select the corresponding ORCA `molecule.hess` file.
4. Confirm the output `molecule.fmu` path.
5. Click **Convert**.
6. Output files are written next to the selected output path:
   - `molecule.fmu`
   - `molecule.mpo` by default, or `molecule_gen.mpo` if an existing `.mpo` is kept

## What the converter does

1. Reads the ORCA `.out` file and locates the last **VIBRATIONAL FREQUENCIES** section.
2. Selects the nearest Cartesian coordinate block around that frequency section.
3. Reads charge and multiplicity from the lines nearest to the selected coordinate block.
4. Reads coordinates from the `.hess` `$atoms` block when available, because this is usually the most direct coordinate source for the Hessian.
5. Parses the `.hess` `$hessian` block and reconstructs the full Cartesian Hessian matrix.
6. Optionally unweights a mass-weighted Hessian and/or converts Hessian units to Hartree/Bohr².
7. Writes a minimal Gaussian fchk-like `.fmu` file containing atomic numbers, coordinates, and Cartesian force constants.
8. Optionally generates a helper `.mpo` connectivity file using conservative distance-based rules.

## Connectivity / `.mpo` behavior

VEDA can infer connectivity automatically, but metal complexes may be over-connected. This can create too many rings and lead to downstream VEDA problems. The generated `.mpo` therefore uses conservative heuristics:

- non-metal bonds are assigned from covalent radii plus a distance margin;
- metal-ligand bonds are limited to typical donor atoms (`N`, `O`, `S`, `P`, halogens);
- metal coordination is capped by a maximum coordination number and distance cutoff;
- existing `.mpo` files are not overwritten without confirmation.

For difficult cases, regenerate or repair the connectivity with `adjust_mpo.py` / `.pyw`.

## Quick start (VEDA)

VEDA typically prompts for `.fmu` / `.fmt`.

- Provide the generated `molecule.fmu`.
- If VEDA also expects a log-like file (`.fmt`), you can optionally provide ORCA output as `.fmt`
  (for example, copy/rename `molecule.out` to `molecule.fmt`) depending on your workflow.
- If VEDA asks for connectivity, provide the generated `.mpo` and inspect/edit it when needed.

## Files

- `orca_hess_to_veda_fmu.py` — console-capable GUI entry point
- `orca_hess_to_veda_fmu.pyw` — no-console Windows GUI entry point
- `adjust_mpo.py` — console-capable `.mpo` generator/repair tool
- `adjust_mpo.pyw` — no-console Windows GUI entry point
- `docs/MANUAL_EN.md` — detailed English manual
- `docs/MANUAL_JA.md` — 詳細マニュアル（日本語）
- `CHANGELOG.md` — release history
- `CITATION.cff` — citation metadata
- `LICENSE` — MIT License

## Notes / limitations

- The `.fmu` produced here is **not** a Gaussian `.log`. It is a **minimal fchk-like** text container
  that includes atoms, coordinates, and Cartesian force constants (Hessian).
- `.mpo` connectivity is distance-based and heuristic. For challenging systems, especially metal complexes,
  inspect the generated `.mpo` and use `adjust_mpo.py` when necessary.
- The frequency sanity check is intended as a diagnostic, not as a full replacement for ORCA's vibrational analysis.

## Citation

If this program is used for results reported in a manuscript, cite the exact released version used in the calculations. Recommended procedure:

1. Update `CITATION.cff` with the released version, release date, repository URL, and DOI if available.
2. Create a GitHub release tag, for example `v1.1.0`.
3. Archive the release in Zenodo or another repository that assigns a DOI.
4. Cite the archived release DOI in the manuscript and mention the version number in the computational details.

Example manuscript wording:

> ORCA Hessian output was converted to VEDA-readable `.fmu` files using ORCA Hessian → VEDA .fmu Converter v1.1.0 (Yamamoto, 2026).

## Support

- Please report bugs and feature requests via **GitHub Issues**. Public discussion helps other users reproduce and solve similar problems.
- For private inquiries, contact:  
  `k-yamamoto@ous.ac.jp`

## License

MIT License (see `LICENSE`).
