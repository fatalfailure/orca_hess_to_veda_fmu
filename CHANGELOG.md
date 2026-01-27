## [Unreleased]

- Rename `build_mpo` to `adjust_mpo` (wrappers kept for backward compatibility).
- `adjust_mpo` no longer overwrites input `.mpo`; default output is `*_adj.mpo`.
- Increase default max coordination for non-metal heavy atoms to 6 (was 4); metals to 8 (was 6).
- Main converter avoids overwriting existing `.mpo` by writing `*_gen.mpo` when needed.

# Changelog


## 1.0.1 (2026-01-26)
- UI: Moved advanced Hessian options (mass-weighted toggle and unit selection) into a collapsed “Advanced settings” panel to avoid confusing default users.
- Kept advanced options available for troubleshooting uniform frequency scaling issues.

## 1.0.0 (2026-01-26)
- Public release preparation
- English-only UI/messages/comments
- Added English/Japanese manuals and README
- Provided both `.py` (console) and `.pyw` (no-console GUI) entry points
