#! python3
# -*- coding: utf-8 -*-
# orca_hess_to_veda_fmu.pyw — ORCA (.out/.hess) → for VEDA minimal fchk (.fmu)
#
# Specification overview:
# - Convert ORCA .out/.hess into a Gaussian fchk-like .fmu readable by VEDA
# - Works even if frequencies appear mid-file: anchor on the last frequency header
#   and select the nearest coordinate block around it
# - Read charge/multiplicity from the lines nearest to the coordinate-block header
# (Comment translated to English for public release)
# - Prefer coordinates from the .hess $atoms block when available
# - Parse ORCA $hessian and reconstruct a full dim x dim matrix
# - Output is ASCII with CRLF line endings; arrays use Gaussian fchk-style "N=" headers

import sys, os, re, math, traceback
from pathlib import Path

# optional in-window DnD
_has_tkdnd = False
try:
    from tkinterdnd2 import DND_FILES, TkinterDnD
    _has_tkdnd = True
except Exception:
    _has_tkdnd = False

# ---------- GUI ----------
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, simpledialog

APP_TITLE = "ORCA -> VEDA .fmu Converter (GUI)"
DEFAULT_UNITS = "hartree_bohr2"
UNIT_CHOICES = ("hartree_bohr2", "mdyn_a2", "j_m2")

# ---------- chemistry basics ----------
SYMBOL_TO_Z = {
    "H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
    "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,
    "K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,"Cu":29,"Zn":30,
    "Ga":31,"Ge":32,"As":33,"Se":34,"Br":35,"Kr":36,
    "Rb":37,"Sr":38,"Y":39,"Zr":40,"Nb":41,"Mo":42,"Tc":43,"Ru":44,"Rh":45,"Pd":46,"Ag":47,"Cd":48,
    "In":49,"Sn":50,"Sb":51,"Te":52,"I":53,"Xe":54,
}
# Reverse lookup (atomic number -> symbol)
Z_TO_SYMBOL = {v: k for k, v in SYMBOL_TO_Z.items()}

# Covalent radii in Å (mostly Pyykkö 2009 style values; coarse is fine for connectivity heuristics)
# Only elements likely to appear in typical ORCA/DFT jobs are listed; others fall back to a safe default.
COVALENT_RADIUS_A = {
    "H":0.31,"B":0.85,"C":0.76,"N":0.71,"O":0.66,"F":0.57,
    "Si":1.11,"P":1.07,"S":1.05,"Cl":1.02,"Br":1.20,"I":1.39,
    "Li":1.28,"Na":1.66,"K":2.03,"Mg":1.41,"Ca":1.76,
    "Al":1.21,"Ga":1.22,"In":1.42,"Sn":1.39,
    "Fe":1.24,"Co":1.18,"Ni":1.17,"Cu":1.22,"Zn":1.20,
    "Mn":1.39,"Cr":1.28,
}

# --- MPO generation heuristics ---
# VEDA can *guess* MPO from geometry, but for metal complexes it often over-connects (ring explosion → DD2 shortage).
# Here we generate a conservative MPO alongside .fmu.
DEFAULT_DONORS = {"N", "O", "S", "P", "F", "Cl", "Br", "I"}

# Metal-specific rules: only bond to donor atoms, keep the nearest max_coord within cutoff (Å)
# Metal-specific rules: only bond to donor atoms, keep the nearest max_coord within cutoff (Å)
#
# Rationale:
# - VEDA's automatic connectivity inference can over-connect around metals.
# - We therefore generate a conservative MPO: metals only bond to typical donor atoms (N/O/S/P/halogens)
#   and we keep only the nearest neighbors up to max_coord.
# - If a metal you use is not listed, it will be treated as a normal element (non-metal pass),
#   and the 'metal' coordination limits will not apply.
DEFAULT_METAL_MAX_COORD = 6
DEFAULT_METAL_CUTOFF_A = 2.60

METAL_SYMBOLS = {
    # 3d series (common)
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    # 4d series (common)
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    # 5d series (common)
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
    # main-group metals often seen in coordination chemistry
    "Al","Ga","In","Sn","Pb","Bi",
}

METAL_RULES = {
    sym: {"donors": DEFAULT_DONORS, "max_coord": DEFAULT_METAL_MAX_COORD, "cutoff": DEFAULT_METAL_CUTOFF_A}
    for sym in METAL_SYMBOLS
}

def _sym_from_Z(z: int) -> str:
    return Z_TO_SYMBOL.get(int(z), "X")

def _covrad(sym: str) -> float:
    return COVALENT_RADIUS_A.get(sym, 0.77)  # fallback ~C

def _dist(a, b) -> float:
    dx = a[0]-b[0]; dy = a[1]-b[1]; dz = a[2]-b[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def build_connectivity(symbols, coords_ang,
                       scale: float = 1.20,
                       add_margin: float = 0.25):
    """Return adjacency list (1-based indices) for MPO.

    - Non-metals: bond if d <= scale*(r_i+r_j)+add_margin.
    - Metals in METAL_RULES: bond only to donors and keep nearest max_coord within cutoff.
    """
    n = len(symbols)
    adj = [set() for _ in range(n)]

    # --- First: non-metal bonds (including ligand framework) ---
    metals = set(METAL_RULES.keys())
    for i in range(n):
        si = symbols[i]
        for j in range(i+1, n):
            sj = symbols[j]
            # defer any metal-involving pair to metal pass (prevents metal-C/H overbonding)
            if si in metals or sj in metals:
                continue
            d = _dist(coords_ang[i], coords_ang[j])
            thr = scale * (_covrad(si) + _covrad(sj)) + add_margin
            if d <= thr:
                adj[i].add(j+1); adj[j].add(i+1)

    # --- Second: metal bonds (conservative) ---
    for i in range(n):
        si = symbols[i]
        if si not in METAL_RULES:
            continue
        rule = METAL_RULES[si]
        donors = rule["donors"]
        cutoff = float(rule["cutoff"])
        maxc = int(rule["max_coord"])

        cand = []
        for j in range(n):
            if j == i:
                continue
            sj = symbols[j]
            if sj not in donors:
                continue
            d = _dist(coords_ang[i], coords_ang[j])
            if d <= cutoff:
                cand.append((d, j))
        cand.sort(key=lambda x: x[0])
        keep = cand[:maxc]

        for d, j in keep:
            adj[i].add(j+1); adj[j].add(i+1)

    return adj

def write_mpo(path: Path, Z_list, coords_ang):
    """Write a VEDA-style .mpo based on conservative connectivity."""
    symbols = [_sym_from_Z(z) for z in Z_list]
    adj = build_connectivity(symbols, coords_ang)

    lines = []
    # Header is optional; keep minimal and VEDA-friendly.
    for i, neigh in enumerate(adj, start=1):
        sym_i = symbols[i-1]
        neigh_sorted = sorted(neigh)
        parts = []
        for j in neigh_sorted:
            parts.append(f"{symbols[j-1]} {j}")
        deg = len(parts)
        if parts:
            line = f"{deg:2d}. {sym_i} {i}: " + ", ".join(parts)
        else:
            line = f"{deg:2d}. {sym_i} {i}:"
        lines.append(line)

    # Use CRLF for Windows friendliness (VEDA ecosystem)
    path.write_text("\r\n".join(lines) + "\r\n", encoding="utf-8")

Z_TO_MASS = {  # amu
    1:1.00784, 2:4.002602, 3:6.94, 4:9.0121831, 5:10.81, 6:12.011, 7:14.0067, 8:15.999,
    9:18.998403, 10:20.1797, 11:22.989769, 12:24.305, 13:26.981538, 14:28.085, 15:30.973762,
    16:32.06, 17:35.45, 18:39.948, 19:39.0983, 20:40.078, 26:55.845, 29:63.546, 30:65.38,
    35:79.904, 53:126.90447
}
AMU_TO_ME = 1822.888486209  # m_e per amu

# ---------- units ----------
def unit_factor_to_hartree_bohr2(units: str) -> float:
    u = (units or "").lower()
    if u in ("hartree_bohr2","hartree/bohr^2","hartree_bohr^2","hartree/bohr2"):
        return 1.0
    HARTREE = 4.3597447222071e-18  # J
    BOHR = 5.29177210903e-11       # m
    if u in ("j_m2","j/m^2","joule/m^2"):
        return 1.0 / HARTREE * (BOHR**2)
    if u in ("mdyn_a2","mdyn/ang^2","mdyne/ang^2","mdyn/å^2"):
        return unit_factor_to_hartree_bohr2("j_m2") * 1e12  # 1 mdyn/Å^2 = 1e12 J/m^2
    raise ValueError(f"Unsupported units: {units}")

# ---------- fchk writer ----------
def pack_lower_triangular(H):
    dim = len(H)
    data = []
    for i in range(dim):
        for j in range(i+1):
            data.append(H[i][j])
    return data

def write_minimal_fchk(path: Path, Z_list, coords_ang, charge, mult, H_full):
    """Write .fmu in a minimal Gaussian fchk-compatible format (ASCII + CRLF)."""
    nat = len(Z_list)
    dim = 3*nat
    H_pack = pack_lower_triangular(H_full)

    def w_count(f, label, typ, count):
        # (Comment translated to English for public release)
        f.write(f"{label:<40s} {typ:>1s}  N={count:>12d}\r\n")
    def w_scalar(f, label, typ, val):
        f.write(f"{label:<40s} {typ:>1s} {val:>12d}\r\n")
    def w_arrR(f, data, per_line=5):
        for i in range(0, len(data), per_line):
            chunk = data[i:i+per_line]
            f.write(" ".join(f"{x:>16.8f}" for x in chunk) + "\r\n")
    def w_arrI(f, data, per_line=6):
        for i in range(0, len(data), per_line):
            chunk = data[i:i+per_line]
            f.write(" ".join(f"{int(x):>12d}" for x in chunk) + "\r\n")

    # (Comment translated to English for public release)
    with path.open("w", newline="", encoding="ascii") as f:
        f.write("Generated by orca_hess_to_veda_fmu.pyw (minimal fchk for VEDA)\r\n")
        # (Comment translated to English for public release)
        w_scalar(f, "Number of atoms", "I", int(nat))
        w_scalar(f, "Charge", "I", int(charge))
        w_scalar(f, "Multiplicity", "I", int(mult))
        # (Comment translated to English for public release)
        w_count(f, "Atomic numbers", "I", nat)
        w_arrI(f, Z_list)
        flat = []
        for (x,y,z) in coords_ang:
            flat.extend([x,y,z])
        w_count(f, "Current cartesian coordinates", "R", 3*nat)
        w_arrR(f, flat)
        w_count(f, "Cartesian Force Constants", "R", dim*(dim+1)//2)
        w_arrR(f, H_pack)

# ---------- helpers for parsing ----------
_HDR_COORD = re.compile(
    r"^(?P<hdr>(?:CARTESIAN\s+COORDINATES|COORDINATES)\s*\((?P<unit>ANGSTROEM|ANGSTROM|A\.U\.|AU|BOHR)\))\s*$",
    re.I | re.M
)
_HDR_FREQ = re.compile(r"VIBRATIONAL\s+FREQUENC(?:IES|Y)", re.I)

def _normalize_out_text(txt: str) -> str:
    """Remove ORCA line-number prefixes (e.g., '| 51> ...')."""
    return re.sub(r'(?m)^\s*\|\s*\d+>\s*', '', txt)

def _find_freq_anchor(txt: str) -> int:
    """Position of the last frequency header (or len(txt) if not found)."""
    ms = list(_HDR_FREQ.finditer(txt))
    return ms[-1].start() if ms else len(txt)

def _scan_coord_table(txt: str, start_pos: int) -> str:
    """Extract the coordinate block following a header."""
    lines = txt[start_pos:].splitlines()
    out = []
    started = False

    def looks_like_coord_line(s: str) -> bool:
        s = s.strip()
        if not s:
            return False
        if re.match(r"^(Atom|Nr\.?|No\.?|x\s+y\s+z|[-=]{3,})$", s, re.I):
            return False
        toks = s.replace(",", " ").split()
        floats = sum(1 for t in toks if _is_float(t))
        return floats >= 3

    for ln in lines:
        s = ln.rstrip("\r\n")
        if not s.strip():
            if started:
                break
            else:
                continue
        if looks_like_coord_line(s):
            started = True
            out.append(s)
        else:
            if started:
                break
            continue
    return "\n".join(out)

def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except Exception:
        return False

def _token_to_atom(toks):
    """Parse various coordinate line formats into (Z, x, y, z)."""
    p = toks
    def F(i): return float(p[i])

    # A: Sym x y z
    if len(p) >= 4 and p[0].isalpha() and all(_is_float(t) for t in p[1:4]):
        sym = p[0][0] + p[0][1:].lower()
        if sym in SYMBOL_TO_Z:
            return SYMBOL_TO_Z[sym], F(1), F(2), F(3)
    # B: idx Sym x y z
    if len(p) >= 5 and p[0].isdigit() and p[1].isalpha() and all(_is_float(t) for t in p[2:5]):
        sym = p[1][0] + p[1][1:].lower()
        if sym in SYMBOL_TO_Z:
            return SYMBOL_TO_Z[sym], F(2), F(3), F(4)
    # C: idx Sym Z x y z
    if len(p) >= 6 and p[0].isdigit() and p[1].isalpha() and p[2].isdigit() and all(_is_float(t) for t in p[3:6]):
        return int(p[2]), F(3), F(4), F(5)
    # D: idx Z x y z
    if len(p) >= 5 and p[0].isdigit() and p[1].isdigit() and all(_is_float(t) for t in p[2:5]):
        return int(p[1]), F(2), F(3), F(4)
    # E: Sym Z x y z
    if len(p) >= 5 and p[0].isalpha() and p[1].isdigit() and all(_is_float(t) for t in p[2:5]):
        return int(p[1]), F(2), F(3), F(4)
    return None

def _find_nearest_coord_block(txt: str, anchor: int):
    """Select the coordinate header nearest to the frequency anchor and return (unit, block, header_pos)."""
    before = None
    after = None
    for m in _HDR_COORD.finditer(txt):
        if m.start() < anchor:
            before = m
        elif after is None:
            after = m
    if not before and not after:
        raise ValueError("Coordinates block header not found around frequency section")

    def pick(m):
        unit = m.group("unit").upper()
        block = _scan_coord_table(txt, m.end())
        return unit, block, m.start()

    if before and after:
        d_before = anchor - before.start()
        d_after = after.start() - anchor
        return pick(before) if d_before <= d_after else pick(after)
    if before:
        return pick(before)
    return pick(after)

# --- Charge/Multiplicity detector (header-pos based) ---
_RE_COUPLED = re.compile(r"Charge\s*=\s*([-\d]+)\s+Multiplicity\s*=\s*([-\d]+)", re.I)
_RE_STAR = re.compile(
    r"^\s*(?:\|\s*\d+>\s*)?\*\s+(?:xyz|xyzfile|int|zmat|zmatfile)\s+([-\d]+)\s+([-\d]+)\b",
    re.I | re.M
)
_RE_CHARGE1 = re.compile(r"Total\s+Charge\s*[:=]\s*([-\d]+)", re.I)
_RE_CHARGE2 = re.compile(r"Overall\s+charge.*?([-\d]+)", re.I)
_RE_CHARGE3 = re.compile(r"^\s*Charge\s*[:=]?\s*([-\d]+)\s*$", re.I | re.M)
_RE_MULT1 = re.compile(r"Multiplicity\s*(?:\(2S\+1\))?\s*[:=]\s*([-\d]+)", re.I)
_RE_MULT2 = re.compile(r"Spin\s+multiplicity\s*[:=]\s*([-\d]+)", re.I)
_RE_MULT3 = re.compile(r"^\s*Multiplicity\s*[:=]?\s*([-\d]+)\s*$", re.I | re.M)

def _find_charge_mult_near(txt: str, target_pos: int):
    """Return the charge/multiplicity closest to the coordinate header position target_pos."""
    # (Comment translated to English for public release)
    pairs = []
    for m in _RE_COUPLED.finditer(txt):
        pairs.append((abs(m.start() - target_pos), int(m.group(1)), int(m.group(2))))
    for m in _RE_STAR.finditer(txt):
        pairs.append((abs(m.start() - target_pos), int(m.group(1)), int(m.group(2))))
    if pairs:
        pairs.sort(key=lambda x: x[0])
        return pairs[0][1], pairs[0][2]

    # (Comment translated to English for public release)
    charges = []
    for rx in (_RE_CHARGE1, _RE_CHARGE2, _RE_CHARGE3):
        for m in rx.finditer(txt):
            charges.append((abs(m.start() - target_pos), int(m.group(1))))
    mults = []
    for rx in (_RE_MULT1, _RE_MULT2, _RE_MULT3):
        for m in rx.finditer(txt):
            mults.append((abs(m.start() - target_pos), int(m.group(1))))

    charge = min(charges, key=lambda x: x[0])[1] if charges else None
    mult = min(mults, key=lambda x: x[0])[1] if mults else None
    return charge, mult

# ---------- .out robust parser ----------
def parse_out_simple(out_path: Path):
    """Using the last frequency header as an anchor, get the nearest coordinate block and nearby charge/multiplicity."""
    with open(out_path, "r", encoding="utf-8", errors="ignore") as fh:
        txt = fh.read()
    txt = _normalize_out_text(txt)

    anchor = _find_freq_anchor(txt)
    unit, block, hdr_pos = _find_nearest_coord_block(txt, anchor)

    # Charge / Multiplicity
    charge, mult = _find_charge_mult_near(txt, hdr_pos)

    # coordinates
    Z_list = []
    coords = []
    for ln in block.splitlines():
        s = ln.strip()
        if not s or re.match(r"^(Atom|Nr\.?|No\.?|x\s+y\s+z|[-=]{3,})$", s, re.I):
            continue
        toks = s.replace(",", " ").split()
        got = _token_to_atom(toks)
        if got is None:
            continue
        znum, x, y, z = got
        if 1 <= int(znum) <= 118:
            Z_list.append(int(znum))
            coords.append((x, y, z))

    if not Z_list:
        raise ValueError("Failed to parse any atoms in coordinates block (nearest to freq section)")

    # BOHR/AU → Å
    if unit in ("BOHR", "A.U.", "AU"):
        BOHR = 0.529177210903
        coords = [(x * BOHR, y * BOHR, z * BOHR) for (x, y, z) in coords]

    return Z_list, coords, charge, mult

# ---------- .hess: coordinates from $atoms ----------


def parse_out_frequencies(out_path: Path):
    """
    Extract harmonic frequencies (cm^-1) from the last 'VIBRATIONAL FREQUENCIES' block in ORCA .out.

    Returns [] if not found.
    """
    try:
        txt = out_path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return []
    txt = _normalize_out_text(txt)

    # (Comment translated to English for public release)
    anchors = [m.start() for m in _HDR_FREQ.finditer(txt)]
    if not anchors:
        return []
    a = anchors[-1]
    tail = txt[a:]
    lines = tail.splitlines()

    freqs = []
    # (Comment translated to English for public release)
    # (Comment translated to English for public release)
    #   1      100.23      0.00 ...
    #   2     1500.12      0.00 ...
    num_re = re.compile(r"^[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[Ee][+-]?\d+)?$")

    # (Comment translated to English for public release)
    for line in lines[1:800]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("----") or s.lower().startswith("mode"):
            continue
        # (Comment translated to English for public release)
        if _HDR_COORD.match(s) or s.upper().startswith("NORMAL MODES") or s.startswith("$"):
            break

        toks = s.split()
        if len(toks) < 2:
            continue
        # (Comment translated to English for public release)
        if re.fullmatch(r"\d+", toks[0]) and num_re.match(toks[1]):
            try:
                f = float(toks[1])
            except Exception:
                continue
            freqs.append(f)

    return freqs


def compute_freqs_from_hessian_cm1(H_cart, Z_list, assume_massweighted_input: bool):
    """
    Roughly recompute vibrational frequencies (cm^-1) from the Hessian for sanity checking.

    - H_cart: (3N x 3N) Hessian. Units assumed Hartree/Bohr^2.
    - If assume_massweighted_input=True: interpret H_cart as already mass-weighted and
      do not apply additional mass-weighting.

    Note: translational/rotational projection is approximate; we drop the smallest 6
    eigenvalues (5 for linear molecules) as a simple treatment.
    """
    try:
        import numpy as _np
    except Exception as e:
        raise RuntimeError("numpy is required for frequency check but is not available") from e

    nat = len(Z_list)
    dim = 3 * nat
    if len(H_cart) != dim or len(H_cart[0]) != dim:
        raise ValueError("Hessian dimension mismatch")

    H = _np.array(H_cart, dtype=float)

    # (Comment translated to English for public release)
    m_amu = _np.array([Z_TO_MASS.get(z, float(z)) for z in Z_list], dtype=float)
    m_me = m_amu * AMU_TO_ME
    m3 = _np.repeat(m_me, 3)

    if assume_massweighted_input:
        Hw = H
    else:
        inv_sqrt = 1.0 / _np.sqrt(m3)
        Hw = (inv_sqrt[:, None] * H) * inv_sqrt[None, :]

    # (Comment translated to English for public release)
    Hw = 0.5 * (Hw + Hw.T)

    # (Comment translated to English for public release)
    w2, _ = _np.linalg.eigh(Hw)

    # (Comment translated to English for public release)
    HARTREE = 4.3597447222071e-18  # J
    ME = 9.1093837015e-31          # kg
    BOHR = 5.29177210903e-11       # m
    c_cm = 2.99792458e10           # cm/s
    K = math.sqrt(HARTREE / (ME * BOHR * BOHR)) / (2.0 * math.pi * c_cm)

    freqs = []
    for lam in w2:
        # (Comment translated to English for public release)
        if lam >= 0:
            freqs.append(math.sqrt(lam) * K)
        else:
            freqs.append(-math.sqrt(-lam) * K)

    freqs_sorted = sorted(freqs, key=lambda x: abs(x))

    # (Comment translated to English for public release)
    drop = 6
    # (Comment translated to English for public release)
    cand5 = freqs_sorted[5:]
    cand6 = freqs_sorted[6:]
    # (Comment translated to English for public release)
    # (Comment translated to English for public release)
    use = cand6 if len(cand6) >= max(0, len(cand5) - 1) else cand5

    return sorted(use)

def _rms_diff(a, b, n=20):
    if not a or not b:
        return None
    m = min(len(a), len(b), n)
    if m <= 0:
        return None
    s = 0.0
    for i in range(m):
        d = a[i] - b[i]
        s += d * d
    return math.sqrt(s / m)
def parse_atoms_from_hess(hess_path: Path):
    """Read (Z_list, coords [Angstrom]) from the $atoms block."""
    with open(hess_path, "r", encoding="utf-8", errors="ignore") as fh:
        txt = fh.read()

    m = re.search(r"^\s*\$atoms\s*$", txt, flags=re.M | re.I)
    if not m:
        raise ValueError("$atoms block not found in .hess")

    rest = txt[m.end():]
    lines = rest.splitlines()

    i = 0
    # (Comment translated to English for public release)
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines):
        raise ValueError("No atom count after $atoms")

    # (Comment translated to English for public release)
    try:
        nat = int(lines[i].split()[0])
    except Exception as e:
        raise ValueError(f"Failed to read atom count after $atoms: {e}")
    i += 1

    BOHR = 0.529177210903
    Z_list = []
    coords_ang = []

    for k in range(nat):
        if i + k >= len(lines):
            raise ValueError("Not enough atom lines in $atoms block")
        parts = lines[i + k].split()
        if len(parts) < 5:
            raise ValueError(f"Bad $atoms line: {lines[i+k]!r}")
        sym = parts[0]
        x_bohr = float(parts[2])
        y_bohr = float(parts[3])
        z_bohr = float(parts[4])

        Z = SYMBOL_TO_Z.get(sym.capitalize())
        if Z is None:
            raise ValueError(f"Unknown element symbol in $atoms: {sym!r}")
        Z_list.append(Z)
        coords_ang.append((x_bohr * BOHR, y_bohr * BOHR, z_bohr * BOHR))

    return Z_list, coords_ang

# (Comment translated to English for public release)
def parse_hess_simple(hess_path: Path, natoms: int):
    """
    Robustly parse the $hessian block in an ORCA .hess file.

    Typical format:
      $hessian
      dim
      <col header>    0    1    2    3    4
      <row>  v00  v01  v02  v03  v04
      ...

    Strategy:
    - Detect integer-only column-header lines and map values to those columns.
    - If headers are missing/broken, fall back to filling values left-to-right.
    - Raise ValueError if data are missing (fail loudly rather than silently).
    - Symmetrize at the end using (H + H^T)/2.
    """
    with open(hess_path, "r", encoding="utf-8", errors="ignore") as fh:
        txt = fh.read()

    m = re.search(r"^\s*\$hessian\s*$", txt, flags=re.M | re.I)
    if not m:
        raise ValueError("$hessian block not found in .hess")

    rest = txt[m.end():]
    lines = rest.splitlines()

    # (Comment translated to English for public release)
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines):
        raise ValueError("No dimension line after $hessian")

    try:
        dim = int(lines[i].split()[0])
    except Exception as e:
        raise ValueError(f"Failed to read dimension after $hessian: {e}")
    i += 1

    # (Comment translated to English for public release)
    if dim % 3 == 0:
        inferred_nat = dim // 3
        if natoms <= 0 or natoms != inferred_nat:
            natoms = inferred_nat
    elif not (natoms > 0 and 3 * natoms == dim):
        raise ValueError(f"Dimension {dim} is not a multiple of 3; cannot infer natoms properly")

    H = [[0.0] * dim for _ in range(dim)]

    # (Comment translated to English for public release)
    current_cols = None  # type: list[int] | None

    # (Comment translated to English for public release)
    row_fill = [0] * dim

    seen = set()  # (row, col)

    def _is_int_tok(t: str) -> bool:
        return re.fullmatch(r"[+-]?\d+", t) is not None

    for line in lines[i:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("$"):
            break

        tokens = line.split()
        if not tokens:
            continue

        # (Comment translated to English for public release)
        if all(_is_int_tok(t) for t in tokens):
            # (Comment translated to English for public release)
            cols = [int(t) for t in tokens]
            # (Comment translated to English for public release)
            cols = [c for c in cols if 0 <= c < dim]
            current_cols = cols if cols else None
            continue

        # (Comment translated to English for public release)
        try:
            row_idx = int(tokens[0])
        except Exception:
            continue
        if not (0 <= row_idx < dim):
            continue

        values = tokens[1:]
        if not values:
            continue

        if current_cols is not None and len(current_cols) >= 1:
            # (Comment translated to English for public release)
            # (Comment translated to English for public release)
            n = min(len(values), len(current_cols))
            for k in range(n):
                col_idx = current_cols[k]
                try:
                    H[row_idx][col_idx] = float(values[k])
                    seen.add((row_idx, col_idx))
                except Exception:
                    pass
        else:
            # (Comment translated to English for public release)
            col = row_fill[row_idx]
            for v in values:
                if col >= dim:
                    break
                try:
                    H[row_idx][col] = float(v)
                    seen.add((row_idx, col))
                except Exception:
                    pass
                col += 1
            row_fill[row_idx] = col

    # (Comment translated to English for public release)
    # (Comment translated to English for public release)
    expected_full = dim * dim
    expected_lower = dim * (dim + 1) // 2
    # (Comment translated to English for public release)
    if len(seen) < expected_lower * 0.90:
        raise ValueError(
            f"Hessian parse seems incomplete: filled {len(seen)}/{expected_lower} (lower-tri expected). "
            "(.hess format mismatch? column headers not handled?)"
        )
    # (Comment translated to English for public release)
    for ir in range(dim):
        for jc in range(ir + 1, dim):
            v = 0.5 * (H[ir][jc] + H[jc][ir])
            H[ir][jc] = H[jc][ir] = v

    return H

def unweight_massweighted(H_mw, Z_list):
    dim = len(H_mw)
    nat = len(Z_list)
    assert dim == 3 * nat
    ms = []
    for z in Z_list:
        m = Z_TO_MASS.get(z, float(z)) * AMU_TO_ME
        s = math.sqrt(m)
        ms.extend([s, s, s])
    H = [[0.0] * dim for _ in range(dim)]
    for i in range(dim):
        si = ms[i]
        row = H_mw[i]
        for j in range(dim):
            H[i][j] = si * row[j] * ms[j]
    return H

# ---------- App ----------
class App:
    def __init__(self, root):
        self.root = root
        root.title(APP_TITLE)
        root.geometry("820x560")
        root.minsize(760, 500)

        main = ttk.Frame(root, padding=10)
        main.pack(fill="both", expand=True)

        self.var_out = tk.StringVar()
        self.var_hess = tk.StringVar()
        self.var_fmu = tk.StringVar()

        row = 0
        ttk.Label(main, text=".out file").grid(row=row, column=0, sticky="e")
        e_out = ttk.Entry(main, textvariable=self.var_out, width=86)
        e_out.grid(row=row, column=1, sticky="we", padx=6)
        ttk.Button(main, text="Browse...", command=self.browse_out).grid(row=row, column=2, padx=2)
        row += 1

        ttk.Label(main, text=".hess file").grid(row=row, column=0, sticky="e")
        e_hess = ttk.Entry(main, textvariable=self.var_hess, width=86)
        e_hess.grid(row=row, column=1, sticky="we", padx=6)
        ttk.Button(main, text="Browse...", command=self.browse_hess).grid(row=row, column=2, padx=2)
        row += 1

        ttk.Label(main, text="Output .fmu").grid(row=row, column=0, sticky="e")
        e_fmu = ttk.Entry(main, textvariable=self.var_fmu, width=86)
        e_fmu.grid(row=row, column=1, sticky="we", padx=6)
        ttk.Button(main, text="Change...", command=self.browse_fmu).grid(row=row, column=2, padx=2)
        row += 1

        # Advanced settings (collapsed by default)
        # Most users should not need to change these.
        self.var_mass = tk.BooleanVar(value=False)
        self.var_units = tk.StringVar(value=DEFAULT_UNITS)

        self._adv_open = tk.BooleanVar(value=False)
        self.adv_button = ttk.Button(main, text="Advanced settings ▶", command=self.toggle_advanced)
        self.adv_button.grid(row=row, column=0, columnspan=3, sticky="w", pady=(8, 2))
        row += 1

        self._adv_row = row

        self.adv_frame = ttk.Frame(main)
        # Widgets inside the advanced frame
        ttk.Checkbutton(
            self.adv_frame,
            text="Hessian is mass-weighted (unweight before writing .fmu)",
            variable=self.var_mass,
        ).grid(row=0, column=0, columnspan=2, sticky="w", pady=(0, 6))

        ttk.Label(self.adv_frame, text="Hessian input unit").grid(row=1, column=0, sticky="e", padx=(0, 6))
        self.adv_units_cmb = ttk.Combobox(
            self.adv_frame,
            textvariable=self.var_units,
            values=UNIT_CHOICES,
            width=20,
            state="readonly",
        )
        self.adv_units_cmb.grid(row=1, column=1, sticky="w")

        ttk.Label(
            self.adv_frame,
            text="Tip: Defaults work for standard ORCA .hess files. Change only if frequencies are uniformly scaled.",
            foreground="#555555",
        ).grid(row=2, column=0, columnspan=2, sticky="w", pady=(6, 0))


        self.var_gen_mpo = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            main,
            text="Also generate helper .mpo (connectivity) [recommended]. For some metal complexes, you can further refine connectivity with adjust_mpo.",
            variable=self.var_gen_mpo,
        ).grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 6))
        row += 1

        btns = ttk.Frame(main)
        btns.grid(row=row, column=0, columnspan=3, sticky="we", pady=8)
        ttk.Button(btns, text="Convert", command=self.convert).pack(side="left")
        ttk.Button(btns, text="Close", command=root.destroy).pack(side="right")
        row += 1

        self.log = tk.Text(main, height=16)
        self.log.grid(row=row, column=0, columnspan=3, sticky="nsew")
        main.rowconfigure(row, weight=1)
        main.columnconfigure(1, weight=1)

        if _has_tkdnd:
            for widget in (e_out, e_hess, e_fmu, self.log):
                widget.drop_target_register(DND_FILES)
                widget.dnd_bind("<<Drop>>", self.on_drop)

        self.consume_argv()

    def logprint(self, *a):
        msg = " ".join(str(x) for x in a) + "\n"
        self.log.insert("end", msg)
        self.log.see("end")
        self.root.update_idletasks()

    def on_drop(self, event):
        paths = self._parse_dnd_paths(event.data)
        for p in paths:
            self._accept_path(Path(p))

    def _parse_dnd_paths(self, data: str):
        out = []
        cur = ""
        in_brace = False
        for ch in data:
            if ch == "{":
                in_brace = True
                cur = ""
            elif ch == "}":
                in_brace = False
                out.append(cur)
                cur = ""
            elif ch == " " and not in_brace:
                if cur:
                    out.append(cur)
                    cur = ""
            else:
                cur += ch
        if cur:
            out.append(cur)
        return out

    def _accept_path(self, p: Path):
        if not p.exists():
            return
        sfx = p.suffix.lower()
        if sfx in (".out", ".log"):
            self.var_out.set(str(p))
            hess = p.with_suffix(".hess")
            if hess.exists():
                self.var_hess.set(str(hess))
            self._suggest_fmu()
        elif sfx == ".hess":
            self.var_hess.set(str(p))
            out = p.with_suffix(".out")
            if out.exists():
                self.var_out.set(str(out))
            self._suggest_fmu()
        elif sfx in (".fmu", ".fchk", ".fch"):
            self.var_fmu.set(str(p))
        else:
            self.logprint("Ignored:", p)

    def _suggest_fmu(self):
        base = None
        if self.var_out.get():
            base = Path(self.var_out.get()).with_suffix(".fmu")
        elif self.var_hess.get():
            base = Path(self.var_hess.get()).with_suffix(".fmu")
        if base:
            self.var_fmu.set(str(base))

    def consume_argv(self):
        argv = sys.argv[1:]
        if not argv:
            return
        for a in argv:
            self._accept_path(Path(a))

    def browse_out(self):
        path = filedialog.askopenfilename(
            title="Select .out",
            filetypes=[("ORCA out/log", "*.out *.log"), ("All", "*.*")],
        )
        if path:
            self.var_out.set(path)
            h = Path(path).with_suffix(".hess")
            if h.exists():
                self.var_hess.set(str(h))
            self._suggest_fmu()

    def browse_hess(self):
        path = filedialog.askopenfilename(
            title="Select .hess",
            filetypes=[("ORCA hess", "*.hess"), ("All", "*.*")],
        )
        if path:
            self.var_hess.set(path)
            o = Path(path).with_suffix(".out")
            if o.exists():
                self.var_out.set(str(o))
            self._suggest_fmu()

    
    def toggle_advanced(self):
        """Show/hide the advanced settings panel."""
        if not self._adv_open.get():
            # Open
            self.adv_frame.grid(row=self._adv_row, column=0, columnspan=3, sticky="w", pady=(0, 10))
            self.adv_button.config(text="Advanced settings ▼")
            self._adv_open.set(True)
        else:
            # Close
            self.adv_frame.grid_forget()
            self.adv_button.config(text="Advanced settings ▶")
            self._adv_open.set(False)

    def browse_fmu(self):
        initial = self.var_fmu.get() or (self.var_out.get() or self.var_hess.get())
        if initial:
            initial = str(Path(initial).with_suffix(".fmu"))
        path = filedialog.asksaveasfilename(
            title="Select output .fmu",
            defaultextension=".fmu",
            filetypes=[
                ("VEDA fmu (formatted checkpoint)", "*.fmu"),
                ("All", "*.*"),
            ],
            initialfile=os.path.basename(initial) if initial else None,
            initialdir=os.path.dirname(initial) if initial else None,
        )
        if path:
            self.var_fmu.set(path)

    def convert(self):
        try:
            self.logprint("python exe:", sys.executable)

            outp = Path(self.var_out.get())
            hessp = Path(self.var_hess.get())
            fmup = Path(self.var_fmu.get()) if self.var_fmu.get() else None
            if not outp.exists():
                raise FileNotFoundError(".out file not found.")
            if not hessp.exists():
                raise FileNotFoundError(".hess file not found.")
            if not fmup:
                fmup = outp.with_suffix(".fmu")
                self.var_fmu.set(str(fmup))
            units = self.var_units.get()
            massw = self.var_mass.get()

            self.logprint("=== Conversion start ===")
            self.logprint("out :", outp)
            self.logprint("hess:", hessp)
            self.logprint("fmu :", fmup)
            self.logprint("units:", units, " mass-weighted?:", massw)

            # (Comment translated to English for public release)
            Z_list, coords, charge, mult = parse_out_simple(outp)
            self.logprint(
                f".out: natoms={len(Z_list)}  charge={charge}  multiplicity={mult}"
            )

            # (Comment translated to English for public release)
            if charge is None or mult is None:
                self.logprint("Warning: charge/multiplicity not found; please enter them in the GUI.")
                if charge is None:
                    charge = simpledialog.askinteger(
                        APP_TITLE = "ORCA -> VEDA .fmu Converter (GUI)"
                        "Please enter the charge.",
                        initialvalue=0,
                        minvalue=-100,
                        maxvalue=100,
                    )
                    if charge is None:
                        raise ValueError("Aborted because charge was not provided.")
                if mult is None:
                    mult = simpledialog.askinteger(
                        APP_TITLE = "ORCA -> VEDA .fmu Converter (GUI)"
                        "Please enter the multiplicity.",
                        initialvalue=1,
                        minvalue=1,
                        maxvalue=21,
                    )
                    if mult is None:
                        raise ValueError("Aborted because multiplicity was not provided.")

            # (Comment translated to English for public release)
            try:
                Z_h, coords_h = parse_atoms_from_hess(hessp)
                Z_list = Z_h
                coords = coords_h
                self.logprint("coords: using $atoms block from .hess (bohr -> Å)")
            except Exception as e:
                self.logprint("coords from .hess failed, keep .out coords:", e)

            # .hess → Hessian
            self.logprint("parse $hessian block from .hess")
            H_full = parse_hess_simple(hessp, natoms=len(Z_list))
            # (Line translated/removed for public release)

            # (Comment translated to English for public release)
            if massw:
                self.logprint("Mass-weighted Hessian -> unweight")
                H_full = unweight_massweighted(H_full, Z_list)

            # (Comment translated to English for public release)
            fac = unit_factor_to_hartree_bohr2(units)
            if abs(fac - 1.0) > 1e-20:
                self.logprint(f"Unit conversion: x {fac:g}  (-> Hartree/Bohr^2)")
                dim = len(H_full)
                for i in range(dim):
                    for j in range(dim):
                        H_full[i][j] *= fac

            

            # (Comment translated to English for public release)
            try:
                out_freqs = parse_out_frequencies(outp)
                if out_freqs:
                    # (Comment translated to English for public release)
                    # (Comment translated to English for public release)
                    f_as_cart = compute_freqs_from_hessian_cm1(H_full, Z_list, assume_massweighted_input=False)
                    f_as_mw   = compute_freqs_from_hessian_cm1(H_full, Z_list, assume_massweighted_input=True)

                    out_sorted = sorted(out_freqs)

                    rms_cart = _rms_diff(out_sorted, f_as_cart, n=20)
                    rms_mw   = _rms_diff(out_sorted, f_as_mw,   n=20)

                    self.logprint(f"freq check (RMS over first modes): as_cart={rms_cart}  as_massweighted={rms_mw}")

                    # (Comment translated to English for public release)
                    if (rms_cart is not None) and (rms_mw is not None):
                        # (Comment translated to English for public release)
                        if rms_mw + 1e-9 < rms_cart * 0.7:
                            self.logprint("Warning: the Hessian appears to be already mass-weighted."
                                          "Turning on the 'mass-weighted -> unweight' option may improve results.")
                        elif rms_cart + 1e-9 < rms_mw * 0.7:
                            self.logprint("Info: the Hessian appears to be Cartesian (not mass-weighted).")

                    # (Comment translated to English for public release)
                    if f_as_cart and out_sorted:
                        # (Comment translated to English for public release)
                        m = min(len(f_as_cart), len(out_sorted))
                        if m >= 10:
                            # (Comment translated to English for public release)
                            a = [abs(x) for x in out_sorted[4:min(14, m)]]
                            b = [abs(x) for x in f_as_cart[4:min(14, m)]]
                            if all(x > 1e-6 for x in b):
                                ratios = [ai/bi for ai,bi in zip(a,b)]
                                med = sorted(ratios)[len(ratios)//2]
                                spread = max(ratios) - min(ratios)
                                if med > 1.5 or med < 0.67:
                                    self.logprint(f"Warning: frequencies seem scaled by about {med:.3g}. The selected Hessian input unit may be incorrect."
                                                  "The selected Hessian input unit may be incorrect.")
            except Exception as e:
                self.logprint("freq check skipped:", e)
# (Comment translated to English for public release)
            write_minimal_fchk(
                fmup,
                Z_list,
                coords,
                _force_int(charge),
                _force_int(mult),
                H_full,
            )
            self.logprint("OK: wrote output ->", fmup)
            # (Comment translated to English for public release)
            if getattr(self, "var_gen_mpo", None) is not None and self.var_gen_mpo.get():
                mpop = fmup.with_suffix('.mpo')
                if mpop.exists():
                    ans = messagebox.askyesno(
                        "Overwrite .mpo?",
                        "A .mpo file with the same name already exists.\n\n"
                        f"Existing file:\n{mpop}\n\n"
                        "Do you want to overwrite it?\n"
                        "Choose 'No' to create a new file with suffix '_gen.mpo' instead."
                    )
                    if not ans:
                        mpop = mpop.with_name(fmup.stem + '_gen.mpo')
                try:
                    write_mpo(mpop, Z_list, coords)
                    self.logprint("OK: generated .mpo ->", mpop)
                except Exception as e:
                    self.logprint("Warning: failed to generate .mpo:", e)

            messagebox.showinfo(APP_TITLE, f"Conversion succeeded.\n{fmup}")
        except Exception as e:
            self.logprint("ERROR:", e)
            self.logprint(traceback.format_exc())
            messagebox.showerror(APP_TITLE, f"Conversion failed:\n{e}")

def _force_int(v):
    try:
        return int(v)
    except Exception:
        return int(float(v))

# ---------- entry point ----------

def main():
    if _has_tkdnd:
        root = TkinterDnD.Tk()
    else:
        root = tk.Tk()
    try:
        root.tk.call("tk", "scaling", 1.25)
    except Exception:
        pass
    App(root)
    root.mainloop()


# Backward-compatible alias
_launch = main

if __name__ == "__main__":
    main()