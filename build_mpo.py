# build_mpo.pyw  (metal-aware defaults)
# ------------------------------------------------------------
# (Comment translated to English for public release)
# (Comment translated to English for public release)
# (Comment translated to English for public release)
# (Comment translated to English for public release)
# (Comment translated to English for public release)
# ------------------------------------------------------------

from __future__ import annotations
import math, os, re, sys
from typing import List, Tuple, Dict, Set
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

# (Comment translated to English for public release)
Z2SYM = {
    1:"H",2:"He",3:"Li",4:"Be",5:"B",6:"C",7:"N",8:"O",9:"F",10:"Ne",
    11:"Na",12:"Mg",13:"Al",14:"Si",15:"P",16:"S",17:"Cl",18:"Ar",
    19:"K",20:"Ca",21:"Sc",22:"Ti",23:"V",24:"Cr",25:"Mn",26:"Fe",27:"Co",28:"Ni",29:"Cu",30:"Zn",
    31:"Ga",32:"Ge",33:"As",34:"Se",35:"Br",36:"Kr",
    37:"Rb",38:"Sr",39:"Y",40:"Zr",41:"Nb",42:"Mo",43:"Tc",44:"Ru",45:"Rh",46:"Pd",47:"Ag",48:"Cd",
    49:"In",50:"Sn",51:"Sb",52:"Te",53:"I",54:"Xe",
    72:"Hf",73:"Ta",74:"W",75:"Re",76:"Os",77:"Ir",78:"Pt",79:"Au",80:"Hg"
}
SYM2Z = {v:k for k,v in Z2SYM.items()}

# (Comment translated to English for public release)
COV_RAD = {
    1:0.31, 5:0.84, 6:0.76, 7:0.71, 8:0.66, 9:0.57,
    14:1.11, 15:1.07, 16:1.05, 17:1.02, 35:1.20, 53:1.39,
    # (Comment translated to English for public release)
    21:1.44,22:1.36,23:1.34,24:1.28,25:1.27,26:1.24,27:1.18,28:1.17,29:1.22,30:1.20,
    39:1.62,40:1.48,41:1.46,42:1.39,44:1.34,45:1.34,46:1.28,47:1.34,48:1.48,
    72:1.50,73:1.38,74:1.30,75:1.28,76:1.26,77:1.27,78:1.36,79:1.36
}
DEFAULT_RAD = 1.20

# (Comment translated to English for public release)
METALS = {
    # (Line translated/removed for public release)
    21,22,23,24,25,26,27,28,29,30, 39,40,41,42,43,44,45,46,47,48,
    57,72,73,74,75,76,77,78,79,80
}

# (Comment translated to English for public release)
def read_xyz(path:str) -> Tuple[List[int], List[Tuple[float,float,float]]]:
    with open(path,'r',encoding='utf-8',errors='ignore') as f:
        lines = f.read().strip().splitlines()
    if len(lines) < 3: raise ValueError("XYZ too short")
    n = int(lines[0].strip())
    Z, R = [], []
    for line in lines[2:2+n]:
        parts = line.split()
        if len(parts) < 4: continue
        sym = parts[0].capitalize()
        znum = SYM2Z.get(sym)
        if not znum: raise ValueError(f"Unknown element symbol in XYZ: {sym}")
        x,y,z = map(float, parts[1:4])
        Z.append(znum); R.append((x,y,z))
    if len(Z) != n: raise ValueError("XYZ atom count mismatch")
    return Z, R

def _read_fchk_block_numbers(lines:List[str], start:int, count:int) -> List[int]:
    out=[]; i=start
    while len(out)<count and i<len(lines):
        for t in lines[i].split():
            out.append(int(t)); 
            if len(out)==count: break
        i+=1
    return out

def _read_fchk_block_reals(lines:List[str], start:int, count:int) -> List[float]:
    out=[]; i=start
    while len(out)<count and i<len(lines):
        for t in lines[i].split():
            out.append(float(t)); 
            if len(out)==count: break
        i+=1
    return out

def read_fmu(path:str) -> Tuple[List[int], List[Tuple[float,float,float]]]:
    with open(path,'r',encoding='ascii',errors='ignore') as f:
        lines = f.read().splitlines()
    pat_nums   = re.compile(r'^Atomic numbers\s+I\s+N=\s*(\d+)\s*$', re.I)
    pat_coords = re.compile(r'^Current cartesian coordinates\s+R\s+N=\s*(\d+)\s*$', re.I)

    nat=None; idx=None
    for i,ln in enumerate(lines):
        m=pat_nums.match(ln)
        if m: nat=int(m.group(1)); idx=i+1; break
    if nat is None: raise ValueError("FMU: 'Atomic numbers' block not found")
    Z = _read_fchk_block_numbers(lines, idx, nat)
    if len(Z)!=nat: raise ValueError("FMU: atomic numbers truncated")

    ncoord=None; idx=None
    for i,ln in enumerate(lines):
        m=pat_coords.match(ln)
        if m: ncoord=int(m.group(1)); idx=i+1; break
    if ncoord is None or ncoord!=3*nat: raise ValueError("FMU: coords block missing/length error")
    flat = _read_fchk_block_reals(lines, idx, 3*nat)
    R=[(flat[3*a],flat[3*a+1],flat[3*a+2]) for a in range(nat)]
    return Z,R

# (Comment translated to English for public release)
def cov_r(z:int)->float: return COV_RAD.get(z, DEFAULT_RAD)
def is_heavy(z:int)->bool: return z!=1
def is_metal(z:int)->bool: return z in METALS
def dist(a,b)->float:
    dx=a[0]-b[0]; dy=a[1]-b[1]; dz=a[2]-b[2]
    return (dx*dx+dy*dy+dz*dz)**0.5

# (Comment translated to English for public release)
def normalized_dist(Zi:int, Zj:int, Ri, Rj)->float:
    rij = cov_r(Zi) + cov_r(Zj)
    return dist(Ri,Rj)/(rij if rij>1e-9 else 1e9)

def cutoff_for_pair(norm:float, Zi:int, Zj:int,
                    fc1_general:float, fc2_general:float,
                    fc1_metal:float,   fc2_metal:float) -> Tuple[bool,bool]:
    """(Docstring translated to English for public release)"""
    if is_metal(Zi) or is_metal(Zj):
        return (norm <= fc1_metal, norm <= fc2_metal)
    else:
        return (norm <= fc1_general, norm <= fc2_general)

# (Comment translated to English for public release)
def build_adjacency(
    Z:List[int], R:List[Tuple[float,float,float]],
    fc1_g:float, fc2_g:float, fc1_m:float, fc2_m:float,
    enable_secondary:bool,
    max_heavy:int|None, max_h:int|None, max_metal:int|None
)->List[Set[int]]:
    n=len(Z)
    adj=[set() for _ in range(n+1)]
    pairs=[]
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            norm = normalized_dist(Z[i-1],Z[j-1],R[i-1],R[j-1])
            pairs.append((norm,i,j))
    pairs.sort()

    # primary
    for norm,i,j in pairs:
        p, s = cutoff_for_pair(norm, Z[i-1],Z[j-1], fc1_g,fc2_g, fc1_m,fc2_m)
        if p:
            adj[i].add(j); adj[j].add(i)

    # (Comment translated to English for public release)
    for i in range(1,n+1):
        if Z[i-1]==1 and not any(is_heavy(Z[k-1]) for k in adj[i]):
            best=None; bestn=1e9
            for j in range(1,n+1):
                if i==j or not is_heavy(Z[j-1]): continue
                norm = normalized_dist(Z[i-1],Z[j-1],R[i-1],R[j-1])
                _, s = cutoff_for_pair(norm, Z[i-1],Z[j-1], fc1_g,fc2_g, fc1_m,fc2_m)
                if s and norm<bestn:
                    bestn=norm; best=j
            if best:
                adj[i].add(best); adj[best].add(i)

    # (Comment translated to English for public release)
    if enable_secondary:
        for norm,i,j in pairs:
            if j in adj[i]: continue
            _, s = cutoff_for_pair(norm, Z[i-1],Z[j-1], fc1_g,fc2_g, fc1_m,fc2_m)
            if s:
                adj[i].add(j); adj[j].add(i)

    # (Comment translated to English for public release)
    comp=[-1]*(n+1); cid=0
    for v in range(1,n+1):
        if comp[v]!=-1: continue
        stack=[v]; comp[v]=cid
        while stack:
            x=stack.pop()
            for y in adj[x]:
                if comp[y]==-1: comp[y]=cid; stack.append(y)
        cid+=1
    if cid>1:
        for _,i,j in pairs:
            if comp[i]!=comp[j]:
                adj[i].add(j); adj[j].add(i)
                old=comp[j]; new=comp[i]
                for v in range(1,n+1):
                    if comp[v]==old: comp[v]=new
                if len(set(comp[1:]))==1: break

    # (Comment translated to English for public release)
    for i in range(1,n+1):
        z=Z[i-1]
        if z==1: lim=max_h
        elif is_metal(z): lim=max_metal
        else: lim=max_heavy
        if lim is None or len(adj[i])<=lim: continue
        # (Comment translated to English for public release)
        neigh=sorted(list(adj[i]), key=lambda j: dist(R[i-1],R[j-1]))
        keep=set(neigh[:lim]); drop=set(neigh[lim:])
        for j in drop:
            adj[i].discard(j); adj[j].discard(i)

    return adj[1:]

# ==== MPO I/O ====
def write_mpo(path:str, Z:List[int], adj:List[Set[int]]):
    with open(path,'w',encoding='ascii',newline='\r\n') as f:
        n=len(adj)
        for i in range(1,n+1):
            neigh=sorted(adj[i-1])
            sym = Z2SYM.get(Z[i-1], f"Z{Z[i-1]}")
            line = f"{len(neigh)} . {sym} {i} :"
            if neigh:
                items=[]
                for j in neigh:
                    sj = Z2SYM.get(Z[j-1], f"Z{Z[j-1]}")
                    items.append(f"{sj} {j}")
                line += " " + " , ".join(items)
            f.write(line+"\r\n")

def read_mpo(path:str)->List[List[int]]:
    with open(path,'r',encoding='utf-8',errors='ignore') as f:
        lines=[ln.strip() for ln in f if ln.strip()]
    neighbors:Dict[int,List[int]]={}
    pat=re.compile(r'^\s*(\d+)\s*\.\s*([A-Za-z]{1,2})\s+(\d+)\s*:\s*(.*)$')
    for ln in lines:
        m=pat.match(ln)
        if not m:
            m2=re.match(r'^\s*0\s*\.\s*([A-Za-z]{1,2})\s+(\d+)\s*:\s*$', ln)
            if m2: neighbors[int(m2.group(2))]=[]; continue
            continue
        idx=int(m.group(3)); rest=m.group(4).strip()
        lst=[]
        if rest:
            for t in rest.split(','):
                p=t.strip().split()
                if p and p[-1].isdigit(): lst.append(int(p[-1]))
        neighbors[idx]=lst
    if not neighbors: return []
    n=max(neighbors.keys())
    out=[[] for _ in range(n)]
    for i,lst in neighbors.items(): out[i-1]=sorted(set(lst))
    return out

def repair_mpo(
    Z:List[int], R:List[Tuple[float,float,float]], mpo_in:str,
    fc1_g:float, fc2_g:float, fc1_m:float, fc2_m:float,
    enable_secondary:bool, max_heavy:int|None, max_h:int|None, max_metal:int|None
)->List[Set[int]]:
    base = build_adjacency(Z,R, fc1_g,fc2_g, fc1_m,fc2_m, enable_secondary, max_heavy,max_h,max_metal)
    exist = read_mpo(mpo_in)
    if not exist: return base
    n=len(Z); adj=[set(s) for s in base]
    # (Comment translated to English for public release)
    for i in range(1,n+1):
        if i-1>=len(exist): continue
        for j in exist[i-1]:
            if not (1<=j<=n) or j==i: continue
            norm = normalized_dist(Z[i-1],Z[j-1],R[i-1],R[j-1])
            p,s = cutoff_for_pair(norm, Z[i-1],Z[j-1], fc1_g,fc2_g, fc1_m,fc2_m)
            if p or (enable_secondary and s):
                adj[i-1].add(j); adj[j-1].add(i)
    # (Comment translated to English for public release)
    for i in range(1,n+1):
        z=Z[i-1]
        if z==1: lim=max_h
        elif is_metal(z): lim=max_metal
        else: lim=max_heavy
        if lim is None or len(adj[i-1])<=lim: continue
        neigh=sorted(list(adj[i-1]), key=lambda j: dist(R[i-1],R[j-1]))
        keep=set(neigh[:lim]); drop=set(neigh[lim:])
        for j in drop:
            adj[i-1].discard(j); adj[j-1].discard(i)
    return adj

# ==== GUI ====
APP_TITLE = "ORCA -> VEDA .fmu Converter (GUI)"

# (Comment translated to English for public release)
DEF_PRIMARY_GENERAL  = 1.25
DEF_SECONDARY_GENERAL= 1.45
DEF_PRIMARY_METAL    = 1.40
DEF_SECONDARY_METAL  = 1.85
DEF_MAX_HEAVY        = 4
DEF_MAX_H            = 2
DEF_MAX_METAL        = 6
DEF_ENABLE_SECONDARY = True

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE); self.geometry("780x560"); self.minsize(760,560)

        self.var_src=tk.StringVar(); self.var_mpo_in=tk.StringVar(); self.var_out=tk.StringVar()
        self.var_fc1g=tk.DoubleVar(value=DEF_PRIMARY_GENERAL)
        self.var_fc2g=tk.DoubleVar(value=DEF_SECONDARY_GENERAL)
        self.var_fc1m=tk.DoubleVar(value=DEF_PRIMARY_METAL)
        self.var_fc2m=tk.DoubleVar(value=DEF_SECONDARY_METAL)
        self.var_sec=tk.BooleanVar(value=DEF_ENABLE_SECONDARY)
        self.var_maxHvy=tk.StringVar(value=str(DEF_MAX_HEAVY))
        self.var_maxH  =tk.StringVar(value=str(DEF_MAX_H))
        self.var_maxMet=tk.StringVar(value=str(DEF_MAX_METAL))

        frm=ttk.Frame(self,padding=10); frm.pack(fill="both",expand=True)
        # (Comment translated to English for public release)
        ttk.Label(frm,text="Source (.fmu/.fchk/.fch/.fmt or .xyz):").grid(row=0,column=0,sticky="w")
        ttk.Entry(frm,textvariable=self.var_src).grid(row=0,column=1,sticky="ew",padx=6)
        ttk.Button(frm,text="Browse...",command=self.browse_src).grid(row=0,column=2,sticky="ew")
        ttk.Label(frm,text="Repair from .mpo (optional):").grid(row=1,column=0,sticky="w")
        ttk.Entry(frm,textvariable=self.var_mpo_in).grid(row=1,column=1,sticky="ew",padx=6)
        ttk.Button(frm,text="Browse...",command=self.browse_mpo).grid(row=1,column=2,sticky="ew")
        ttk.Label(frm,text="Output .mpo:").grid(row=2,column=0,sticky="w")
        ttk.Entry(frm,textvariable=self.var_out).grid(row=2,column=1,sticky="ew",padx=6)
        ttk.Button(frm,text="Save as...",command=self.browse_out).grid(row=2,column=2,sticky="ew")
        frm.grid_columnconfigure(1,weight=1)

        ttk.Separator(frm).grid(row=3,column=0,columnspan=3,sticky="ew",pady=8)
        # (Comment translated to English for public release)
        box=ttk.LabelFrame(frm,text="Cutoffs (normalized distance = d / (r_covA + r_covB))")
        box.grid(row=4,column=0,columnspan=3,sticky="ew")
        g=ttk.Frame(box); g.pack(fill="x",padx=6,pady=4)
        ttk.Label(g,text="General primary").grid(row=0,column=0,sticky="w")
        ttk.Entry(g,width=8,textvariable=self.var_fc1g).grid(row=0,column=1,sticky="w",padx=6)
        ttk.Label(g,text="General secondary").grid(row=0,column=2,sticky="w")
        ttk.Entry(g,width=8,textvariable=self.var_fc2g).grid(row=0,column=3,sticky="w",padx=6)
        ttk.Label(g,text="Metal primary").grid(row=1,column=0,sticky="w")
        ttk.Entry(g,width=8,textvariable=self.var_fc1m).grid(row=1,column=1,sticky="w",padx=6)
        ttk.Label(g,text="Metal secondary").grid(row=1,column=2,sticky="w")
        ttk.Entry(g,width=8,textvariable=self.var_fc2m).grid(row=1,column=3,sticky="w",padx=6)
        ttk.Checkbutton(g,text="Enable secondary near-contacts",variable=self.var_sec).grid(row=2,column=0,columnspan=4,sticky="w",pady=4)

        # (Comment translated to English for public release)
        box2=ttk.LabelFrame(frm,text="Max neighbors (pruning)")
        box2.grid(row=5,column=0,columnspan=3,sticky="ew",pady=6)
        h=ttk.Frame(box2); h.pack(fill="x",padx=6,pady=4)
        ttk.Label(h,text="Heavy (non-H, non-metal)").grid(row=0,column=0,sticky="w")
        ttk.Entry(h,width=6,textvariable=self.var_maxHvy).grid(row=0,column=1,sticky="w",padx=6)
        ttk.Label(h,text="H").grid(row=0,column=2,sticky="w",padx=(12,0))
        ttk.Entry(h,width=6,textvariable=self.var_maxH).grid(row=0,column=3,sticky="w",padx=6)
        ttk.Label(h,text="Metal").grid(row=0,column=4,sticky="w",padx=(12,0))
        ttk.Entry(h,width=6,textvariable=self.var_maxMet).grid(row=0,column=5,sticky="w",padx=6)

        # (Comment translated to English for public release)
        btns=ttk.Frame(frm); btns.grid(row=6,column=0,columnspan=3,sticky="ew",pady=8)
        ttk.Button(btns,text="Build/Repair MPO",command=self.run).pack(side="left")
        ttk.Button(btns,text="Close",command=self.destroy).pack(side="right")

        # (Comment translated to English for public release)
        ttk.Label(frm,text="Log:").grid(row=7,column=0,sticky="w")
        self.txt=tk.Text(frm,height=14); self.txt.grid(row=8,column=0,columnspan=3,sticky="nsew")
        frm.grid_rowconfigure(8,weight=1)

        self.after(50,self.init_from_argv)

    def log(self,s:str):
        self.txt.insert("end",s+"\n"); self.txt.see("end"); self.update_idletasks()

    def browse_src(self):
        p=filedialog.askopenfilename(title="Select .fmu/.fchk/.fch/.fmt or .xyz",
            filetypes=[("FMU/FCHK/FMT","*.fmu;*.fchk;*.fch;*.fmt"),("XYZ","*.xyz"),("All","*.*")])
        if p: self.var_src.set(p); self.suggest_out()

    def browse_mpo(self):
        p=filedialog.askopenfilename(title="Select .mpo", filetypes=[("MPO","*.mpo"),("All","*.*")])
        if p: self.var_mpo_in.set(p)

    def browse_out(self):
        p=filedialog.asksaveasfilename(title="Save .mpo",defaultextension=".mpo",
            filetypes=[("MPO","*.mpo"),("All","*.*")])
        if p: self.var_out.set(p)

    def suggest_out(self):
        src=self.var_src.get().strip()
        if not src: return
        root,_=os.path.splitext(src)
        self.var_out.set(root+".mpo")

    def init_from_argv(self):
        files=[a for a in sys.argv[1:] if os.path.isfile(a)]
        if not files: return
        fmu=next((p for p in files if os.path.splitext(p)[1].lower() in (".fmu",".fchk",".fch",".fmt")),None)
        xyz=next((p for p in files if os.path.splitext(p)[1].lower()==".xyz"),None)
        mpo=next((p for p in files if os.path.splitext(p)[1].lower()==".mpo"),None)
        if fmu or xyz: self.var_src.set(fmu or xyz); self.suggest_out()
        if mpo: self.var_mpo_in.set(mpo)

    def read_src(self)->Tuple[List[int],List[Tuple[float,float,float]]]:
        src=self.var_src.get().strip()
        if not src: raise ValueError("Source not set")
        ext=os.path.splitext(src)[1].lower()
        if ext in (".fmu",".fchk",".fch",".fmt"): 
            Z,R=read_fmu(src); self.log(f"read FMU: natoms={len(Z)}"); return Z,R
        if ext==".xyz":
            Z,R=read_xyz(src); self.log(f"read XYZ: natoms={len(Z)}"); return Z,R
        raise ValueError("Source must be .fmu/.fchk/.fch/.fmt or .xyz")

    def run(self):
        try:
            Z,R=self.read_src()
            outp=self.var_out.get().strip()
            if not outp: self.suggest_out(); outp=self.var_out.get().strip()
            if not outp: raise ValueError("Output .mpo not set")

            # (Comment translated to English for public release)
            fc1g=float(self.var_fc1g.get()); fc2g=float(self.var_fc2g.get())
            fc1m=float(self.var_fc1m.get()); fc2m=float(self.var_fc2m.get())
            sec =bool(self.var_sec.get())
            maxHvy=int(self.var_maxHvy.get()) if self.var_maxHvy.get().strip() else None
            maxH  =int(self.var_maxH.get())   if self.var_maxH.get().strip()   else None
            maxMet=int(self.var_maxMet.get()) if self.var_maxMet.get().strip() else None

            mpoin=self.var_mpo_in.get().strip()
            if mpoin:
                self.log(f"repair from: {mpoin}")
                adj=repair_mpo(Z,R, mpoin, fc1g,fc2g, fc1m,fc2m, sec, maxHvy,maxH,maxMet)
            else:
                adj=build_adjacency(Z,R, fc1g,fc2g, fc1m,fc2m, sec, maxHvy,maxH,maxMet)

            write_mpo(outp, Z, adj)
            self.log(f"OK: wrote MPO -> {outp}")
            messagebox.showinfo(APP_TITLE, f"OK: wrote MPO ->\n{outp}")
        except Exception as e:
            self.log(f"ERROR: {type(e).__name__}: {e}")
            messagebox.showerror(APP_TITLE, f"ERROR: {e}")

def main():
    app=App(); app.mainloop()

if __name__=="__main__":
    main()
