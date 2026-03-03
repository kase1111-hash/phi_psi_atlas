#!/usr/bin/env python3
"""
Extended Gradient Test: 5 Additional NMR-Characterized Proteins
================================================================
Doubles the dataset from 4 to 9 proteins to pin down the crossover threshold.

New proteins chosen for:
  - Spread across the Gini spectrum  
  - Well-characterized NMR dynamics (S² data published)
  - Clean apo / bound PDB pairs (same chain, same organism)
  - Diverse fold classes and perturbation types

New set:
  1. Barnase (110 res, α+β ribonuclease)
     Apo: 1A2P (X-ray, 1.5 Å)
     Bound: 1BRS chain A (barnase-barstar complex, 2.0 Å)
     NMR: S² ~ 0.85 overall (Sahu et al. 2000)
     
  2. BPTI — bovine pancreatic trypsin inhibitor (58 res, mixed α/β)
     Apo: 5PTI (X-ray, 1.0 Å) — one of the most studied proteins
     Bound: 3BTK chain I (BPTI-trypsin complex, 1.85 Å)
     NMR: S² ~ 0.85 average (Bae et al. 1998, Otting et al. 1993)
     
  3. RNase A — bovine ribonuclease A (124 res, α+β)
     Apo: 7RSA (X-ray, 1.26 Å)
     Bound: 1RPH (RNase A with phosphate bound, 1.5 Å)
     Also: 1AFU (RNase A + uridine vanadate inhibitor)
     NMR: S² ~ 0.87 average (Cole & Loria 2002)
     
  4. FKBP12 — FK506 binding protein (107 res, β-rich)
     Apo: 2PPN (X-ray, 1.7 Å, unliganded FKBP12)
     Bound: 1FKJ (FKBP12 + FK506, 1.7 Å) 
     NMR: S² ~ 0.86 average (Cheng et al. 1994)
     
  5. Thioredoxin — E. coli (108 res, α/β/α sandwich)
     Apo: 2TRX (X-ray, 1.68 Å, oxidized)
     Bound: 1THO (thioredoxin-peptide complex) — or use reduced form
     Alt: 1XOB (reduced thioredoxin, 1.2 Å) — redox as "perturbation"
     NMR: S² ~ 0.88 average (Stone et al. 1993)

Usage:
    python gradient_test_extended.py
    
    Requires: numpy. Downloads PDB files automatically.

Author: Kase / True North Construction LLC  
License: CC0 1.0 Universal
"""

import numpy as np
import os, sys, json
from urllib.request import urlretrieve
from pathlib import Path
from collections import OrderedDict

DATA_DIR = Path("gradient_test_extended")


# ═══════════════════════════════════════════════════════════════════
# STRUCTURE DEFINITIONS
# ═══════════════════════════════════════════════════════════════════

PROTEINS = OrderedDict([

    # ── From the first round (for combined analysis) ──
    ("calmodulin", {
        "apo": {"pdb": "1CFD", "chain": "A", "is_nmr": True,
                "label": "Apo CaM (NMR)"},
        "bound": {"pdb": "1CLL", "chain": "A", "is_nmr": False,
                  "label": "Holo CaM (X-ray, Ca²⁺)",
                  "perturbation": "Ca²⁺ binding"},
    }),
    ("ubiquitin", {
        "apo": {"pdb": "1UBQ", "chain": "A", "is_nmr": False,
                "label": "Apo Ubiquitin (1.8 Å)"},
        "bound": {"pdb": "1S1Q", "chain": "B", "is_nmr": False,
                  "label": "Ub-TSG101 Complex",
                  "perturbation": "TSG101 UEV binding"},
    }),
    ("gb1", {
        "apo": {"pdb": "1PGA", "chain": "A", "is_nmr": False,
                "label": "Apo GB1 (2.07 Å)"},
        "bound": {"pdb": "1PGB", "chain": "A", "is_nmr": False,
                  "label": "GB1 alt crystal form",
                  "perturbation": "Crystal repacking"},
    }),
    ("lysozyme", {
        "apo": {"pdb": "1AKI", "chain": "A", "is_nmr": False,
                "label": "Apo HEWL (1.5 Å)"},
        "bound": {"pdb": "1HEW", "chain": "A", "is_nmr": False,
                  "label": "HEWL + tri-NAG",
                  "perturbation": "tri-NAG inhibitor"},
    }),

    # ── NEW: 5 additional proteins ──

    ("barnase", {
        "apo": {"pdb": "1A2P", "chain": "A", "is_nmr": False,
                "label": "Apo Barnase (1.5 Å)"},
        "bound": {"pdb": "1BRS", "chain": "A", "is_nmr": False,
                  "label": "Barnase-Barstar Complex (2.0 Å)",
                  "perturbation": "Barstar binding (Kd ~ 10⁻¹⁴ M)"},
        "extra_bound": {"pdb": "1BRN", "chain": "A", "is_nmr": False,
                  "label": "Barnase + d(CGAC) inhibitor",
                  "perturbation": "DNA tetranucleotide inhibitor"},
        "notes": "α+β ribonuclease, 110 res. Ultra-tight barstar complex.",
    }),

    ("bpti", {
        "apo": {"pdb": "5PTI", "chain": "A", "is_nmr": False,
                "label": "Apo BPTI (1.0 Å)"},
        "bound": {"pdb": "3BTK", "chain": "I", "is_nmr": False,
                  "label": "BPTI-Trypsin Complex (1.85 Å)",
                  "perturbation": "Trypsin binding (Kunitz-type inhibition)"},
        "notes": "58 res, mixed α/β, one of most studied small proteins.",
    }),

    ("rnase_a", {
        "apo": {"pdb": "7RSA", "chain": "A", "is_nmr": False,
                "label": "Apo RNase A (1.26 Å)"},
        "bound": {"pdb": "1RPH", "chain": "A", "is_nmr": False,
                  "label": "RNase A + phosphate (1.5 Å)",
                  "perturbation": "Phosphate ion in active site"},
        "extra_bound": {"pdb": "1AFU", "chain": "A", "is_nmr": False,
                  "label": "RNase A + uridine vanadate",
                  "perturbation": "Uridine vanadate transition-state analog"},
        "notes": "124 res, α+β. Classic enzyme. Phosphate is minimal perturbation.",
    }),

    ("fkbp12", {
        "apo": {"pdb": "2PPN", "chain": "A", "is_nmr": False,
                "label": "Apo FKBP12 (1.7 Å)"},
        "bound": {"pdb": "1FKJ", "chain": "A", "is_nmr": False,
                  "label": "FKBP12 + FK506 (1.7 Å)",
                  "perturbation": "FK506 immunosuppressant binding"},
        "notes": "107 res, β-rich (β-sheet wrap). Drug target for immunosuppression.",
    }),

    ("thioredoxin", {
        "apo": {"pdb": "2TRX", "chain": "A", "is_nmr": False,
                "label": "Apo Thioredoxin oxidized (1.68 Å)"},
        "bound": {"pdb": "1XOB", "chain": "A", "is_nmr": False,
                  "label": "Thioredoxin reduced (1.2 Å)",
                  "perturbation": "Disulfide reduction (redox switch)"},
        "notes": "108 res, α/β/α sandwich. Redox perturbation rather than binding.",
    }),
])


# ═══════════════════════════════════════════════════════════════════
# CORE COMPUTATION (identical to previous scripts)
# ═══════════════════════════════════════════════════════════════════

def download_pdb(pdb_id, data_dir=DATA_DIR):
    data_dir.mkdir(exist_ok=True)
    fp = data_dir / f"{pdb_id.upper()}.pdb"
    if fp.exists(): return fp
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    print(f"    Downloading {pdb_id} ... ", end="", flush=True)
    try:
        urlretrieve(url, fp); print("ok")
    except Exception as e:
        print(f"FAIL ({e})"); return None
    return fp

def parse_backbone(filepath, chain_id="A", model_num=1):
    residues = {}; cm = 0; inm = False
    with open(filepath) as f:
        for line in f:
            if line.startswith("MODEL"):
                cm = int(line[10:14].strip()); inm = (cm == model_num); continue
            if line.startswith("ENDMDL"):
                if inm: break
                inm = False; continue
            if cm == 0: inm = True
            if not line.startswith("ATOM") or not inm: continue
            an = line[12:16].strip()
            if an not in ("N","CA","C"): continue
            ch = line[21].strip()
            if ch != chain_id: continue
            try: rn = int(line[22:26].strip())
            except: continue
            if line[26].strip(): continue
            x,y,z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if rn not in residues: residues[rn] = {}
            residues[rn][an] = np.array([x,y,z])
    return residues

def dihedral(p1,p2,p3,p4):
    b1=p2-p1; b2=p3-p2; b3=p4-p3
    n1=np.cross(b1,b2); n2=np.cross(b2,b3)
    nn1=np.linalg.norm(n1); nn2=np.linalg.norm(n2)
    if nn1<1e-10 or nn2<1e-10: return np.nan
    n1/=nn1; n2/=nn2
    m1=np.cross(n1, b2/np.linalg.norm(b2))
    return np.arctan2(np.dot(m1,n2), np.dot(n1,n2))

def get_dihedrals(residues):
    sr = sorted(residues.keys()); out = {}
    for i,rn in enumerate(sr):
        r = residues[rn]
        if not all(k in r for k in ("N","CA","C")): continue
        phi = np.nan
        if i>0:
            prn = sr[i-1]; pr = residues.get(prn,{})
            if "C" in pr and (rn-prn)==1:
                phi = dihedral(pr["C"],r["N"],r["CA"],r["C"])
        psi = np.nan
        if i<len(sr)-1:
            nrn = sr[i+1]; nr = residues.get(nrn,{})
            if "N" in nr and (nrn-rn)==1:
                psi = dihedral(r["N"],r["CA"],r["C"],nr["N"])
        out[rn] = (phi,psi)
    return out

def adiff(a,b): return np.arctan2(np.sin(a-b), np.cos(a-b))

def get_kappa_sq(dih):
    sr = sorted(dih.keys()); ksq = {}
    for i in range(1, len(sr)-1):
        rp,rc,rn = sr[i-1],sr[i],sr[i+1]
        if (rc-rp)!=1 or (rn-rc)!=1: continue
        pp,sp = dih[rp]; pc,sc = dih[rc]; pn,sn = dih[rn]
        if any(np.isnan(x) for x in [pp,sp,pc,sc,pn,sn]): continue
        dp1=adiff(pc,pp); ds1=adiff(sc,sp)
        dp2=adiff(pn,pc); ds2=adiff(sn,sc)
        a1=np.sqrt(dp1**2+ds1**2); a2=np.sqrt(dp2**2+ds2**2)
        if a1<1e-10 or a2<1e-10: continue
        dTp=dp2/a2-dp1/a1; dTs=ds2/a2-ds1/a1
        dm = 0.5*(a1+a2)
        ksq[rc] = (dTp**2+dTs**2)/(dm**2)
    return ksq

def gini(values):
    v = np.sort(np.array(values,dtype=float))
    v = v[~np.isnan(v)]
    if len(v)==0: return np.nan
    n=len(v); s=np.sum(v)
    if s<1e-15: return 0.0
    return (2*np.sum(np.arange(1,n+1)*v)/(n*s)) - (n+1)/n

def rama_class(phi,psi):
    if np.isnan(phi) or np.isnan(psi): return "unk"
    pd=np.degrees(phi); sd=np.degrees(psi)
    if -180<=pd<=-30 and -80<=sd<=-10: return "alpha"
    if -180<=pd<=-30 and (80<=sd<=180 or -10<=sd<=80): return "beta"
    if 30<=pd<=120 and -60<=sd<=60: return "alphaL"
    return "coil"


# ═══════════════════════════════════════════════════════════════════
# PIPELINE
# ═══════════════════════════════════════════════════════════════════

def analyze(filepath, chain, label, is_nmr=False):
    models = [1]
    if is_nmr:
        mc=0
        with open(filepath) as f:
            for line in f:
                if line.startswith("MODEL"): mc+=1
        if mc>1: models = list(range(1,min(mc+1,21)))
    
    all_g = []; primary = None
    for m in models:
        res = parse_backbone(filepath, chain, m)
        if len(res)<5: continue
        d = get_dihedrals(res); k = get_kappa_sq(d)
        if len(k)<3: continue
        g = gini(list(k.values())); all_g.append(g)
        if m==1:
            ss = {"alpha":0,"beta":0,"coil":0}
            for rn in k:
                if rn in d:
                    c = rama_class(*d[rn])
                    if c in ss: ss[c]+=1
            primary = {"gini":g, "n_res":len(res), "n_kappa":len(k),
                       "kappa_sq":k, "ss":ss}
    
    if primary is None:
        print(f"    [ERROR] {label}"); return None
    
    vals = list(primary["kappa_sq"].values())
    total_ss = sum(primary["ss"].values())
    ss_str = ""
    if total_ss > 0:
        ss_str = (f"  α={100*primary['ss']['alpha']/total_ss:.0f}% "
                  f"β={100*primary['ss']['beta']/total_ss:.0f}%")
    
    ens_str = ""
    if len(all_g)>1:
        primary["ens_mean"] = np.mean(all_g)
        primary["ens_sd"] = np.std(all_g)
        ens_str = f"  (ens: {np.mean(all_g):.4f}±{np.std(all_g):.4f})"
    
    print(f"    {label}")
    print(f"      n={primary['n_res']:>3d}  κ²={primary['n_kappa']:>3d}  "
          f"Σκ²={np.sum(vals):.0f}  Gini={primary['gini']:.4f}{ens_str}{ss_str}")
    
    return primary


def main():
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  EXTENDED GRADIENT TEST: 9 Proteins (4 original + 5 new)       ║")
    print("║  Goal: Pin down the crossover threshold                        ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    
    # Download
    print("\n── Downloads ──")
    paths = {}
    for pn,pd in PROTEINS.items():
        for sk in ["apo","bound","extra_bound"]:
            if sk not in pd: continue
            pid = pd[sk]["pdb"]
            if pid not in paths:
                p = download_pdb(pid)
                if p: paths[pid] = p
    
    # Analyze
    print("\n── Analysis ──")
    results = {}
    for pn,pd in PROTEINS.items():
        print(f"\n  {pn.upper()}" + (f"  ({pd.get('notes','')})" if 'notes' in pd else ""))
        results[pn] = {}
        for sk in ["apo","bound","extra_bound"]:
            if sk not in pd: continue
            st = pd[sk]; pid = st["pdb"]
            if pid not in paths: 
                print(f"    [SKIP] {pid}"); continue
            r = analyze(paths[pid], st["chain"], f"{pid}: {st['label']}", st.get("is_nmr",False))
            if r:
                results[pn][sk] = r
                results[pn][sk]["pdb"] = pid
    
    # ═══════════════════════════════════════════════════════════════
    # GRADIENT TEST
    # ═══════════════════════════════════════════════════════════════
    
    tests = []
    for pn,pd in PROTEINS.items():
        if "apo" not in results[pn]: continue
        g_apo = results[pn]["apo"]["gini"]
        for sk in ["bound","extra_bound"]:
            if sk not in results[pn]: continue
            g_bnd = results[pn][sk]["gini"]
            dg = g_bnd - g_apo
            pert = pd[sk].get("perturbation","?")
            tests.append({
                "protein": pn,
                "pdb_apo": pd["apo"]["pdb"],
                "pdb_bnd": pd[sk]["pdb"],
                "perturbation": pert,
                "g_apo": g_apo,
                "g_bound": g_bnd,
                "dg": dg,
            })
    
    # Threshold scan
    print("\n\n" + "="*70)
    print("THRESHOLD SCAN")
    print("="*70)
    
    best_score = 0; best_range = [0,0]
    scores_by_thresh = {}
    
    for ti in range(700, 950):
        th = ti/1000
        score = 0
        for t in tests:
            pred = "+" if t["g_apo"] < th else "-"
            act = "+" if t["dg"] > 0 else "-"
            if pred == act: score += 1
        scores_by_thresh[th] = score
        if score > best_score:
            best_score = score
            best_range = [th, th]
        elif score == best_score and th == best_range[1] + 0.001:
            best_range[1] = th
    
    # Find perfect range
    perfect_start = None; perfect_end = None
    for ti in range(700, 950):
        th = ti/1000
        if scores_by_thresh[th] == best_score:
            if perfect_start is None: perfect_start = th
            perfect_end = th
    
    print(f"\n  Best score: {best_score}/{len(tests)}")
    print(f"  Perfect range: [{perfect_start:.3f}, {perfect_end:.3f}]")
    if perfect_start and perfect_end:
        print(f"  Midpoint: {(perfect_start+perfect_end)/2:.3f}")
        print(f"  Width: {perfect_end-perfect_start:.3f}")
    
    # Show results at key thresholds
    for label, th in [("Original (0.810)", 0.810), 
                       ("AlphaFold mean (0.835)", 0.835),
                       ("Midpoint of perfect range", (perfect_start+perfect_end)/2 if perfect_start else 0.86)]:
        score = scores_by_thresh.get(th, scores_by_thresh.get(round(th,3), 0))
        # Recount
        sc = 0; misses = []
        for t in tests:
            pred = "+" if t["g_apo"] < th else "-"
            act = "+" if t["dg"] > 0 else "-"
            if pred == act: sc += 1
            else: misses.append(f"{t['protein']}({t['pdb_bnd']})")
        miss_str = ", ".join(misses) if misses else "none"
        print(f"  At {label}: {sc}/{len(tests)}  misses: {miss_str}")
    
    # Full results table
    print("\n\n" + "="*70)
    print("FULL RESULTS TABLE")
    print("="*70)
    print(f"\n  {'Protein':12s} {'Apo':>5s} {'G_apo':>7s} │ {'Bnd':>5s} {'G_bnd':>7s} {'ΔGini':>8s} {'Dir':>5s}")
    print("  " + "─"*60)
    
    for t in sorted(tests, key=lambda x: x["g_apo"]):
        d = "STIFF" if t["dg"]>0 else "LOOSE"
        print(f"  {t['protein']:12s} {t['pdb_apo']:>5s} {t['g_apo']:7.4f} │ "
              f"{t['pdb_bnd']:>5s} {t['g_bound']:7.4f} {t['dg']:+8.4f} {d:>5s}")
    
    # Sorted by baseline — the gradient visualization
    print("\n\n" + "="*70)
    print("GRADIENT VISUALIZATION (sorted by baseline Gini)")
    print("="*70)
    
    baselines = {}
    for t in tests:
        if t["protein"] not in baselines:
            baselines[t["protein"]] = {"g": t["g_apo"], "directions": []}
        baselines[t["protein"]]["directions"].append("↑" if t["dg"]>0 else "↓")
    
    print(f"\n  {'Protein':12s} {'Baseline':>8s}  {'Dir':>4s}  Bar")
    print("  " + "─"*55)
    
    for pn, info in sorted(baselines.items(), key=lambda x: x[1]["g"]):
        dirs = "".join(info["directions"])
        g = info["g"]
        # Visual bar
        pos = int((g - 0.65) / 0.30 * 40)
        bar = "·" * pos + "█" + "·" * (40 - pos)
        print(f"  {pn:12s} {g:8.4f}  {dirs:>4s}  |{bar}|")
    
    # Draw threshold
    if perfect_start and perfect_end:
        mid = (perfect_start + perfect_end) / 2
        pos = int((mid - 0.65) / 0.30 * 40)
        thresh_bar = " " * pos + "▲"
        print(f"  {'THRESHOLD':12s} {mid:8.3f}        |{thresh_bar}")
    
    # Export
    export = {"tests": tests, "best_score": best_score, 
              "n_tests": len(tests),
              "perfect_range": [perfect_start, perfect_end],
              "midpoint": (perfect_start+perfect_end)/2 if perfect_start else None}
    
    # Add per-residue data
    export["proteins"] = {}
    for pn in results:
        export["proteins"][pn] = {}
        for sk in results[pn]:
            r = results[pn][sk]
            export["proteins"][pn][sk] = {
                "pdb": r.get("pdb","?"),
                "gini": r["gini"],
                "n_res": r["n_res"],
                "n_kappa": r["n_kappa"],
                "kappa_sq": {str(k):v for k,v in r["kappa_sq"].items()},
            }
    
    out_path = DATA_DIR / "extended_gradient_results.json"
    with open(out_path, "w") as f:
        json.dump(export, f, indent=2)
    print(f"\nExported to: {out_path}")


if __name__ == "__main__":
    main()
