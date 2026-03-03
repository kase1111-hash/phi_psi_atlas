#!/usr/bin/env python3
"""
Large-Scale Gradient Test v3 — Bulletproof Edition
====================================================
Uses ONLY:
  - RCSB Search API (find structures)
  - RCSB REST entry API (get resolution + UniProt mapping)
  - Direct PDB file download + parsing (get chains, ligands, curvature)

No GraphQL. No polymer_entities endpoint. No nonpolymer_entities endpoint.

Usage:
    python gradient_test_v3.py [--max-pairs 200] [--resolution 2.0] [--mode apo-holo|temperature|both]
"""

import numpy as np
import os, sys, json, time, gzip, re
from urllib.request import urlopen, Request, urlretrieve
from urllib.error import HTTPError
from pathlib import Path
from collections import defaultdict
import argparse

DATA_DIR = Path("gradient_v3")
PDB_DIR = DATA_DIR / "pdbs"

SOLVENT_LIGANDS = {
    "HOH","DOD","WAT","SO4","PO4","ACT","GOL","EDO","PEG","DMS","MPD",
    "CL","NA","K","MG","CA","ZN","MN","FE","CU","NI","CO",
    "BME","DTT","TRS","CIT","TAR","FMT","IMD",
    "EPE","MES","HED","PGE","1PE","P6G","12P",
    "IOD","BR","F","SCN","NO3",
    "CD","HG","SR","BA","CS","PB","PT",
    "NH4","LI","RB","AZI","BOG","C8E","LDA","SDS","UNX",
    "OCS","CME","CSO","KCX","LLP","TPO","SEP","PTR","MSE",
}

# ═══════════════════════════════════════════════════════════════════
# RCSB SEARCH API
# ═══════════════════════════════════════════════════════════════════

def rcsb_search(query_json):
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    req = Request(url, data=json.dumps(query_json).encode(),
                  headers={"Content-Type": "application/json"})
    try:
        with urlopen(req, timeout=60) as resp:
            return json.loads(resp.read())
    except HTTPError as e:
        body = ""
        try: body = e.read().decode('utf-8', errors='replace')[:500]
        except: pass
        print(f"  Search API HTTP {e.code}: {body[:200]}")
        return None
    except Exception as e:
        print(f"  Search API error: {e}")
        return None

def search_xray(resolution=2.0, max_results=3000):
    """Search for X-ray structures."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "exptl.method",
                                "operator": "exact_match",
                                "value": "X-RAY DIFFRACTION"}},
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                                "operator": "range",
                                "value": {"from": 0, "to": resolution,
                                          "include_lower": True, "include_upper": True}}}
            ]
        },
        "request_options": {
            "results_content_type": ["experimental"],
            "paginate": {"start": 0, "rows": max_results}
        },
        "return_type": "entry"
    }
    result = rcsb_search(query)
    if result and "result_set" in result:
        return [r["identifier"] for r in result["result_set"]]
    return []

def search_xray_by_temp(temp_low, temp_high, resolution=2.5, max_results=3000):
    """Search for X-ray structures within a temperature range."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "exptl.method",
                                "operator": "exact_match",
                                "value": "X-RAY DIFFRACTION"}},
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "diffrn.ambient_temp",
                                "operator": "range",
                                "value": {"from": temp_low, "to": temp_high,
                                          "include_lower": True, "include_upper": True}}},
                {"type": "terminal", "service": "text",
                 "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                                "operator": "range",
                                "value": {"from": 0, "to": resolution,
                                          "include_lower": True, "include_upper": True}}}
            ]
        },
        "request_options": {
            "results_content_type": ["experimental"],
            "paginate": {"start": 0, "rows": max_results}
        },
        "return_type": "entry"
    }
    result = rcsb_search(query)
    if result and "result_set" in result:
        return [r["identifier"] for r in result["result_set"]]
    return []

# ═══════════════════════════════════════════════════════════════════
# RCSB REST API — entry only
# ═══════════════════════════════════════════════════════════════════

def get_entry_info(pdb_id):
    """Get resolution, UniProt IDs, and temperature from REST entry endpoint."""
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        with urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read())
    except:
        return None

    info = {"pdb": pdb_id}

    # Resolution
    try:
        r = data.get("rcsb_entry_info", {}).get("resolution_combined", [])
        info["resolution"] = r[0] if isinstance(r, list) and r else (r if isinstance(r, (int, float)) else None)
    except:
        info["resolution"] = None

    # Temperature
    try:
        d = data.get("diffrn", [])
        if isinstance(d, list) and d:
            info["temp"] = d[0].get("ambient_temp")
        else:
            info["temp"] = None
    except:
        info["temp"] = None

    # Entity count info (to know how many polymer entities exist)
    try:
        info["polymer_entity_count"] = data.get("rcsb_entry_info", {}).get("polymer_entity_count_protein", 0)
    except:
        info["polymer_entity_count"] = 0

    return info


def get_uniprot_for_entity(pdb_id, entity_id):
    """Get UniProt mapping for a specific polymer entity."""
    try:
        url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
        with urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read())
        container = data.get("rcsb_polymer_entity_container_identifiers", {})
        uids = container.get("uniprot_ids", [])
        chains = container.get("auth_asym_ids", [])
        seq = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        etype = data.get("entity_poly", {}).get("rcsb_entity_polymer_type", "")
        return {
            "uniprot": uids[0] if uids else None,
            "chains": chains,
            "seq_len": len(seq) if seq else 0,
            "type": etype,
        }
    except:
        return None

# ═══════════════════════════════════════════════════════════════════
# PDB FILE PARSING
# ═══════════════════════════════════════════════════════════════════

def download_pdb(pdb_id):
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    fp = PDB_DIR / f"{pdb_id.upper()}.pdb"
    if fp.exists():
        return fp
    try:
        gz = PDB_DIR / f"{pdb_id.upper()}.pdb.gz"
        urlretrieve(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb.gz", gz)
        with gzip.open(gz, 'rb') as f:
            fp.write_bytes(f.read())
        gz.unlink()
        return fp
    except:
        try:
            urlretrieve(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", fp)
            return fp
        except:
            return None


def parse_pdb_info(filepath):
    """Parse PDB file for chains, residue counts, and HETATM ligands."""
    chains = defaultdict(set)  # chain_id -> set of residue numbers
    hetatm_ligands = set()

    with open(filepath) as f:
        for line in f:
            if line.startswith("ATOM"):
                ch = line[21]
                try:
                    rn = int(line[22:26].strip())
                    chains[ch].add(rn)
                except:
                    pass
            elif line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname.upper() not in SOLVENT_LIGANDS:
                    hetatm_ligands.add(resname)

    return {
        "chains": {ch: len(residues) for ch, residues in chains.items()},
        "has_ligand": len(hetatm_ligands) > 0,
        "ligands": list(hetatm_ligands),
    }


def parse_backbone(filepath, chain_id):
    residues = {}
    with open(filepath) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            an = line[12:16].strip()
            if an not in ("N", "CA", "C"):
                continue
            if line[21] != chain_id:
                continue
            try:
                rn = int(line[22:26].strip())
            except:
                continue
            if line[26].strip():
                continue
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            if rn not in residues:
                residues[rn] = {}
            residues[rn][an] = np.array([x, y, z])
    return residues

# ═══════════════════════════════════════════════════════════════════
# CURVATURE
# ═══════════════════════════════════════════════════════════════════

def dihedral(p1, p2, p3, p4):
    b1, b2, b3 = p2-p1, p3-p2, p4-p3
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    nn1, nn2 = np.linalg.norm(n1), np.linalg.norm(n2)
    if nn1 < 1e-10 or nn2 < 1e-10:
        return np.nan
    n1, n2 = n1/nn1, n2/nn2
    m1 = np.cross(n1, b2/np.linalg.norm(b2))
    return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))

def get_dihedrals(residues):
    sr = sorted(residues.keys())
    out = {}
    for i, rn in enumerate(sr):
        r = residues[rn]
        if not all(k in r for k in ("N", "CA", "C")):
            continue
        phi = np.nan
        if i > 0:
            prn = sr[i-1]; pr = residues.get(prn, {})
            if "C" in pr and (rn - prn) == 1:
                phi = dihedral(pr["C"], r["N"], r["CA"], r["C"])
        psi = np.nan
        if i < len(sr) - 1:
            nrn = sr[i+1]; nr = residues.get(nrn, {})
            if "N" in nr and (nrn - rn) == 1:
                psi = dihedral(r["N"], r["CA"], r["C"], nr["N"])
        out[rn] = (phi, psi)
    return out

def adiff(a, b):
    return np.arctan2(np.sin(a-b), np.cos(a-b))

def get_kappa_sq(dih):
    sr = sorted(dih.keys())
    ksq = {}
    for i in range(1, len(sr)-1):
        rp, rc, rn = sr[i-1], sr[i], sr[i+1]
        if (rc-rp) != 1 or (rn-rc) != 1:
            continue
        pp, sp = dih[rp]; pc, sc = dih[rc]; pn, sn = dih[rn]
        if any(np.isnan(x) for x in [pp, sp, pc, sc, pn, sn]):
            continue
        dp1, ds1 = adiff(pc, pp), adiff(sc, sp)
        dp2, ds2 = adiff(pn, pc), adiff(sn, sc)
        a1 = np.sqrt(dp1**2 + ds1**2)
        a2 = np.sqrt(dp2**2 + ds2**2)
        if a1 < 1e-10 or a2 < 1e-10:
            continue
        dTp = dp2/a2 - dp1/a1
        dTs = ds2/a2 - ds1/a1
        dm = 0.5*(a1+a2)
        ksq[rc] = (dTp**2 + dTs**2) / (dm**2)
    return ksq

def gini(values):
    v = np.sort(np.array(values, dtype=float))
    v = v[~np.isnan(v)]
    if len(v) == 0: return np.nan
    n = len(v); s = np.sum(v)
    if s < 1e-15: return 0.0
    return (2*np.sum(np.arange(1,n+1)*v) / (n*s)) - (n+1)/n

def compute_gini(filepath, chain_id):
    res = parse_backbone(filepath, chain_id)
    if len(res) < 10: return None
    dih = get_dihedrals(res)
    ksq = get_kappa_sq(dih)
    if len(ksq) < 5: return None
    return {"gini": gini(list(ksq.values())), "n_res": len(res), "n_kappa": len(ksq)}

# ═══════════════════════════════════════════════════════════════════
# METADATA COLLECTION — Two-phase approach
# ═══════════════════════════════════════════════════════════════════

def collect_metadata(pdb_ids, need_temp=False):
    """
    Phase 1: Download PDB files and parse for chains + ligands (fast, no API)
    Phase 2: Get UniProt mapping via REST entity endpoint (slower, but targeted)
    """
    print(f"  Phase 1: Download and parse {len(pdb_ids)} PDB files...")

    entries = {}
    failed = 0
    for idx, pdb_id in enumerate(pdb_ids):
        if idx % 200 == 0 and idx > 0:
            print(f"    Downloaded {idx}/{len(pdb_ids)} ({len(entries)} valid, {failed} failed)...")

        filepath = download_pdb(pdb_id)
        if not filepath:
            failed += 1
            continue

        pdb_info = parse_pdb_info(filepath)
        if not pdb_info["chains"]:
            failed += 1
            continue

        # Get resolution (and temp if needed) from REST
        entry_info = get_entry_info(pdb_id)
        if not entry_info:
            failed += 1
            continue

        if need_temp and entry_info.get("temp") is None:
            failed += 1
            continue

        # Find the longest protein chain
        best_chain = max(pdb_info["chains"].items(), key=lambda x: x[1])
        if best_chain[1] < 50 or best_chain[1] > 500:
            continue

        entries[pdb_id] = {
            "filepath": str(filepath),
            "chain": best_chain[0],
            "n_res": best_chain[1],
            "has_ligand": pdb_info["has_ligand"],
            "ligands": pdb_info["ligands"],
            "resolution": entry_info.get("resolution"),
            "temp": entry_info.get("temp"),
            "entity_count": entry_info.get("polymer_entity_count", 0),
        }

        time.sleep(0.15)  # rate limit for REST calls

    print(f"    Downloaded {len(pdb_ids)}/{len(pdb_ids)} ({len(entries)} valid, {failed} failed)")

    # Phase 2: Get UniProt IDs
    print(f"  Phase 2: Fetch UniProt mappings for {len(entries)} entries...")
    entries_with_uniprot = {}
    for idx, (pdb_id, info) in enumerate(entries.items()):
        if idx % 200 == 0 and idx > 0:
            print(f"    UniProt: {idx}/{len(entries)} ({len(entries_with_uniprot)} mapped)...")

        # Try entity 1 first (most common)
        for eid in range(1, min(info["entity_count"] + 1, 5)):
            edata = get_uniprot_for_entity(pdb_id, eid)
            if not edata or not edata.get("uniprot"):
                continue
            if edata.get("type", "") != "Protein":
                continue
            if info["chain"] in edata.get("chains", []):
                info["uniprot"] = edata["uniprot"]
                entries_with_uniprot[pdb_id] = info
                break
            # If chain not matched, still use if it's the only protein entity
            if eid == 1 and edata.get("seq_len", 0) >= 50:
                info["uniprot"] = edata["uniprot"]
                info["chain"] = edata["chains"][0] if edata.get("chains") else info["chain"]
                entries_with_uniprot[pdb_id] = info
                break

        time.sleep(0.15)

    print(f"  Got UniProt for {len(entries_with_uniprot)} entries")
    return entries_with_uniprot


# ═══════════════════════════════════════════════════════════════════
# APO-HOLO PAIRS
# ═══════════════════════════════════════════════════════════════════

def build_apo_holo_pairs(entries, max_pairs=200):
    by_uniprot = defaultdict(lambda: {"apo": [], "holo": []})
    for pdb_id, info in entries.items():
        key = "holo" if info["has_ligand"] else "apo"
        by_uniprot[info["uniprot"]][key].append({"pdb": pdb_id, **info})

    pairs = []
    for uniprot, groups in by_uniprot.items():
        if not groups["apo"] or not groups["holo"]:
            continue
        best_apo = min(groups["apo"], key=lambda x: x.get("resolution") or 99)
        best_holo = min(groups["holo"], key=lambda x: x.get("resolution") or 99)
        if best_apo["pdb"] == best_holo["pdb"]:
            continue
        pairs.append({
            "uniprot": uniprot,
            "apo_pdb": best_apo["pdb"], "apo_chain": best_apo["chain"],
            "holo_pdb": best_holo["pdb"], "holo_chain": best_holo["chain"],
            "apo_filepath": best_apo["filepath"], "holo_filepath": best_holo["filepath"],
            "apo_res": best_apo.get("resolution"), "holo_res": best_holo.get("resolution"),
        })
    pairs.sort(key=lambda x: ((x["apo_res"] or 99) + (x["holo_res"] or 99)) / 2)
    return pairs[:max_pairs]


# ═══════════════════════════════════════════════════════════════════
# TEMPERATURE PAIRS
# ═══════════════════════════════════════════════════════════════════

def build_temp_pairs(cryo_entries, rt_entries, max_pairs=200):
    by_uniprot = defaultdict(lambda: {"cryo": [], "rt": []})
    for pdb_id, info in cryo_entries.items():
        by_uniprot[info["uniprot"]]["cryo"].append({"pdb": pdb_id, **info})
    for pdb_id, info in rt_entries.items():
        by_uniprot[info["uniprot"]]["rt"].append({"pdb": pdb_id, **info})

    pairs = []
    for uniprot, groups in by_uniprot.items():
        if not groups["cryo"] or not groups["rt"]:
            continue
        best_cryo = min(groups["cryo"], key=lambda x: x.get("resolution") or 99)
        best_rt = min(groups["rt"], key=lambda x: x.get("resolution") or 99)
        if best_cryo["pdb"] == best_rt["pdb"]:
            continue
        pairs.append({
            "uniprot": uniprot,
            "cryo_pdb": best_cryo["pdb"], "cryo_chain": best_cryo["chain"],
            "rt_pdb": best_rt["pdb"], "rt_chain": best_rt["chain"],
            "cryo_filepath": best_cryo["filepath"], "rt_filepath": best_rt["filepath"],
            "cryo_temp": best_cryo.get("temp"), "rt_temp": best_rt.get("temp"),
            "cryo_res": best_cryo.get("resolution"), "rt_res": best_rt.get("resolution"),
        })
    pairs.sort(key=lambda x: ((x["cryo_res"] or 99) + (x["rt_res"] or 99)) / 2)
    return pairs[:max_pairs]


# ═══════════════════════════════════════════════════════════════════
# ANALYSIS (shared)
# ═══════════════════════════════════════════════════════════════════

def analyze_gradient(results, label_base, label_pert, results_file):
    g_base = np.array([r["g_base"] for r in results])
    g_pert = np.array([r["g_pert"] for r in results])
    dgs = np.array([r["dg"] for r in results])

    print(f"\n  Dataset: {len(results)} pairs")
    print(f"  {label_base} Gini: mean={np.mean(g_base):.4f} ± {np.std(g_base):.4f}  "
          f"range=[{np.min(g_base):.4f}, {np.max(g_base):.4f}]")
    print(f"  {label_pert} Gini: mean={np.mean(g_pert):.4f} ± {np.std(g_pert):.4f}  "
          f"range=[{np.min(g_pert):.4f}, {np.max(g_pert):.4f}]")
    print(f"  ΔGini:     mean={np.mean(dgs):+.4f} ± {np.std(dgs):.4f}")
    print(f"  % stiffening (ΔG > 0): {100*np.mean(dgs > 0):.1f}%")
    print(f"  % loosening  (ΔG < 0): {100*np.mean(dgs < 0):.1f}%")

    # Variance convergence
    var_base = np.var(g_base)
    var_pert = np.var(g_pert)
    vr = var_pert / var_base if var_base > 0 else float('inf')
    print(f"\n  Variance({label_base}): {var_base:.6f}")
    print(f"  Variance({label_pert}): {var_pert:.6f}")
    print(f"  Ratio: {vr:.4f} {'← converging' if vr < 1 else '← diverging'}")

    # Correlation
    corr = np.corrcoef(g_base, dgs)[0, 1]
    print(f"\n  Pearson r(G_base, ΔG) = {corr:.4f}")

    # Threshold scan
    best_score = 0; best_ths = []
    for ti in range(600, 950):
        th = ti/1000
        c = sum(1 for r in results if (("+" if r["g_base"]<th else "-") == ("+" if r["dg"]>0 else "-")))
        if c > best_score: best_score = c; best_ths = [th]
        elif c == best_score: best_ths.append(th)
    mid = (best_ths[0]+best_ths[-1])/2
    print(f"\n  Best accuracy: {best_score}/{len(results)} ({100*best_score/len(results):.1f}%)")
    print(f"  Threshold: [{best_ths[0]:.3f}, {best_ths[-1]:.3f}] midpoint={mid:.3f}")

    for label, th in [("0.810 (original)", 0.810), ("0.835 (AlphaFold)", 0.835),
                       ("0.860 (pilot)", 0.860), ("Best midpoint", mid)]:
        c = sum(1 for r in results if (("+" if r["g_base"]<th else "-") == ("+" if r["dg"]>0 else "-")))
        print(f"    {label:>25s}: {c}/{len(results)} ({100*c/len(results):.1f}%)")

    # Binned
    print(f"\n  {'Bin':>12s} {'N':>5s} {'mean ΔG':>8s} {'%stiff':>7s} {'%loose':>7s}")
    print("  " + "─"*42)
    for lo, hi in [(0.5,.70),(.70,.75),(.75,.80),(.80,.85),(.85,.90),(.90,.95),(.95,1.0)]:
        m = (g_base >= lo) & (g_base < hi); n = np.sum(m)
        if n < 2: continue
        print(f"  [{lo:.2f},{hi:.2f}) {n:>5d} {np.mean(dgs[m]):>+8.4f} {100*np.mean(dgs[m]>0):>6.1f}% {100*np.mean(dgs[m]<0):>6.1f}%")

    # Distance-based accuracy
    print(f"\n  Accuracy by distance from threshold:")
    dist = np.abs(g_base - mid)
    for d_lo, d_hi, dl in [(0,.02,"near"),(0.02,.05,"mid"),(0.05,.10,"far"),(0.10,.50,"v.far")]:
        m = (dist>=d_lo)&(dist<d_hi); n = np.sum(m)
        if n < 2: continue
        c = sum(1 for i,r in enumerate(results) if m[i] and (("+" if r["g_base"]<mid else "-")==("+" if r["dg"]>0 else "-")))
        print(f"    {dl:>6s}: {c}/{n} ({100*c/n:.1f}%)")

    # Save
    export = {"n_pairs": len(results), "best_accuracy": best_score/len(results),
              "threshold": mid, "pearson_r": corr, "var_ratio": vr,
              "mean_g_base": float(np.mean(g_base)), "mean_dg": float(np.mean(dgs)),
              "results": results}
    with open(results_file, "w") as f:
        json.dump(export, f, indent=2)
    print(f"\n  Saved: {results_file}")
    return mid


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-pairs", type=int, default=200)
    parser.add_argument("--max-search", type=int, default=3000)
    parser.add_argument("--resolution", type=float, default=2.0)
    parser.add_argument("--mode", choices=["apo-holo", "temperature", "both"], default="both")
    args = parser.parse_args()

    DATA_DIR.mkdir(exist_ok=True)
    PDB_DIR.mkdir(parents=True, exist_ok=True)

    thresholds = {}

    # ── APO-HOLO TEST ──
    if args.mode in ("apo-holo", "both"):
        print("╔══════════════════════════════════════════════════════════════╗")
        print("║  APO-HOLO GRADIENT TEST                                    ║")
        print(f"║  Target: {args.max_pairs} pairs, resolution ≤ {args.resolution} Å               ║")
        print("╚══════════════════════════════════════════════════════════════╝")

        pairs_file = DATA_DIR / "apo_holo_pairs.json"
        if pairs_file.exists():
            print("\n  Loading cached pairs...")
            with open(pairs_file) as f: pairs = json.load(f)
        else:
            pdb_ids = search_xray(resolution=args.resolution, max_results=args.max_search)
            print(f"  Found {len(pdb_ids)} structures")
            if not pdb_ids: sys.exit(1)

            entries = collect_metadata(pdb_ids, need_temp=False)
            pairs = build_apo_holo_pairs(entries, max_pairs=args.max_pairs)
            print(f"  Built {len(pairs)} apo-holo pairs")
            with open(pairs_file, "w") as f: json.dump(pairs, f, indent=2)

        # Compute Gini
        print(f"\n── Computing Gini for {len(pairs)} apo-holo pairs ──")
        results = []
        for idx, p in enumerate(pairs):
            if idx % 50 == 0: print(f"  {idx}/{len(pairs)}...")
            g_apo = compute_gini(p["apo_filepath"], p["apo_chain"])
            g_holo = compute_gini(p["holo_filepath"], p["holo_chain"])
            if g_apo and g_holo:
                results.append({
                    "uniprot": p["uniprot"], "base_pdb": p["apo_pdb"], "pert_pdb": p["holo_pdb"],
                    "g_base": g_apo["gini"], "g_pert": g_holo["gini"],
                    "dg": g_holo["gini"] - g_apo["gini"], "n_res": g_apo["n_res"],
                })

        print(f"\n{'='*60}")
        print("APO-HOLO RESULTS")
        print('='*60)
        th = analyze_gradient(results, "Apo", "Holo", DATA_DIR / "apo_holo_results.json")
        thresholds["apo_holo"] = th

    # ── TEMPERATURE TEST ──
    if args.mode in ("temperature", "both"):
        print("\n╔══════════════════════════════════════════════════════════════╗")
        print("║  TEMPERATURE GRADIENT TEST                                 ║")
        print("║  Cryo (100K) → Room Temperature (295K)                    ║")
        print("╚══════════════════════════════════════════════════════════════╝")

        pairs_file = DATA_DIR / "temp_pairs.json"
        if pairs_file.exists():
            print("\n  Loading cached pairs...")
            with open(pairs_file) as f: pairs = json.load(f)
        else:
            print("\n── Searching for RT structures (270-330K) ──")
            rt_ids = search_xray_by_temp(270, 330, resolution=2.5, max_results=args.max_search)
            print(f"  Found {len(rt_ids)} RT structures")

            print("\n── Searching for cryo structures (80-120K) ──")
            cryo_ids = search_xray_by_temp(80, 120, resolution=2.5, max_results=args.max_search)
            print(f"  Found {len(cryo_ids)} cryo structures")

            if not rt_ids or not cryo_ids: sys.exit(1)

            print("\n── Collecting RT metadata ──")
            rt_entries = collect_metadata(rt_ids, need_temp=True)
            print(f"\n── Collecting cryo metadata ──")
            cryo_entries = collect_metadata(cryo_ids, need_temp=True)

            pairs = build_temp_pairs(cryo_entries, rt_entries, max_pairs=args.max_pairs)
            print(f"  Built {len(pairs)} cryo-RT pairs")
            with open(pairs_file, "w") as f: json.dump(pairs, f, indent=2)

        # Compute Gini
        print(f"\n── Computing Gini for {len(pairs)} temperature pairs ──")
        results = []
        for idx, p in enumerate(pairs):
            if idx % 50 == 0: print(f"  {idx}/{len(pairs)}...")
            g_cryo = compute_gini(p["cryo_filepath"], p["cryo_chain"])
            g_rt = compute_gini(p["rt_filepath"], p["rt_chain"])
            if g_cryo and g_rt:
                results.append({
                    "uniprot": p["uniprot"], "base_pdb": p["cryo_pdb"], "pert_pdb": p["rt_pdb"],
                    "g_base": g_cryo["gini"], "g_pert": g_rt["gini"],
                    "dg": g_rt["gini"] - g_cryo["gini"], "n_res": g_cryo["n_res"],
                    "cryo_temp": p.get("cryo_temp"), "rt_temp": p.get("rt_temp"),
                })

        print(f"\n{'='*60}")
        print("TEMPERATURE RESULTS")
        print('='*60)
        th = analyze_gradient(results, "Cryo", "RT", DATA_DIR / "temperature_results.json")
        thresholds["temperature"] = th

    # ── CROSS-COMPARISON ──
    if len(thresholds) == 2:
        print(f"\n{'='*60}")
        print("CROSS-COMPARISON")
        print('='*60)
        print(f"  Apo-holo threshold:    {thresholds['apo_holo']:.3f}")
        print(f"  Temperature threshold: {thresholds['temperature']:.3f}")
        diff = abs(thresholds['apo_holo'] - thresholds['temperature'])
        print(f"  Difference:            {diff:.3f}")
        if diff < 0.03:
            print(f"  → CONVERGENT: Both perturbations point to same attractor")
        else:
            print(f"  → DIVERGENT: Different perturbations find different thresholds")


if __name__ == "__main__":
    main()
