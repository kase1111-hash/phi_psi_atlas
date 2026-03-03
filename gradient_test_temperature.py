#!/usr/bin/env python3
"""
Large-Scale Temperature Gradient Test
======================================
Tests whether thermal energy drives Gini toward the attractor.

Prediction: heating (cryo → RT) should push Gini toward ~0.86:
  - Proteins with G_cryo > attractor should LOOSEN (ΔG < 0)
  - Proteins with G_cryo < attractor should STIFFEN (ΔG > 0)

Strategy:
  1. Search RCSB for room-temperature X-ray structures (T > 270K)
  2. For each, find the same protein at cryo (T ~ 100K) via UniProt
  3. Compute Gini for both, measure ΔGini

This is cleaner than apo-holo because:
  - No ligand ambiguity (same chemical state)
  - No binding mechanics (just thermal energy)
  - Same crystal form in many cases
  - The perturbation (temperature) is a scalar with known direction

Usage:
    python gradient_test_temperature.py [--max-pairs 200]
    
    Requires: numpy, urllib (stdlib). Network access to RCSB.

Author: Kase / True North Construction LLC
License: CC0 1.0 Universal (unpatentable prior art)
"""

import numpy as np
import os, sys, json, time, gzip
from urllib.request import urlopen, Request, urlretrieve
from urllib.error import URLError, HTTPError
from pathlib import Path
from collections import defaultdict
import argparse

DATA_DIR = Path("gradient_temperature")
PDB_DIR = DATA_DIR / "pdbs"
RESULTS_FILE = DATA_DIR / "temperature_results.json"

# ═══════════════════════════════════════════════════════════════════
# RCSB PDB SEARCH API
# ═══════════════════════════════════════════════════════════════════

def rcsb_search(query_json):
    """Execute RCSB PDB search API query."""
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    req = Request(url, data=json.dumps(query_json).encode(),
                  headers={"Content-Type": "application/json"})
    try:
        with urlopen(req, timeout=60) as resp:
            return json.loads(resp.read())
    except HTTPError as e:
        body = ""
        try:
            body = e.read().decode('utf-8', errors='replace')[:500]
        except:
            pass
        print(f"  Search API HTTP {e.code}: {e.reason}")
        if body:
            print(f"  Response body: {body}")
        return None
    except Exception as e:
        print(f"  Search API error: {e}")
        return None


def search_by_temperature(temp_low, temp_high, max_results=2000, resolution=2.5):
    """Search for X-ray structures within a temperature range."""
    
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "diffrn.ambient_temp",
                        "operator": "range",
                        "value": {
                            "from": temp_low,
                            "to": temp_high,
                            "include_lower": True,
                            "include_upper": True
                        }
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "range",
                        "value": {
                            "from": 0,
                            "to": resolution,
                            "include_lower": True,
                            "include_upper": True
                        }
                    }
                }
            ]
        },
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "paginate": {"start": 0, "rows": max_results}
        },
        "return_type": "entry"
    }
    
    result = rcsb_search(query)
    if result and "result_set" in result:
        return [r["identifier"] for r in result["result_set"]]
    
    # Fallback: simpler query without resolution
    print("  Primary query failed, trying without resolution filter...")
    query["query"]["nodes"] = query["query"]["nodes"][:2]  # keep only method + temp
    result = rcsb_search(query)
    if result and "result_set" in result:
        return [r["identifier"] for r in result["result_set"]]
    
    return []


def graphql_batch_with_temp(pdb_ids):
    """Fetch entry data including temperature via GraphQL."""
    url = "https://data.rcsb.org/graphql"
    query = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        diffrn {
          ambient_temp
        }
        rcsb_entry_info {
          resolution_combined
        }
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            uniprot_ids
            auth_asym_ids
          }
          entity_poly {
            pdbx_seq_one_letter_code_can
            rcsb_entity_polymer_type
          }
        }
      }
    }
    """
    body = json.dumps({"query": query, "variables": {"ids": pdb_ids}})
    req = Request(url, data=body.encode(),
                  headers={"Content-Type": "application/json"})
    try:
        with urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
        return data.get("data", {}).get("entries", [])
    except Exception as e:
        print(f"    GraphQL error: {e}")
        return []


def get_entries_metadata_batch(pdb_ids, batch_size=40):
    """Batch-fetch metadata for many PDB entries via GraphQL."""
    print(f"  Fetching metadata for {len(pdb_ids)} structures via GraphQL...")
    
    results = []
    
    for batch_start in range(0, len(pdb_ids), batch_size):
        batch = pdb_ids[batch_start:batch_start + batch_size]
        entries = graphql_batch_with_temp(batch)
        
        for entry_data in entries:
            if not entry_data:
                continue
            pdb_id = entry_data.get("rcsb_id", "")
            if not pdb_id:
                continue
            
            # Temperature
            temp = None
            try:
                diffrn = entry_data.get("diffrn") or []
                if isinstance(diffrn, list) and diffrn:
                    temp = diffrn[0].get("ambient_temp")
                elif isinstance(diffrn, dict):
                    temp = diffrn.get("ambient_temp")
            except:
                pass
            
            # Resolution
            resolution = None
            try:
                res_list = (entry_data.get("rcsb_entry_info") or {}).get("resolution_combined", [])
                if isinstance(res_list, list) and res_list:
                    resolution = res_list[0]
                elif isinstance(res_list, (int, float)):
                    resolution = res_list
            except:
                pass
            
            # Chains
            chains = []
            poly_entities = entry_data.get("polymer_entities") or []
            for entity in poly_entities:
                if not entity:
                    continue
                try:
                    etype = (entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type", "")
                    if etype and etype != "Protein":
                        continue
                except:
                    pass
                
                uniprot = None
                try:
                    container = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                    uids = container.get("uniprot_ids") or []
                    if uids:
                        uniprot = uids[0]
                except:
                    pass
                
                chain_ids = []
                try:
                    container = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                    chain_ids = container.get("auth_asym_ids") or []
                except:
                    pass
                
                seq_len = 0
                try:
                    seq = (entity.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can", "")
                    seq_len = len(seq) if seq else 0
                except:
                    pass
                
                if uniprot and chain_ids and 50 <= seq_len <= 500:
                    chains.append({
                        "uniprot": uniprot,
                        "chain": chain_ids[0],
                        "seq_len": seq_len,
                    })
            
            if chains and temp is not None:
                results.append({
                    "pdb": pdb_id,
                    "temp": float(temp),
                    "resolution": resolution,
                    "chains": chains,
                })
        
        count = min(batch_start + batch_size, len(pdb_ids))
        if count % 200 == 0 or count == len(pdb_ids):
            print(f"    Processed {count}/{len(pdb_ids)} ({len(results)} valid)...")
        time.sleep(0.3)
    
    print(f"  Got {len(results)} valid entries")
    return results


# ═══════════════════════════════════════════════════════════════════
# PAIR BUILDING
# ═══════════════════════════════════════════════════════════════════

def build_cryo_rt_pairs(rt_entries, cryo_entries, max_pairs=300):
    """
    Match RT and cryo structures by UniProt ID.
    Pick best resolution for each temperature.
    """
    
    # Index by UniProt
    by_uniprot = defaultdict(lambda: {"rt": [], "cryo": []})
    
    for entry in rt_entries:
        for chain in entry["chains"]:
            by_uniprot[chain["uniprot"]]["rt"].append({
                "pdb": entry["pdb"],
                "chain": chain["chain"],
                "temp": entry["temp"],
                "resolution": entry["resolution"] or 99.0,
                "seq_len": chain["seq_len"],
            })
    
    for entry in cryo_entries:
        for chain in entry["chains"]:
            by_uniprot[chain["uniprot"]]["cryo"].append({
                "pdb": entry["pdb"],
                "chain": chain["chain"],
                "temp": entry["temp"],
                "resolution": entry["resolution"] or 99.0,
                "seq_len": chain["seq_len"],
            })
    
    # Build pairs
    pairs = []
    for uniprot, entries in by_uniprot.items():
        if not entries["rt"] or not entries["cryo"]:
            continue
        
        best_rt = min(entries["rt"], key=lambda x: x["resolution"])
        best_cryo = min(entries["cryo"], key=lambda x: x["resolution"])
        
        if best_rt["pdb"] == best_cryo["pdb"]:
            continue
        
        pairs.append({
            "uniprot": uniprot,
            "cryo_pdb": best_cryo["pdb"],
            "cryo_chain": best_cryo["chain"],
            "cryo_temp": best_cryo["temp"],
            "cryo_res": best_cryo["resolution"],
            "rt_pdb": best_rt["pdb"],
            "rt_chain": best_rt["chain"],
            "rt_temp": best_rt["temp"],
            "rt_res": best_rt["resolution"],
            "seq_len": best_cryo["seq_len"],
        })
    
    # Sort by average resolution
    pairs.sort(key=lambda x: (x["cryo_res"] + x["rt_res"]) / 2)
    return pairs[:max_pairs]


# ═══════════════════════════════════════════════════════════════════
# PDB + CURVATURE (same as apo-holo script)
# ═══════════════════════════════════════════════════════════════════

def download_pdb(pdb_id, data_dir=PDB_DIR):
    data_dir.mkdir(parents=True, exist_ok=True)
    filepath = data_dir / f"{pdb_id.upper()}.pdb"
    if filepath.exists():
        return filepath
    gz_url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb.gz"
    try:
        gz_path = data_dir / f"{pdb_id.upper()}.pdb.gz"
        urlretrieve(gz_url, gz_path)
        with gzip.open(gz_path, 'rb') as gz:
            with open(filepath, 'wb') as out:
                out.write(gz.read())
        gz_path.unlink()
        return filepath
    except:
        pass
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        urlretrieve(url, filepath)
        return filepath
    except:
        return None


def parse_backbone(filepath, chain_id="A"):
    residues = {}
    with open(filepath) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            an = line[12:16].strip()
            if an not in ("N", "CA", "C"):
                continue
            ch = line[21].strip()
            if ch != chain_id:
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
    if len(v) == 0:
        return np.nan
    n = len(v)
    s = np.sum(v)
    if s < 1e-15:
        return 0.0
    return (2*np.sum(np.arange(1,n+1)*v) / (n*s)) - (n+1)/n

def compute_gini(filepath, chain_id):
    res = parse_backbone(filepath, chain_id)
    if len(res) < 10:
        return None
    dih = get_dihedrals(res)
    ksq = get_kappa_sq(dih)
    if len(ksq) < 5:
        return None
    g = gini(list(ksq.values()))
    return {"gini": g, "n_res": len(res), "n_kappa": len(ksq)}


# ═══════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Temperature gradient test")
    parser.add_argument("--max-pairs", type=int, default=200)
    parser.add_argument("--max-search", type=int, default=3000)
    parser.add_argument("--resolution", type=float, default=2.5,
                        help="Resolution cutoff (Å) — wider than apo-holo since RT data is lower res")
    parser.add_argument("--skip-search", action="store_true")
    args = parser.parse_args()
    
    DATA_DIR.mkdir(exist_ok=True)
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  TEMPERATURE GRADIENT TEST                                     ║")
    print("║  Cryo (100K) → Room Temperature (295K)                        ║")
    print(f"║  Target: {args.max_pairs} pairs, resolution ≤ {args.resolution} Å                       ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    
    pairs_file = DATA_DIR / "cryo_rt_pairs.json"
    
    if args.skip_search and pairs_file.exists():
        print("\n── Loading cached pairs ──")
        with open(pairs_file) as f:
            pairs = json.load(f)
        print(f"  Loaded {len(pairs)} pairs")
    else:
        # Step 1: Find room-temperature structures
        print("\n── Step 1: Search for room-temperature structures (270-330K) ──")
        rt_ids = search_by_temperature(270, 330, max_results=args.max_search,
                                        resolution=args.resolution)
        print(f"  Found {len(rt_ids)} RT structures")
        
        if not rt_ids:
            print("No RT structures found. Check network.")
            sys.exit(1)
        
        # Step 2: Get metadata for RT structures (batched)
        print("\n── Step 2: Fetch RT structure metadata ──")
        rt_entries = get_entries_metadata_batch(rt_ids)
        
        # Step 3: Collect UniProt IDs from RT entries
        rt_uniprots = set()
        for entry in rt_entries:
            for chain in entry["chains"]:
                rt_uniprots.add(chain["uniprot"])
        print(f"  {len(rt_uniprots)} unique UniProt accessions in RT set")
        
        # Step 4: Find cryo structures for these same proteins
        print("\n── Step 3: Search for cryo structures (80-120K) ──")
        cryo_ids = search_by_temperature(80, 120, max_results=min(args.max_search * 2, 5000),
                                          resolution=args.resolution)
        print(f"  Found {len(cryo_ids)} cryo structures")
        
        # Step 5: Get metadata for cryo structures (batched)
        print("\n── Step 4: Fetch cryo structure metadata ──")
        cryo_entries = get_entries_metadata_batch(cryo_ids)
        
        # Step 6: Build pairs
        print("\n── Step 5: Build cryo-RT pairs ──")
        pairs = build_cryo_rt_pairs(rt_entries, cryo_entries, max_pairs=args.max_pairs)
        print(f"  Built {len(pairs)} matched pairs")
        
        with open(pairs_file, "w") as f:
            json.dump(pairs, f, indent=2)
    
    # Step 7: Download and compute
    print(f"\n── Step 6: Download PDBs and compute Gini ({len(pairs)} pairs) ──")
    
    results = []
    failed = 0
    
    for idx, pair in enumerate(pairs):
        if idx % 20 == 0:
            print(f"  Processing pair {idx+1}/{len(pairs)}...")
        
        cryo_path = download_pdb(pair["cryo_pdb"])
        rt_path = download_pdb(pair["rt_pdb"])
        
        if not cryo_path or not rt_path:
            failed += 1
            continue
        
        cryo_result = compute_gini(cryo_path, pair["cryo_chain"])
        rt_result = compute_gini(rt_path, pair["rt_chain"])
        
        if cryo_result is None or rt_result is None:
            failed += 1
            continue
        
        dg = rt_result["gini"] - cryo_result["gini"]
        
        results.append({
            "uniprot": pair["uniprot"],
            "cryo_pdb": pair["cryo_pdb"],
            "rt_pdb": pair["rt_pdb"],
            "cryo_temp": pair["cryo_temp"],
            "rt_temp": pair["rt_temp"],
            "g_cryo": cryo_result["gini"],
            "g_rt": rt_result["gini"],
            "dg": dg,
            "n_res": cryo_result["n_res"],
        })
        
        time.sleep(0.1)
    
    print(f"\n  Completed: {len(results)} pairs ({failed} failed)")
    
    if len(results) < 10:
        print("Too few results for analysis.")
        sys.exit(1)
    
    # ═══════════════════════════════════════════════════════════════
    # ANALYSIS
    # ═══════════════════════════════════════════════════════════════
    
    print("\n" + "="*70)
    print("ANALYSIS: TEMPERATURE GRADIENT")
    print("="*70)
    
    g_cryos = np.array([r["g_cryo"] for r in results])
    g_rts = np.array([r["g_rt"] for r in results])
    dgs = np.array([r["dg"] for r in results])
    
    print(f"\n  Dataset: {len(results)} cryo-RT pairs")
    print(f"  Cryo Gini: mean={np.mean(g_cryos):.4f} ± {np.std(g_cryos):.4f}  "
          f"range=[{np.min(g_cryos):.4f}, {np.max(g_cryos):.4f}]")
    print(f"  RT Gini:   mean={np.mean(g_rts):.4f} ± {np.std(g_rts):.4f}  "
          f"range=[{np.min(g_rts):.4f}, {np.max(g_rts):.4f}]")
    print(f"  ΔGini:     mean={np.mean(dgs):+.4f} ± {np.std(dgs):.4f}")
    print(f"  % stiffening (ΔG > 0): {100*np.mean(dgs > 0):.1f}%")
    print(f"  % loosening  (ΔG < 0): {100*np.mean(dgs < 0):.1f}%")
    
    # KEY TEST: Does Gini converge toward a single value?
    print("\n── Convergence Test ──")
    print(f"  Variance(cryo):  {np.var(g_cryos):.6f}")
    print(f"  Variance(RT):    {np.var(g_rts):.6f}")
    variance_ratio = np.var(g_rts) / np.var(g_cryos)
    print(f"  Var ratio (RT/cryo): {variance_ratio:.4f}")
    if variance_ratio < 1:
        print(f"  → RT distribution is NARROWER — consistent with attractor convergence")
    else:
        print(f"  → RT distribution is WIDER — inconsistent with attractor")
    
    # Correlation: does baseline predict direction?
    print("\n── Correlation: G_cryo vs ΔGini ──")
    corr = np.corrcoef(g_cryos, dgs)[0, 1]
    print(f"  Pearson r = {corr:.4f}")
    if corr < -0.1:
        print(f"  → Negative correlation: high-Gini proteins loosen, low-Gini stiffen")
    elif corr > 0.1:
        print(f"  → Positive correlation: unexpected direction")
    else:
        print(f"  → Weak correlation")
    
    # Threshold scan — same as apo-holo
    print("\n── Threshold Scan ──")
    best_score = 0
    best_thresholds = []
    
    for ti in range(600, 950):
        th = ti / 1000
        correct = sum(1 for r in results 
                      if (("+" if r["g_cryo"] < th else "-") == ("+" if r["dg"] > 0 else "-")))
        if correct > best_score:
            best_score = correct
            best_thresholds = [th]
        elif correct == best_score:
            best_thresholds.append(th)
    
    best_acc = best_score / len(results)
    mid = (best_thresholds[0] + best_thresholds[-1]) / 2
    print(f"  Best accuracy: {best_score}/{len(results)} ({100*best_acc:.1f}%)")
    print(f"  Best threshold range: [{best_thresholds[0]:.3f}, {best_thresholds[-1]:.3f}]")
    print(f"  Midpoint: {mid:.3f}")
    
    for label, th in [("0.810 (original)", 0.810),
                       ("0.835 (AlphaFold)", 0.835),
                       ("0.860 (pilot)", 0.860),
                       ("Best midpoint", mid)]:
        correct = sum(1 for r in results 
                      if (("+" if r["g_cryo"] < th else "-") == ("+" if r["dg"] > 0 else "-")))
        print(f"    At {label:>25s}: {correct}/{len(results)} ({100*correct/len(results):.1f}%)")
    
    # Binned analysis
    print("\n── Binned Analysis ──")
    bins = [(0.5, 0.70), (0.70, 0.75), (0.75, 0.80), (0.80, 0.85), 
            (0.85, 0.90), (0.90, 0.95), (0.95, 1.0)]
    print(f"  {'Bin':>12s} {'N':>5s} {'mean ΔG':>8s} {'% stiff':>8s} {'% loose':>8s}")
    print("  " + "─"*45)
    for lo, hi in bins:
        mask = (g_cryos >= lo) & (g_cryos < hi)
        n = np.sum(mask)
        if n < 2:
            continue
        mean_dg = np.mean(dgs[mask])
        pct_stiff = 100*np.mean(dgs[mask] > 0)
        pct_loose = 100*np.mean(dgs[mask] < 0)
        print(f"  [{lo:.2f},{hi:.2f}) {n:>5d} {mean_dg:>+8.4f} {pct_stiff:>7.1f}% {pct_loose:>7.1f}%")
    
    # Effect size by distance from threshold
    print("\n── Prediction Accuracy by Distance from Threshold ──")
    dist = np.abs(g_cryos - mid)
    for d_lo, d_hi, label in [(0, 0.02, "near (0-0.02)"), 
                                (0.02, 0.05, "mid (0.02-0.05)"),
                                (0.05, 0.10, "far (0.05-0.10)"),
                                (0.10, 0.50, "very far (0.10+)")]:
        mask = (dist >= d_lo) & (dist < d_hi)
        n = np.sum(mask)
        if n < 2:
            continue
        correct = sum(1 for i, r in enumerate(results) 
                      if mask[i] and (("+" if r["g_cryo"] < mid else "-") == ("+" if r["dg"] > 0 else "-")))
        acc = correct / n
        mean_abs_dg = np.mean(np.abs(dgs[mask]))
        print(f"  {label:>20s}: {correct}/{n} ({100*acc:.1f}%) mean|ΔG|={mean_abs_dg:.4f}")
    
    # CROSS-COMPARISON: Does temperature threshold match apo-holo?
    print("\n── Cross-Comparison with Apo-Holo Test ──")
    print(f"  Temperature threshold midpoint: {mid:.3f}")
    print(f"  (Compare with apo-holo threshold from large_scale_results.json)")
    
    # Export
    export = {
        "n_pairs": len(results),
        "perturbation": "temperature",
        "cryo_range": "80-120K",
        "rt_range": "270-330K",
        "best_accuracy": best_acc,
        "best_threshold_range": [best_thresholds[0], best_thresholds[-1]],
        "midpoint": mid,
        "pearson_r": corr,
        "variance_ratio_rt_over_cryo": variance_ratio,
        "mean_g_cryo": float(np.mean(g_cryos)),
        "mean_g_rt": float(np.mean(g_rts)),
        "mean_dg": float(np.mean(dgs)),
        "results": results,
    }
    
    with open(RESULTS_FILE, "w") as f:
        json.dump(export, f, indent=2)
    print(f"\nResults exported to: {RESULTS_FILE}")
    
    csv_path = DATA_DIR / "temperature_data.csv"
    with open(csv_path, "w") as f:
        f.write("uniprot,cryo_pdb,rt_pdb,cryo_temp,rt_temp,g_cryo,g_rt,delta_gini,n_res\n")
        for r in sorted(results, key=lambda x: x["g_cryo"]):
            f.write(f"{r['uniprot']},{r['cryo_pdb']},{r['rt_pdb']},"
                    f"{r['cryo_temp']},{r['rt_temp']},"
                    f"{r['g_cryo']:.6f},{r['g_rt']:.6f},{r['dg']:.6f},"
                    f"{r['n_res']}\n")
    print(f"CSV exported to: {csv_path}")


if __name__ == "__main__":
    main()
