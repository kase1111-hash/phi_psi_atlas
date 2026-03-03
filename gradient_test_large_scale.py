#!/usr/bin/env python3
"""
Large-Scale Gradient Test: Automated Apo-Holo Pair Discovery
=============================================================
Instead of hand-picking 9 proteins, this script:
  1. Queries RCSB PDB search API for high-quality X-ray structures
  2. Groups by UniProt ID to find proteins with both apo and holo forms
  3. Computes Gini for each structure
  4. Tests gradient prediction across all pairs
  5. Reports statistics, threshold scan, and full dataset

Target: 100-500 apo-holo pairs, enough to see the real pattern.

Selection criteria:
  - X-ray only, resolution ≤ 2.0 Å (high quality)
  - Single chain, 50-400 residues (avoids multimers and giant complexes)
  - Same UniProt accession for apo and holo
  - Apo = no ligand heavier than water/ions
  - Holo = has at least one non-solvent ligand

Usage:
    python gradient_test_large_scale.py [--max-pairs 200] [--resolution 2.0]
    
    Requires: numpy, urllib/json (stdlib). Network access to RCSB.

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

DATA_DIR = Path("gradient_large_scale")
PDB_DIR = DATA_DIR / "pdbs"
RESULTS_FILE = DATA_DIR / "large_scale_results.json"

# ═══════════════════════════════════════════════════════════════════
# RCSB PDB SEARCH API
# ═══════════════════════════════════════════════════════════════════

def rcsb_search(query_json, return_type="entry"):
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


def find_apo_holo_pairs(max_results=2000, resolution=2.0, min_length=50, max_length=400):
    """
    Find proteins that have both apo (no ligand) and holo (with ligand) structures.
    
    Strategy: 
    1. Search for high-quality X-ray structures with good resolution
    2. Use RCSB data API to get UniProt mapping and ligand info
    3. Group by UniProt and identify apo vs holo
    """
    
    print(f"  Searching RCSB for X-ray structures (res ≤ {resolution} Å, "
          f"{min_length}-{max_length} res)...")
    
    # Use validated RCSB search API v2 query format
    # resolution_combined needs range operator; sequence length filtered later via metadata
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
                        "attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "range",
                        "value": {
                            "from": 0,
                            "to": resolution,
                            "include_lower": True,
                            "include_upper": True
                        }
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater_or_equal",
                        "value": 1
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entry_info.deposited_polymer_monomer_count",
                        "operator": "range",
                        "value": {
                            "from": min_length,
                            "to": max_length,
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
    if not result or "result_set" not in result:
        # Fallback: simplest possible query
        print("  Primary query failed, trying simplified query...")
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
                            "attribute": "rcsb_entry_info.resolution_combined",
                            "operator": "less_or_equal",
                            "value": resolution
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
        if not result or "result_set" not in result:
            # Third fallback: absolute minimum query
            print("  Simplified query also failed. Trying minimal query...")
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION"
                    }
                },
                "request_options": {
                    "paginate": {"start": 0, "rows": max_results}
                },
                "return_type": "entry"
            }
            result = rcsb_search(query)
            if not result or "result_set" not in result:
                print("  All queries failed. Dumping last query for debug:")
                print(json.dumps(query, indent=2))
                return []
    
    pdb_ids = [r["identifier"] for r in result["result_set"]]
    print(f"  Found {len(pdb_ids)} structures")
    return pdb_ids



def graphql_batch(pdb_ids):
    """Fetch entry + polymer entity + nonpolymer entity data via GraphQL."""
    url = "https://data.rcsb.org/graphql"
    query = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
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
        nonpolymer_entities {
          rcsb_nonpolymer_entity_container_identifiers {
            comp_id
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


def get_structure_info_batch(pdb_ids, batch_size=40):
    """
    Get UniProt mapping, chain info, and ligand presence via GraphQL batches.
    ~50x faster than per-entry REST calls.
    """
    print(f"  Fetching metadata for {len(pdb_ids)} structures via GraphQL...")
    
    SOLVENT_LIGANDS = {
        "HOH", "DOD", "WAT",
        "SO4", "PO4", "ACT", "GOL", "EDO", "PEG", "DMS", "MPD",
        "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "CU", "NI", "CO",
        "BME", "DTT", "TRS", "CIT", "TAR", "FMT", "IMD",
        "EPE", "MES", "HED", "PGE", "1PE", "P6G", "12P",
        "IOD", "BR", "F", "SCN", "NO3",
        "CD", "HG", "SR", "BA", "CS", "PB", "PT",
        "NH4", "LI", "RB",
        "AZI", "BOG", "C8E", "LDA", "SDS", "UNX",
    }

    all_info = {}
    
    for batch_start in range(0, len(pdb_ids), batch_size):
        batch = pdb_ids[batch_start:batch_start + batch_size]
        entries = graphql_batch(batch)
        
        for entry_data in entries:
            if not entry_data:
                continue
            pdb_id = entry_data.get("rcsb_id", "")
            if not pdb_id:
                continue
            
            # Resolution
            resolution = None
            try:
                res_list = entry_data.get("rcsb_entry_info", {}).get("resolution_combined", [])
                if isinstance(res_list, list) and res_list:
                    resolution = res_list[0]
                elif isinstance(res_list, (int, float)):
                    resolution = res_list
            except:
                pass
            
            # Polymer entities → chains + UniProt
            chains = []
            poly_entities = entry_data.get("polymer_entities") or []
            for entity in poly_entities:
                if not entity:
                    continue
                
                # Check protein type
                try:
                    etype = (entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type", "")
                    if etype and etype != "Protein":
                        continue
                except:
                    pass
                
                # UniProt
                uniprot = None
                try:
                    container = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                    uids = container.get("uniprot_ids") or []
                    if uids:
                        uniprot = uids[0]
                except:
                    pass
                
                # Chain IDs
                chain_ids = []
                try:
                    container = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                    chain_ids = container.get("auth_asym_ids") or []
                except:
                    pass
                
                # Sequence length
                seq_len = 0
                try:
                    seq = (entity.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can", "")
                    seq_len = len(seq) if seq else 0
                except:
                    pass
                
                if uniprot and chain_ids and seq_len >= 50:
                    chains.append({
                        "uniprot": uniprot,
                        "chain_ids": chain_ids,
                        "seq_len": seq_len,
                    })
            
            # Nonpolymer entities → ligand detection
            has_real_ligand = False
            ligand_ids = []
            nonpoly_entities = entry_data.get("nonpolymer_entities") or []
            for entity in nonpoly_entities:
                if not entity:
                    continue
                try:
                    container = entity.get("rcsb_nonpolymer_entity_container_identifiers") or {}
                    comp_id = container.get("comp_id", "")
                except:
                    continue
                if comp_id and comp_id.upper() not in SOLVENT_LIGANDS:
                    has_real_ligand = True
                    ligand_ids.append(comp_id)
            
            if chains:
                all_info[pdb_id] = {
                    "resolution": resolution,
                    "chains": chains,
                    "has_ligand": has_real_ligand,
                    "ligands": ligand_ids,
                }
        
        count = min(batch_start + batch_size, len(pdb_ids))
        if count % 200 == 0 or count == len(pdb_ids):
            print(f"    Processed {count}/{len(pdb_ids)} ({len(all_info)} valid)...")
        time.sleep(0.3)  # rate limit
    
    print(f"  Got metadata for {len(all_info)} structures")
    return all_info

def build_apo_holo_pairs(structure_info, max_pairs=300):
    """
    Group structures by UniProt ID and create apo-holo pairs.
    For each UniProt, pick best-resolution apo and best-resolution holo.
    """
    
    # Group by UniProt
    by_uniprot = defaultdict(lambda: {"apo": [], "holo": []})
    
    for pdb_id, info in structure_info.items():
        for chain_info in info["chains"]:
            uniprot = chain_info["uniprot"]
            entry = {
                "pdb": pdb_id,
                "chain": chain_info["chain_ids"][0],  # take first chain
                "resolution": info["resolution"] or 99.0,
                "seq_len": chain_info["seq_len"],
            }
            if info["has_ligand"]:
                by_uniprot[uniprot]["holo"].append(entry)
            else:
                by_uniprot[uniprot]["apo"].append(entry)
    
    # Build pairs: UniProts that have BOTH apo and holo
    pairs = []
    for uniprot, entries in by_uniprot.items():
        if not entries["apo"] or not entries["holo"]:
            continue
        
        # Pick best resolution for each
        best_apo = min(entries["apo"], key=lambda x: x["resolution"])
        best_holo = min(entries["holo"], key=lambda x: x["resolution"])
        
        # Skip if same PDB (shouldn't happen, but safety check)
        if best_apo["pdb"] == best_holo["pdb"]:
            continue
        
        pairs.append({
            "uniprot": uniprot,
            "apo_pdb": best_apo["pdb"],
            "apo_chain": best_apo["chain"],
            "apo_res": best_apo["resolution"],
            "holo_pdb": best_holo["pdb"],
            "holo_chain": best_holo["chain"],
            "holo_res": best_holo["resolution"],
            "seq_len": best_apo["seq_len"],
        })
    
    # Sort by average resolution (best first)
    pairs.sort(key=lambda x: (x["apo_res"] + x["holo_res"]) / 2)
    
    return pairs[:max_pairs]


# ═══════════════════════════════════════════════════════════════════
# PDB DOWNLOAD AND PARSING (same as previous scripts)
# ═══════════════════════════════════════════════════════════════════

def download_pdb(pdb_id, data_dir=PDB_DIR):
    """Download PDB file, trying .pdb.gz first for speed."""
    data_dir.mkdir(parents=True, exist_ok=True)
    filepath = data_dir / f"{pdb_id.upper()}.pdb"
    if filepath.exists():
        return filepath
    
    # Try gzipped first
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
    
    # Fallback to uncompressed
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    try:
        urlretrieve(url, filepath)
        return filepath
    except:
        return None


def parse_backbone(filepath, chain_id="A"):
    """Parse backbone N, CA, C atoms from PDB."""
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


# ═══════════════════════════════════════════════════════════════════
# CURVATURE COMPUTATION (identical to atlas methodology)
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
            prn = sr[i-1]
            pr = residues.get(prn, {})
            if "C" in pr and (rn - prn) == 1:
                phi = dihedral(pr["C"], r["N"], r["CA"], r["C"])
        psi = np.nan
        if i < len(sr) - 1:
            nrn = sr[i+1]
            nr = residues.get(nrn, {})
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
        dm = 0.5 * (a1 + a2)
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
    return (2 * np.sum(np.arange(1, n+1) * v) / (n * s)) - (n+1) / n


def compute_gini(filepath, chain_id):
    """Full pipeline: PDB → Gini coefficient."""
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
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(description="Large-scale Gini gradient test")
    parser.add_argument("--max-pairs", type=int, default=200,
                        help="Maximum number of apo-holo pairs (default: 200)")
    parser.add_argument("--max-search", type=int, default=3000,
                        help="Maximum PDB entries to search (default: 3000)")
    parser.add_argument("--resolution", type=float, default=2.0,
                        help="Resolution cutoff in Å (default: 2.0)")
    parser.add_argument("--skip-search", action="store_true",
                        help="Skip search, use cached pairs")
    args = parser.parse_args()
    
    DATA_DIR.mkdir(exist_ok=True)
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  LARGE-SCALE GRADIENT TEST                                     ║")
    print(f"║  Target: {args.max_pairs} apo-holo pairs, resolution ≤ {args.resolution} Å            ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    
    pairs_file = DATA_DIR / "apo_holo_pairs.json"
    
    if args.skip_search and pairs_file.exists():
        print("\n── Loading cached pairs ──")
        with open(pairs_file) as f:
            pairs = json.load(f)
        print(f"  Loaded {len(pairs)} pairs")
    else:
        # Step 1: Search PDB
        print("\n── Step 1: Search RCSB PDB ──")
        pdb_ids = find_apo_holo_pairs(
            max_results=args.max_search,
            resolution=args.resolution
        )
        
        if not pdb_ids:
            print("No structures found. Check network connection.")
            sys.exit(1)
        
        # Step 2: Get metadata
        print("\n── Step 2: Fetch structure metadata ──")
        info = get_structure_info_batch(pdb_ids)
        print(f"  Got metadata for {len(info)} structures")
        
        # Save raw metadata
        with open(DATA_DIR / "structure_metadata.json", "w") as f:
            json.dump(info, f, indent=2)
        
        # Step 3: Build pairs
        print("\n── Step 3: Build apo-holo pairs ──")
        pairs = build_apo_holo_pairs(info, max_pairs=args.max_pairs)
        print(f"  Built {len(pairs)} pairs")
        
        with open(pairs_file, "w") as f:
            json.dump(pairs, f, indent=2)
    
    # Step 4: Download and compute
    print(f"\n── Step 4: Download PDBs and compute Gini ({len(pairs)} pairs) ──")
    
    results = []
    failed = 0
    
    for idx, pair in enumerate(pairs):
        if idx % 20 == 0:
            print(f"  Processing pair {idx+1}/{len(pairs)}...")
        
        # Download apo
        apo_path = download_pdb(pair["apo_pdb"])
        if not apo_path:
            failed += 1
            continue
        
        # Download holo
        holo_path = download_pdb(pair["holo_pdb"])
        if not holo_path:
            failed += 1
            continue
        
        # Compute Gini
        apo_result = compute_gini(apo_path, pair["apo_chain"])
        holo_result = compute_gini(holo_path, pair["holo_chain"])
        
        if apo_result is None or holo_result is None:
            failed += 1
            continue
        
        dg = holo_result["gini"] - apo_result["gini"]
        
        results.append({
            "uniprot": pair["uniprot"],
            "apo_pdb": pair["apo_pdb"],
            "holo_pdb": pair["holo_pdb"],
            "g_apo": apo_result["gini"],
            "g_holo": holo_result["gini"],
            "dg": dg,
            "n_res_apo": apo_result["n_res"],
            "n_res_holo": holo_result["n_res"],
            "apo_res": pair.get("apo_res"),
            "holo_res": pair.get("holo_res"),
        })
        
        # Brief rate limiting for downloads
        time.sleep(0.1)
    
    print(f"\n  Completed: {len(results)} pairs ({failed} failed)")
    
    if len(results) < 10:
        print("Too few results for analysis.")
        sys.exit(1)
    
    # ═══════════════════════════════════════════════════════════════
    # ANALYSIS
    # ═══════════════════════════════════════════════════════════════
    
    print("\n" + "="*70)
    print("ANALYSIS")
    print("="*70)
    
    g_apos = np.array([r["g_apo"] for r in results])
    g_holos = np.array([r["g_holo"] for r in results])
    dgs = np.array([r["dg"] for r in results])
    
    print(f"\n  Dataset: {len(results)} apo-holo pairs")
    print(f"  Apo Gini:  mean={np.mean(g_apos):.4f} ± {np.std(g_apos):.4f}  "
          f"range=[{np.min(g_apos):.4f}, {np.max(g_apos):.4f}]")
    print(f"  Holo Gini: mean={np.mean(g_holos):.4f} ± {np.std(g_holos):.4f}  "
          f"range=[{np.min(g_holos):.4f}, {np.max(g_holos):.4f}]")
    print(f"  ΔGini:     mean={np.mean(dgs):+.4f} ± {np.std(dgs):.4f}")
    print(f"  % stiffening (ΔG > 0): {100*np.mean(dgs > 0):.1f}%")
    print(f"  % loosening  (ΔG < 0): {100*np.mean(dgs < 0):.1f}%")
    
    # Threshold scan
    print("\n── Threshold Scan ──")
    best_score = 0
    best_thresholds = []
    
    for ti in range(600, 950):
        th = ti / 1000
        correct = 0
        for r in results:
            pred = "+" if r["g_apo"] < th else "-"
            actual = "+" if r["dg"] > 0 else "-"
            if pred == actual:
                correct += 1
        acc = correct / len(results)
        if correct > best_score:
            best_score = correct
            best_thresholds = [th]
        elif correct == best_score:
            best_thresholds.append(th)
    
    best_acc = best_score / len(results)
    print(f"  Best accuracy: {best_score}/{len(results)} ({100*best_acc:.1f}%)")
    print(f"  Best threshold range: [{best_thresholds[0]:.3f}, {best_thresholds[-1]:.3f}]")
    print(f"  Midpoint: {(best_thresholds[0] + best_thresholds[-1])/2:.3f}")
    
    # Report at key thresholds
    for label, th in [("0.810 (original)", 0.810),
                       ("0.835 (AlphaFold)", 0.835),
                       ("0.860 (pilot)", 0.860),
                       ("Best midpoint", (best_thresholds[0]+best_thresholds[-1])/2)]:
        correct = sum(1 for r in results 
                      if (("+" if r["g_apo"] < th else "-") == ("+" if r["dg"] > 0 else "-")))
        print(f"  At {label:>25s}: {correct}/{len(results)} ({100*correct/len(results):.1f}%)")
    
    # Correlation analysis
    print("\n── Correlation: G_apo vs ΔGini ──")
    corr = np.corrcoef(g_apos, dgs)[0, 1]
    print(f"  Pearson r = {corr:.4f}")
    
    # Bin analysis
    print("\n── Binned Analysis ──")
    bins = [(0.5, 0.70), (0.70, 0.75), (0.75, 0.80), (0.80, 0.85), 
            (0.85, 0.90), (0.90, 0.95), (0.95, 1.0)]
    print(f"  {'Bin':>12s} {'N':>5s} {'mean ΔG':>8s} {'% stiff':>8s} {'% loose':>8s}")
    print("  " + "─"*45)
    for lo, hi in bins:
        mask = (g_apos >= lo) & (g_apos < hi)
        n = np.sum(mask)
        if n < 2:
            continue
        mean_dg = np.mean(dgs[mask])
        pct_stiff = 100 * np.mean(dgs[mask] > 0)
        pct_loose = 100 * np.mean(dgs[mask] < 0)
        print(f"  [{lo:.2f},{hi:.2f}) {n:>5d} {mean_dg:>+8.4f} {pct_stiff:>7.1f}% {pct_loose:>7.1f}%")
    
    # Effect size by distance from threshold
    print("\n── Effect Size by Distance from Threshold ──")
    mid = (best_thresholds[0] + best_thresholds[-1]) / 2
    dist = np.abs(g_apos - mid)
    for d_lo, d_hi, label in [(0, 0.02, "near (0-0.02)"), 
                                (0.02, 0.05, "mid (0.02-0.05)"),
                                (0.05, 0.10, "far (0.05-0.10)"),
                                (0.10, 0.50, "very far (0.10+)")]:
        mask = (dist >= d_lo) & (dist < d_hi)
        n = np.sum(mask)
        if n < 2:
            continue
        correct = sum(1 for i, r in enumerate(results) 
                      if mask[i] and (("+" if r["g_apo"] < mid else "-") == ("+" if r["dg"] > 0 else "-")))
        acc = correct / n
        print(f"  {label:>20s}: {correct}/{n} ({100*acc:.1f}%) mean|ΔG|={np.mean(np.abs(dgs[mask])):.4f}")
    
    # Export
    export = {
        "n_pairs": len(results),
        "best_accuracy": best_acc,
        "best_threshold_range": [best_thresholds[0], best_thresholds[-1]],
        "midpoint": (best_thresholds[0] + best_thresholds[-1]) / 2,
        "pearson_r": corr,
        "mean_g_apo": float(np.mean(g_apos)),
        "mean_g_holo": float(np.mean(g_holos)),
        "mean_dg": float(np.mean(dgs)),
        "pct_stiffening": float(100 * np.mean(dgs > 0)),
        "results": results,
    }
    
    with open(RESULTS_FILE, "w") as f:
        json.dump(export, f, indent=2)
    print(f"\nResults exported to: {RESULTS_FILE}")
    
    # Also export a simple CSV for plotting
    csv_path = DATA_DIR / "gradient_data.csv"
    with open(csv_path, "w") as f:
        f.write("uniprot,apo_pdb,holo_pdb,g_apo,g_holo,delta_gini,n_res_apo\n")
        for r in sorted(results, key=lambda x: x["g_apo"]):
            f.write(f"{r['uniprot']},{r['apo_pdb']},{r['holo_pdb']},"
                    f"{r['g_apo']:.6f},{r['g_holo']:.6f},{r['dg']:.6f},"
                    f"{r['n_res_apo']}\n")
    print(f"CSV exported to: {csv_path}")


if __name__ == "__main__":
    main()
