#!/usr/bin/env python3
"""
Singularity Annotator — Map mechanical singularities to genes, diseases, and orthologs
======================================================================================

Takes the mechanical_singularities.csv and:
  1. Queries UniProt API for gene names, protein names, disease associations
  2. Cross-references hotspot residues with known pathogenic variants
  3. Finds orthologs in the all_species.csv census to check hotspot conservation
  4. Outputs a ranked therapeutic priority list

Usage:
  python singularity_annotator.py \
    --singularities mechanical_singularities.csv \
    --species all_species.csv \
    --output annotated_singularities.csv \
    --workers 8

Requires network access for UniProt API calls.
"""

import os
import sys
import csv
import json
import time
import argparse
import urllib.request
import urllib.error
from pathlib import Path
from collections import defaultdict
import concurrent.futures
import threading


# ═══════════════════════════════════════════════════════════════
# UNIPROT API
# ═══════════════════════════════════════════════════════════════

def fetch_uniprot_entry(uniprot_id, retries=3):
    """Fetch protein annotation from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                return json.loads(resp.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None  # Not found
            if attempt < retries - 1:
                time.sleep(1 * (attempt + 1))
        except Exception:
            if attempt < retries - 1:
                time.sleep(1 * (attempt + 1))
    return None


def parse_uniprot_entry(entry):
    """Extract relevant fields from UniProt JSON."""
    if not entry:
        return {}
    
    result = {}
    
    # Gene name
    genes = entry.get("genes", [])
    if genes:
        primary = genes[0]
        result["gene"] = primary.get("geneName", {}).get("value", "")
        synonyms = [s.get("value","") for s in primary.get("synonyms", [])]
        result["gene_synonyms"] = "; ".join(synonyms) if synonyms else ""
    
    # Protein name
    desc = entry.get("proteinDescription", {})
    rec = desc.get("recommendedName", {})
    if rec:
        result["protein_name"] = rec.get("fullName", {}).get("value", "")
    else:
        sub = desc.get("submissionNames", [])
        if sub:
            result["protein_name"] = sub[0].get("fullName", {}).get("value", "")
    
    # Organism
    org = entry.get("organism", {})
    result["organism"] = org.get("scientificName", "")
    result["common_name"] = org.get("commonName", "")
    result["taxon_id"] = org.get("taxonId", "")
    
    # Function
    comments = entry.get("comments", [])
    for c in comments:
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            if texts:
                result["function"] = texts[0].get("value", "")[:500]
    
    # Disease associations
    diseases = []
    for c in comments:
        if c.get("commentType") == "DISEASE":
            disease = c.get("disease", {})
            if disease:
                d = {
                    "name": disease.get("diseaseId", ""),
                    "description": disease.get("description", "")[:200],
                    "acronym": disease.get("acronym", ""),
                }
                ref = disease.get("diseaseCrossReference", {})
                if ref:
                    d["mim_id"] = ref.get("id", "")
                diseases.append(d)
    result["diseases"] = diseases
    result["disease_count"] = len(diseases)
    result["disease_names"] = "; ".join(d["name"] for d in diseases) if diseases else ""
    
    # Subcellular location
    for c in comments:
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            locs = []
            for sl in c.get("subcellularLocations", []):
                loc = sl.get("location", {})
                if loc:
                    locs.append(loc.get("value", ""))
            result["location"] = "; ".join(locs) if locs else ""
    
    # Known variants (features)
    variants = []
    features = entry.get("features", [])
    for f in features:
        if f.get("type") == "VARIANT":
            loc = f.get("location", {})
            start = loc.get("start", {}).get("value")
            end = loc.get("end", {}).get("value")
            desc_list = f.get("description", "")
            fev = f.get("featureId", "")
            variants.append({
                "position": start,
                "description": desc_list,
                "id": fev,
            })
    result["variants"] = variants
    result["variant_count"] = len(variants)
    
    # Pathway / keywords
    keywords = entry.get("keywords", [])
    result["keywords"] = "; ".join(kw.get("name","") for kw in keywords[:10])
    
    # Sequence length
    seq = entry.get("sequence", {})
    result["seq_length"] = seq.get("length", 0)
    
    return result


def check_variant_at_hotspot(annotation, hotspot_res, window=5):
    """Check if any known pathogenic variant is at or near the hotspot residue."""
    hotspot = int(hotspot_res) if hotspot_res else 0
    if hotspot == 0:
        return []
    
    hits = []
    for v in annotation.get("variants", []):
        pos = v.get("position")
        if pos is None:
            continue
        try:
            pos = int(pos)
        except:
            continue
        if abs(pos - hotspot) <= window:
            hits.append({
                "position": pos,
                "distance": abs(pos - hotspot),
                "description": v.get("description", ""),
                "id": v.get("id", ""),
            })
    return hits


# ═══════════════════════════════════════════════════════════════
# ORTHOLOG FINDER — using all_species census
# ═══════════════════════════════════════════════════════════════

def load_species_census(species_csv):
    """Load the multi-species census for ortholog comparison."""
    by_species = defaultdict(list)
    
    with open(species_csv) as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                by_species[r["species"]].append({
                    "uid": r["uniprot_id"],
                    "species": r["species"],
                    "length": int(r["length"]),
                    "gini": float(r["kappa_gini"]),
                    "k_mean": float(r["kappa_sq_per_res"]),
                    "k_max": float(r["top_peak_kappa_sq"]),
                    "hotspot_pos": r["top_peak_position"],
                    "f_alpha": float(r["ss_frac_H"]),
                    "f_beta": float(r["ss_frac_E"]),
                })
            except:
                pass
    
    return by_species


def find_mechanical_matches(target, species_census, gini_tol=0.03, size_tol=0.2, k_ratio_tol=0.5):
    """Find proteins in other species with similar mechanical fingerprints.
    
    This is NOT ortholog mapping (that needs sequence alignment).
    This finds proteins with matching mechanical properties — which could
    be orthologs OR convergent mechanical solutions.
    """
    matches = []
    t_gini = float(target["gini"])
    t_size = int(target["n"])
    t_k = float(target["k_mean"])
    t_conc = float(target["concentration_ratio"])
    
    for species, proteins in species_census.items():
        for p in proteins:
            # Size within tolerance
            if abs(p["length"] - t_size) / max(t_size, 1) > size_tol:
                continue
            # Gini within tolerance
            if abs(p["gini"] - t_gini) > gini_tol:
                continue
            # Concentration ratio similar
            p_conc = p["k_max"] / p["k_mean"] if p["k_mean"] > 0 else 0
            if t_conc > 0 and abs(p_conc - t_conc) / t_conc > k_ratio_tol:
                continue
            
            matches.append({
                "uid": p["uid"],
                "species": p["species"],
                "length": p["length"],
                "gini": p["gini"],
                "concentration_ratio": round(p_conc, 0),
                "hotspot_pos": p["hotspot_pos"],
            })
    
    return matches


# ═══════════════════════════════════════════════════════════════
# BATCH ANNOTATOR
# ═══════════════════════════════════════════════════════════════

def annotate_singularities(singularities_csv, species_csv=None, output_csv=None,
                           workers=4, max_entries=None):
    """Main pipeline: annotate all singularities with UniProt data."""
    
    # Load singularities
    targets = []
    with open(singularities_csv) as f:
        for r in csv.DictReader(f):
            targets.append(r)
    
    if max_entries:
        targets = targets[:max_entries]
    
    print(f"  Annotating {len(targets)} singularities with {workers} workers...")
    
    # Load species census if available
    species_census = None
    if species_csv and os.path.exists(species_csv):
        print(f"  Loading species census for ortholog matching...")
        species_census = load_species_census(species_csv)
        total_sp = sum(len(v) for v in species_census.values())
        print(f"  Loaded {total_sp:,} proteins from {len(species_census)} species")
    
    # Thread-safe state
    lock = threading.Lock()
    counter = {"done": 0, "ok": 0, "disease": 0, "hotspot_hit": 0}
    results = []
    
    def _annotate_one(target):
        uid = target["uid"]
        hotspot = target["top1_res"]
        conc = float(target["concentration_ratio"])
        
        # Fetch UniProt
        entry = fetch_uniprot_entry(uid)
        ann = parse_uniprot_entry(entry)
        
        # Check for variants near hotspot
        hotspot_variants = check_variant_at_hotspot(ann, hotspot, window=5)
        
        # Find mechanical matches in other species
        ortho_matches = []
        if species_census:
            ortho_matches = find_mechanical_matches(target, species_census)
        
        row = {
            "uid": uid,
            "gene": ann.get("gene", ""),
            "protein_name": ann.get("protein_name", ""),
            "organism": ann.get("organism", ""),
            "common_name": ann.get("common_name", ""),
            "function": ann.get("function", "")[:300],
            "location": ann.get("location", ""),
            "gini": target["gini"],
            "concentration_ratio": target["concentration_ratio"],
            "hotspot_residue": hotspot,
            "n_residues": target["n"],
            "disease_count": ann.get("disease_count", 0),
            "disease_names": ann.get("disease_names", ""),
            "variant_count": ann.get("variant_count", 0),
            "hotspot_variant_count": len(hotspot_variants),
            "hotspot_variants": "; ".join(
                f"pos{v['position']}(d={v['distance']}): {v['description'][:80]}"
                for v in hotspot_variants
            ),
            "keywords": ann.get("keywords", ""),
            "n_species_matches": len(ortho_matches),
            "species_with_matches": "; ".join(
                sorted(set(m["species"] for m in ortho_matches))
            ),
            "therapeutic_priority": "",  # Will be computed later
        }
        
        with lock:
            counter["done"] += 1
            if ann.get("gene"):
                counter["ok"] += 1
            if ann.get("disease_count", 0) > 0:
                counter["disease"] += 1
            if hotspot_variants:
                counter["hotspot_hit"] += 1
            
            results.append(row)
            
            done = counter["done"]
            total = len(targets)
            if done % 20 == 0 or done == total:
                gene = ann.get("gene", "?")
                dis = f"DIS={ann.get('disease_count',0)}" if ann.get("disease_count",0) else ""
                hv = f"HV={len(hotspot_variants)}" if hotspot_variants else ""
                print(f"  [{done:>4d}/{total}] {uid:15s} {gene:12s} {dis:>6s} {hv}")
        
        return row
    
    # Execute
    t0 = time.time()
    
    if workers <= 1:
        for t in targets:
            _annotate_one(t)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_annotate_one, t) for t in targets]
            concurrent.futures.wait(futures)
    
    elapsed = time.time() - t0
    
    # Compute therapeutic priority score
    for r in results:
        score = 0
        # Has known disease association
        if int(r["disease_count"]) > 0:
            score += 30
        # Has variant at hotspot
        if int(r["hotspot_variant_count"]) > 0:
            score += 40
        # High concentration ratio (more mechanically critical)
        conc = float(r["concentration_ratio"])
        if conc > 500:
            score += 20
        elif conc > 200:
            score += 10
        # Has cross-species matches (replacement candidates exist)
        if int(r["n_species_matches"]) > 0:
            score += 10
        
        r["therapeutic_priority"] = score
    
    # Sort by priority
    results.sort(key=lambda r: -int(r["therapeutic_priority"]))
    
    # Write output
    if output_csv:
        fields = list(results[0].keys()) if results else []
        with open(output_csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in results:
                w.writerow(r)
        print(f"\n  Wrote {len(results)} rows to {output_csv}")
    
    # Summary
    print(f"\n{'='*70}")
    print(f"  ANNOTATION SUMMARY")
    print(f"{'='*70}")
    print(f"  Total annotated:     {counter['ok']}/{len(targets)}")
    print(f"  With disease assoc:  {counter['disease']}")
    print(f"  Variant at hotspot:  {counter['hotspot_hit']}")
    print(f"  Time: {elapsed:.0f}s ({len(targets)/elapsed:.1f} proteins/sec)")
    
    # Top hits
    disease_hits = [r for r in results if int(r["disease_count"]) > 0]
    hotspot_hits = [r for r in results if int(r["hotspot_variant_count"]) > 0]
    
    if disease_hits:
        print(f"\n  {'─'*70}")
        print(f"  TOP DISEASE-ASSOCIATED SINGULARITIES:")
        print(f"  {'─'*70}")
        for r in disease_hits[:30]:
            conc = float(r["concentration_ratio"])
            hv = int(r["hotspot_variant_count"])
            hv_str = f" *** {hv} VARIANT(S) AT HOTSPOT ***" if hv > 0 else ""
            print(f"    {r['uid']:15s} {r['gene']:12s} {conc:>5.0f}× res={r['hotspot_residue']:>5s}  {r['disease_names'][:60]}{hv_str}")
    
    if hotspot_hits:
        print(f"\n  {'─'*70}")
        print(f"  CRITICAL: KNOWN VARIANTS AT MECHANICAL HOTSPOT:")
        print(f"  {'─'*70}")
        for r in hotspot_hits[:20]:
            print(f"    {r['uid']:15s} {r['gene']:12s} hotspot={r['hotspot_residue']}")
            print(f"      Diseases: {r['disease_names'][:80]}")
            print(f"      Variants: {r['hotspot_variants'][:120]}")
            print()
    
    return results


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate mechanical singularities")
    parser.add_argument("--singularities", required=True,
                        help="mechanical_singularities.csv from archetype_mapper")
    parser.add_argument("--species", default=None,
                        help="all_species.csv for cross-species matching")
    parser.add_argument("--output", default="annotated_singularities.csv",
                        help="Output CSV")
    parser.add_argument("--workers", type=int, default=4,
                        help="Parallel threads for API calls")
    parser.add_argument("--max", type=int, default=None,
                        help="Max entries to annotate (for testing)")
    
    args = parser.parse_args()
    
    annotate_singularities(
        args.singularities,
        species_csv=args.species,
        output_csv=args.output,
        workers=args.workers,
        max_entries=args.max,
    )
