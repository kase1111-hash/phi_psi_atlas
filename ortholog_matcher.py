#!/usr/bin/env python3
"""
Ortholog Matcher — Find cross-species replacement candidates for disease singularities
======================================================================================

Takes the 222 human disease singularities and:
  1. Queries UniProt for ortholog/homolog IDs across species
  2. Pulls those orthologs from all_species.csv
  3. Compares mechanical fingerprints (Gini, concentration, hotspot position)
  4. Ranks replacement candidates by mechanical similarity

Usage:
  python ortholog_matcher.py \
    --diseases human_disease_singularities.csv \
    --species all_species.csv \
    --output ortholog_matches.csv \
    --workers 8

Requires network access for UniProt API.
Run time: ~10-30 minutes for 222 genes at 8 workers.
"""

import os
import sys
import csv
import json
import time
import argparse
import urllib.request
import urllib.error
from collections import defaultdict
import concurrent.futures
import threading


# ═══════════════════════════════════════════════════════════════
# UNIPROT ORTHOLOG LOOKUP
# ═══════════════════════════════════════════════════════════════

def fetch_orthologs_by_gene(gene_name, retries=3):
    """Search UniProt for orthologs of a human gene across species.
    Returns list of (uniprot_id, organism, gene_name) tuples."""
    
    # Search for reviewed entries with this gene name across all organisms
    query = f'(gene_exact:"{gene_name}") AND (reviewed:true)'
    url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&format=json&size=50&fields=accession,gene_names,organism_name,organism_id,sequence"
    
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read())
                results = []
                for entry in data.get("results", []):
                    uid = entry.get("primaryAccession", "")
                    org = entry.get("organism", {}).get("scientificName", "")
                    org_id = entry.get("organism", {}).get("taxonId", 0)
                    genes = entry.get("genes", [])
                    gname = genes[0].get("geneName", {}).get("value", "") if genes else ""
                    seq_len = entry.get("sequence", {}).get("length", 0)
                    results.append({
                        "uid": uid,
                        "organism": org,
                        "taxon_id": org_id,
                        "gene": gname,
                        "length": seq_len,
                    })
                return results
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return []
            if attempt < retries - 1:
                time.sleep(1.5 * (attempt + 1))
        except Exception:
            if attempt < retries - 1:
                time.sleep(1.5 * (attempt + 1))
    return []


def fetch_orthologs_by_uid(uniprot_id, retries=3):
    """Get orthologs via UniProt's cross-references (OrthoDB, OMA, etc).
    Falls back to gene name search if no direct ortholog data."""
    
    # First try: get the entry and check for ortholog cross-refs
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=10) as resp:
                entry = json.loads(resp.read())
                
                # Get gene name for fallback search
                genes = entry.get("genes", [])
                gene_name = genes[0].get("geneName", {}).get("value", "") if genes else ""
                
                # Try gene name search across species
                if gene_name:
                    return gene_name, fetch_orthologs_by_gene(gene_name)
                
                return "", []
        except:
            if attempt < retries - 1:
                time.sleep(1 * (attempt + 1))
    return "", []


# ═══════════════════════════════════════════════════════════════
# SPECIES CENSUS LOADER
# ═══════════════════════════════════════════════════════════════

# Map UniProt organism names to all_species.csv species tags
SPECIES_MAP = {
    "Homo sapiens": "human",
    "Mus musculus": "mouse",
    "Rattus norvegicus": "rat",
    "Bos taurus": "cow",
    "Sus scrofa": "pig",
    "Gallus gallus": "chicken",
    "Danio rerio": "zebrafish",
    "Drosophila melanogaster": "fly",
    "Caenorhabditis elegans": "worm",
    "Saccharomyces cerevisiae": "yeast",
    "Schizosaccharomyces pombe": "fission_yeast",
    "Arabidopsis thaliana": "arabidopsis",
    "Oryza sativa": "rice",
    "Escherichia coli": "ecoli",
    "Bacillus subtilis": "bacillus",
    "Mycobacterium tuberculosis": "tuberculosis",
    "Dictyostelium discoideum": "dictyostelium",
    "Candida albicans": "candida",
    "Salmonella typhimurium": "salmonella",
    "Pseudomonas aeruginosa": "pseudomonas",
}

def load_species_index(species_csv):
    """Load all_species.csv indexed by (species, uniprot_id)."""
    index = {}
    by_species_uid = {}
    
    with open(species_csv) as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                uid = r["uniprot_id"]
                species = r["species"]
                entry = {
                    "uid": uid,
                    "species": species,
                    "length": int(r["length"]),
                    "gini": float(r["kappa_gini"]),
                    "k_per_res": float(r["kappa_sq_per_res"]),
                    "k_max": float(r["top_peak_kappa_sq"]),
                    "hotspot_pos": r["top_peak_position"],
                    "f_alpha": float(r["ss_frac_H"]),
                    "f_beta": float(r["ss_frac_E"]),
                    "plddt": float(r.get("mean_plddt", 0)),
                }
                entry["concentration"] = entry["k_max"] / entry["k_per_res"] if entry["k_per_res"] > 0 else 0
                
                # Relative hotspot position (0-1 scale)
                try:
                    entry["hotspot_rel"] = int(entry["hotspot_pos"]) / entry["length"] if entry["length"] > 0 else 0
                except:
                    entry["hotspot_rel"] = 0
                
                by_species_uid[(species, uid)] = entry
                
            except:
                pass
    
    print(f"  Loaded {len(by_species_uid):,} entries from species census")
    return by_species_uid


# ═══════════════════════════════════════════════════════════════
# MECHANICAL SIMILARITY SCORING
# ═══════════════════════════════════════════════════════════════

def mechanical_similarity(human_entry, ortholog_entry):
    """Score how mechanically similar two proteins are.
    Returns score 0-100 and breakdown."""
    
    score = 0
    details = []
    
    # Gini similarity (0-30 points)
    gini_diff = abs(human_entry["gini"] - ortholog_entry["gini"])
    if gini_diff < 0.01:
        score += 30
        details.append(f"Gini: excellent (Δ={gini_diff:.4f})")
    elif gini_diff < 0.03:
        score += 25
        details.append(f"Gini: good (Δ={gini_diff:.4f})")
    elif gini_diff < 0.05:
        score += 15
        details.append(f"Gini: moderate (Δ={gini_diff:.4f})")
    elif gini_diff < 0.10:
        score += 5
        details.append(f"Gini: weak (Δ={gini_diff:.4f})")
    else:
        details.append(f"Gini: poor (Δ={gini_diff:.4f})")
    
    # Concentration ratio similarity (0-25 points)
    h_conc = human_entry.get("concentration", 0)
    o_conc = ortholog_entry.get("concentration", 0)
    if h_conc > 0 and o_conc > 0:
        ratio = min(h_conc, o_conc) / max(h_conc, o_conc)
        if ratio > 0.7:
            score += 25
            details.append(f"Concentration: excellent ({ratio:.2f})")
        elif ratio > 0.4:
            score += 15
            details.append(f"Concentration: good ({ratio:.2f})")
        elif ratio > 0.2:
            score += 8
            details.append(f"Concentration: moderate ({ratio:.2f})")
        else:
            details.append(f"Concentration: poor ({ratio:.2f})")
    
    # Size similarity (0-20 points)
    h_len = human_entry["length"]
    o_len = ortholog_entry["length"]
    if h_len > 0 and o_len > 0:
        size_ratio = min(h_len, o_len) / max(h_len, o_len)
        if size_ratio > 0.9:
            score += 20
            details.append(f"Size: excellent ({h_len} vs {o_len})")
        elif size_ratio > 0.8:
            score += 15
            details.append(f"Size: good ({h_len} vs {o_len})")
        elif size_ratio > 0.6:
            score += 8
            details.append(f"Size: moderate ({h_len} vs {o_len})")
        else:
            details.append(f"Size: poor ({h_len} vs {o_len})")
    
    # SS composition similarity (0-15 points)
    alpha_diff = abs(human_entry["f_alpha"] - ortholog_entry["f_alpha"])
    beta_diff = abs(human_entry["f_beta"] - ortholog_entry["f_beta"])
    ss_diff = alpha_diff + beta_diff
    if ss_diff < 0.05:
        score += 15
        details.append(f"SS: excellent (Δα={alpha_diff:.3f} Δβ={beta_diff:.3f})")
    elif ss_diff < 0.10:
        score += 10
        details.append(f"SS: good")
    elif ss_diff < 0.20:
        score += 5
        details.append(f"SS: moderate")
    else:
        details.append(f"SS: poor")
    
    # Relative hotspot position (0-10 points)
    h_rel = human_entry.get("hotspot_rel", 0)
    o_rel = ortholog_entry.get("hotspot_rel", 0)
    if h_rel > 0 and o_rel > 0:
        pos_diff = abs(h_rel - o_rel)
        if pos_diff < 0.05:
            score += 10
            details.append(f"Hotspot position: excellent (Δ={pos_diff:.3f})")
        elif pos_diff < 0.10:
            score += 7
            details.append(f"Hotspot position: good")
        elif pos_diff < 0.20:
            score += 3
            details.append(f"Hotspot position: moderate")
    
    return score, details


# ═══════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════

def match_orthologs(diseases_csv, species_csv, output_csv, workers=4):
    """Main pipeline: find cross-species matches for disease singularities."""
    
    # Load disease list
    diseases = []
    with open(diseases_csv) as f:
        for r in csv.DictReader(f):
            diseases.append(r)
    
    print(f"  Loaded {len(diseases)} disease singularities")
    
    # Load species census
    species_index = load_species_index(species_csv)
    
    # Build human entries from species census for comparison
    human_index = {uid: entry for (sp, uid), entry in species_index.items() if sp == "human"}
    print(f"  Human entries in census: {len(human_index):,}")
    
    # Thread-safe state
    lock = threading.Lock()
    counter = {"done": 0, "found": 0, "matched": 0}
    all_matches = []
    
    def _process_one(disease):
        gene = disease.get("gene", "")
        uid = disease.get("uid", "")
        
        if not gene:
            with lock:
                counter["done"] += 1
            return []
        
        # Get human entry from species census
        human_entry = human_index.get(uid)
        if not human_entry:
            # Try to find by scanning
            pass
        
        # Fetch orthologs from UniProt
        gene_searched, orthologs = fetch_orthologs_by_uid(uid)
        
        # For each ortholog, check if it's in our species census
        matches = []
        for orth in orthologs:
            # Skip the human entry itself
            if orth["uid"] == uid:
                continue
            
            # Map organism to species tag
            org = orth["organism"]
            species_tag = None
            for org_name, tag in SPECIES_MAP.items():
                if org_name.lower() in org.lower():
                    species_tag = tag
                    break
            
            if not species_tag:
                continue
            
            # Look up in species census
            census_entry = species_index.get((species_tag, orth["uid"]))
            if not census_entry:
                continue
            
            # Compute mechanical similarity
            if human_entry:
                sim_score, sim_details = mechanical_similarity(human_entry, census_entry)
            else:
                sim_score = -1
                sim_details = ["No human census entry for comparison"]
            
            match = {
                "human_gene": gene,
                "human_uid": uid,
                "human_disease": disease.get("disease_names", "")[:150],
                "human_gini": human_entry["gini"] if human_entry else "",
                "human_concentration": disease.get("concentration_ratio", ""),
                "human_hotspot": disease.get("hotspot_residue", ""),
                "ortholog_uid": orth["uid"],
                "ortholog_species": species_tag,
                "ortholog_organism": orth["organism"],
                "ortholog_gene": orth["gene"],
                "ortholog_length": census_entry["length"],
                "ortholog_gini": round(census_entry["gini"], 4),
                "ortholog_concentration": round(census_entry["concentration"], 1),
                "ortholog_hotspot": census_entry["hotspot_pos"],
                "ortholog_f_alpha": round(census_entry["f_alpha"], 3),
                "ortholog_f_beta": round(census_entry["f_beta"], 3),
                "similarity_score": sim_score,
                "similarity_details": "; ".join(sim_details),
            }
            matches.append(match)
        
        with lock:
            counter["done"] += 1
            if orthologs:
                counter["found"] += 1
            if matches:
                counter["matched"] += 1
                all_matches.extend(matches)
            
            done = counter["done"]
            total = len(diseases)
            n_orth = len(orthologs)
            n_match = len(matches)
            best = max((m["similarity_score"] for m in matches), default=0)
            
            if done % 10 == 0 or done == total:
                print(f"  [{done:>4d}/{total}] {gene:12s} orth={n_orth:>3d} match={n_match:>2d} best={best:>3d}")
        
        time.sleep(0.3)  # Rate limit UniProt API
        return matches
    
    # Execute
    print(f"\n  Searching for orthologs with {workers} workers...")
    t0 = time.time()
    
    if workers <= 1:
        for d in diseases:
            _process_one(d)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_process_one, d) for d in diseases]
            concurrent.futures.wait(futures)
    
    elapsed = time.time() - t0
    
    # Sort by similarity score
    all_matches.sort(key=lambda m: (-m["similarity_score"], m["human_gene"]))
    
    # Write output
    if output_csv and all_matches:
        fields = list(all_matches[0].keys())
        with open(output_csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for m in all_matches:
                w.writerow(m)
        print(f"\n  Wrote {len(all_matches)} matches to {output_csv}")
    
    # Summary
    print(f"\n{'='*70}")
    print(f"  ORTHOLOG MATCHING SUMMARY")
    print(f"{'='*70}")
    print(f"  Disease genes searched:   {len(diseases)}")
    print(f"  Found orthologs:          {counter['found']}")
    print(f"  Matched in census:        {counter['matched']}")
    print(f"  Total matches:            {len(all_matches)}")
    print(f"  Time: {elapsed:.0f}s")
    
    # Best matches per gene
    if all_matches:
        print(f"\n  TOP MATCHES (similarity ≥ 70):")
        print(f"  {'Gene':12s} {'Species':15s} {'Score':>5s} {'H_Gini':>7s} {'O_Gini':>7s} {'H_Conc':>7s} {'O_Conc':>7s}")
        print(f"  {'─'*70}")
        
        seen = set()
        for m in all_matches:
            if m["similarity_score"] >= 70:
                key = (m["human_gene"], m["ortholog_species"])
                if key not in seen:
                    seen.add(key)
                    print(f"  {m['human_gene']:12s} {m['ortholog_species']:15s} {m['similarity_score']:>5d} "
                          f"{m['human_gini']:>7s} {m['ortholog_gini']:>7.4f} "
                          f"{m['human_concentration']:>7s} {m['ortholog_concentration']:>7.1f}")
        
        # Count species coverage
        species_counts = defaultdict(int)
        gene_coverage = defaultdict(set)
        for m in all_matches:
            species_counts[m["ortholog_species"]] += 1
            gene_coverage[m["human_gene"]].add(m["ortholog_species"])
        
        print(f"\n  Species coverage:")
        for sp, n in sorted(species_counts.items(), key=lambda x: -x[1]):
            print(f"    {sp:20s} {n:>5d} matches")
        
        print(f"\n  Genes with orthologs in 5+ species: "
              f"{sum(1 for g,spp in gene_coverage.items() if len(spp)>=5)}")
        print(f"  Genes with orthologs in 10+ species: "
              f"{sum(1 for g,spp in gene_coverage.items() if len(spp)>=10)}")
    
    return all_matches


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find cross-species ortholog matches for disease singularities")
    parser.add_argument("--diseases", required=True,
                        help="human_disease_singularities.csv")
    parser.add_argument("--species", required=True,
                        help="all_species.csv")
    parser.add_argument("--output", default="ortholog_matches.csv",
                        help="Output CSV")
    parser.add_argument("--workers", type=int, default=4,
                        help="Parallel threads (default: 4, be gentle on UniProt API)")
    
    args = parser.parse_args()
    
    match_orthologs(
        args.diseases,
        species_csv=args.species,
        output_csv=args.output,
        workers=args.workers,
    )
