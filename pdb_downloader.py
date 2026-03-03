#!/usr/bin/env python3
"""
PDB & AlphaFold Structure Downloader
=====================================
Downloads crystal structures from RCSB PDB and predicted structures from
AlphaFold DB for cross-species mechanical archetype mapping.

Usage:
  python pdb_downloader.py --mode pdb --ids 4HHB 1HHO 1AKE 4AKE
  python pdb_downloader.py --mode pdb --file pdb_list.txt
  python pdb_downloader.py --mode alphafold --uniprot P69905 P68871
  python pdb_downloader.py --mode alphafold --file uniprot_list.txt
  python pdb_downloader.py --mode family --query "adenylate kinase" --max 50
  python pdb_downloader.py --mode archetype-set
  python pdb_downloader.py --mode species-survey --protein "hemoglobin"

Output: CIF files in ./structures/ (or specified --outdir)
"""

import os
import sys
import json
import time
import argparse
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path


# ═══════════════════════════════════════════════════════════════
# DOWNLOAD HELPERS
# ═══════════════════════════════════════════════════════════════

def download_file(url, outpath, retries=3, delay=0.5):
    """Download a file with retry logic."""
    for attempt in range(retries):
        try:
            urllib.request.urlretrieve(url, outpath)
            size = os.path.getsize(outpath)
            if size > 100:  # Sanity check
                return True
            else:
                os.remove(outpath)
                return False
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return False
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
    return False


def download_pdb_cif(pdb_id, outdir):
    """Download mmCIF from RCSB PDB."""
    pdb_id = pdb_id.strip().upper()
    outpath = os.path.join(outdir, f"{pdb_id}.cif")
    if os.path.exists(outpath):
        print(f"  {pdb_id}: exists, skipping")
        return outpath
    
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    if download_file(url, outpath):
        size_kb = os.path.getsize(outpath) / 1024
        print(f"  {pdb_id}: OK ({size_kb:.0f} KB)")
        return outpath
    else:
        print(f"  {pdb_id}: FAILED")
        return None


def download_alphafold_cif(uniprot_id, outdir):
    """Download predicted structure from AlphaFold DB."""
    uniprot_id = uniprot_id.strip().upper()
    outpath = os.path.join(outdir, f"AF-{uniprot_id}-F1-model_v4.cif")
    if os.path.exists(outpath):
        print(f"  AF-{uniprot_id}: exists, skipping")
        return outpath
    
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"
    if download_file(url, outpath):
        size_kb = os.path.getsize(outpath) / 1024
        print(f"  AF-{uniprot_id}: OK ({size_kb:.0f} KB)")
        return outpath
    
    # Try v3 as fallback
    url_v3 = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v3.cif"
    outpath_v3 = os.path.join(outdir, f"AF-{uniprot_id}-F1-model_v3.cif")
    if download_file(url_v3, outpath_v3):
        size_kb = os.path.getsize(outpath_v3) / 1024
        print(f"  AF-{uniprot_id}: OK v3 ({size_kb:.0f} KB)")
        return outpath_v3
    
    print(f"  AF-{uniprot_id}: FAILED")
    return None


# ═══════════════════════════════════════════════════════════════
# RCSB SEARCH API — Find PDBs by protein family
# ═══════════════════════════════════════════════════════════════

def search_rcsb(query, max_results=50):
    """Search RCSB PDB using their search API."""
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": query
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "exact_match",
                        "value": "X-RAY DIFFRACTION"
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
            "paginate": {"start": 0, "rows": max_results}
        }
    }
    
    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(
        search_url,
        data=data,
        headers={"Content-Type": "application/json"}
    )
    
    try:
        with urllib.request.urlopen(req) as response:
            result = json.loads(response.read())
            ids = [r["identifier"] for r in result.get("result_set", [])]
            print(f"  Found {len(ids)} structures for '{query}'")
            return ids
    except Exception as e:
        print(f"  Search failed: {e}")
        return []


# ═══════════════════════════════════════════════════════════════
# CURATED ARCHETYPE SETS
# ═══════════════════════════════════════════════════════════════

# Cross-species pairs for the mechanical archetype question
ARCHETYPE_PAIRS = {
    # ═══ HINGE RELAY (AdK-like) ═══
    "AdK_Ecoli": {
        "apo": "4AKE", "holo": "1AKE", "chain": "A",
        "organism": "E. coli", "family": "AdK", "archetype": "HINGE_RELAY"
    },
    "AdK_human": {
        "apo": "2C95", "holo": "1Z83", "chain": "A",
        "organism": "Human", "family": "AdK", "archetype": "HINGE_RELAY"
    },
    "AdK_thermophile": {
        "apo": "1ZIN", "holo": "1ZIO", "chain": "A",
        "organism": "A. aeolicus", "family": "AdK", "archetype": "HINGE_RELAY"
    },
    
    # ═══ COOPERATIVE SWITCH (Hemoglobin) ═══
    "Hb_human_T": {
        "apo": "4HHB", "holo": "1HHO", "chain": "A",
        "organism": "Human", "family": "Hemoglobin", "archetype": "COOPERATIVE_SWITCH"
    },
    "Hb_horse": {
        "apo": "2HHB", "holo": "1G09", "chain": "A",
        "organism": "Horse", "family": "Hemoglobin", "archetype": "COOPERATIVE_SWITCH"
    },
    # Lamprey Hb is monomeric — should NOT show cooperativity
    "Hb_lamprey": {
        "apo": "2LHB", "holo": None, "chain": "A",
        "organism": "Lamprey", "family": "Hemoglobin", "archetype": "NO_COOPERATIVITY?"
    },
    # Clam hemoglobin — different quaternary structure
    "Hb_clam": {
        "apo": "1MOH", "holo": None, "chain": "A",
        "organism": "Clam", "family": "Hemoglobin", "archetype": "COOPERATIVE_SWITCH?"
    },
    
    # ═══ ION STIFFENER (Calmodulin-like) ═══
    "CaM_vertebrate_apo": {
        "apo": "1CFD", "holo": "1CLL", "chain": "A",
        "organism": "Vertebrate", "family": "Calmodulin", "archetype": "ION_STIFFENER"
    },
    "CaM_yeast": {
        "apo": "1F54", "holo": "1F55", "chain": "A",
        "organism": "Yeast", "family": "Calmodulin", "archetype": "ION_STIFFENER"
    },
    "TnC_human": {
        "apo": "1SPY", "holo": "1AP4", "chain": "A",
        "organism": "Human", "family": "TroponinC", "archetype": "ION_STIFFENER"
    },
    
    # ═══ ELASTIC SPRING (GCK-like) ═══
    "GCK_human": {
        "apo": "1V4S", "holo": "1V4T", "chain": "A",
        "organism": "Human", "family": "Glucokinase", "archetype": "ELASTIC_SPRING"
    },
    "HK_yeast": {
        "apo": "2YHX", "holo": "1IG8", "chain": "A",
        "organism": "Yeast", "family": "Hexokinase", "archetype": "ELASTIC_SPRING?"
    },
    
    # ═══ LOCKED GORGE (AChE-like) ═══
    "AChE_torpedo": {
        "apo": "1EA5", "holo": "1EVE", "chain": "A",
        "organism": "Torpedo", "family": "AChE", "archetype": "LOCKED_GORGE"
    },
    "AChE_human": {
        "apo": "4EY4", "holo": "4EY7", "chain": "A",
        "organism": "Human", "family": "AChE", "archetype": "LOCKED_GORGE"
    },
    "BChE_human": {
        "apo": "1P0I", "holo": "1P0M", "chain": "A",
        "organism": "Human", "family": "BChE", "archetype": "LOCKED_GORGE?"
    },
    
    # ═══ REDISTRIBUTOR (Myosin-like) ═══
    "Myosin2_dicty": {
        "apo": "1FMV", "holo": "1VOM", "chain": "A",
        "organism": "Dictyostelium", "family": "Myosin", "archetype": "REDISTRIBUTOR"
    },
    "MyosinV_chicken": {
        "apo": "1W8J", "holo": "1W7J", "chain": "A",
        "organism": "Chicken", "family": "MyosinV", "archetype": "REDISTRIBUTOR"
    },
    "Kinesin_human": {
        "apo": "1BG2", "holo": "3KIN", "chain": "A",
        "organism": "Human", "family": "Kinesin", "archetype": "REDISTRIBUTOR?"
    },
    
    # ═══ DOMAIN CLOSURE (Citrate Synthase) ═══
    "CS_pig_open": {
        "apo": "2CTS", "holo": "4CTS", "chain": "A",
        "organism": "Pig", "family": "CitrateSynthase", "archetype": "DOMAIN_CLOSURE"
    },
    "CS_thermo": {
        "apo": "1IOM", "holo": "1IXE", "chain": "A",
        "organism": "T. thermophilus", "family": "CitrateSynthase", "archetype": "DOMAIN_CLOSURE"
    },
    
    # ═══ ALLOSTERIC RELAY (CDK/Kinase) ═══
    "CDK2_human_flavo": {
        "apo": "1HCK", "holo": "1FTJ", "chain": "A",
        "organism": "Human", "family": "CDK2", "archetype": "ALLOSTERIC_RELAY"
    },
    "ERK2_human": {
        "apo": "1ERK", "holo": "2ERK", "chain": "A",
        "organism": "Human", "family": "ERK2", "archetype": "ALLOSTERIC_RELAY"
    },
    
    # ═══ PROTEASE (HIV-PR) ═══
    "HIVPR_subB": {
        "apo": "1HHP", "holo": "1HSG", "chain": "A",
        "organism": "HIV-1 subB", "family": "HIV_PR", "archetype": "MODERATE"
    },
}

# All unique PDB IDs needed
def get_all_pdb_ids():
    """Extract all PDB IDs from archetype pairs."""
    ids = set()
    for pair in ARCHETYPE_PAIRS.values():
        if pair.get("apo"):
            ids.add(pair["apo"])
        if pair.get("holo"):
            ids.add(pair["holo"])
    return sorted(ids)


# ═══ Species survey sets ═══
SPECIES_SURVEYS = {
    "hemoglobin": {
        "description": "Hemoglobin/myoglobin across species — cooperativity evolution",
        "pdbs": [
            # Human
            "4HHB", "1HHO", "2DN1", "2DN2",
            # Horse
            "2HHB", "1G09",
            # Cow
            "1FSX",
            # Lamprey (monomeric!)
            "2LHB",
            # Clam
            "1MOH",
            # Worm (C. elegans globin)
            "1D8U",
            # Myoglobin (sperm whale) — non-cooperative control
            "1MBD", "1BZR", "1BZP",
            # Leghemoglobin (plant!)
            "1GDI", "1BIN",
            # Truncated hemoglobin (bacterial)
            "1DLY", "1IDR",
            # Antarctic fish (no Hb!)
            # Fetal hemoglobin
            "1FDH",
        ],
        "alphafold_uniprots": [
            "P69905",   # Human HBA
            "P68871",   # Human HBB
            "P02144",   # Human myoglobin
            "P01942",   # Horse HBA
            "P02062",   # Horse HBB
            "P02185",   # Horse myoglobin
            "P80946",   # Lamprey Hb
            "P02189",   # Sperm whale Mb
        ]
    },
    "calmodulin": {
        "description": "EF-hand proteins across species — ion stiffening universality",
        "pdbs": [
            "1CFD", "1CLL",     # Vertebrate CaM apo/holo
            "1F54", "1F55",     # Yeast CaM apo/holo
            "1SPY", "1AP4",     # Troponin C
            "2BCX",             # Parvalbumin
            "1RJV",             # Recoverin
            "3CLN",             # Calmodulin 2.2A
        ],
        "alphafold_uniprots": [
            "P0DP23",   # Human CaM
            "P06787",   # Yeast CaM
            "P63098",   # Drosophila CaM
            "P62158",   # Bovine CaM
            "P02585",   # Human TnC slow
        ]
    },
    "lysozyme": {
        "description": "Lysozyme across species — fold conservation vs sequence divergence",
        "pdbs": [
            "1AKI", "1LYZ",    # Hen egg-white
            "6LYZ", "193L",    # Hen variants
            "1HEL",            # Hen high-res
            "1LZA",            # Human
            "2LZM",            # T4 phage (completely different fold!)
            "1GD6",            # Goose
            "1IWT",            # Turkey
            "1JWR",            # Cat
        ],
        "alphafold_uniprots": [
            "P61626",   # Human lysozyme
            "P00698",   # Hen lysozyme
            "P00699",   # Turkey lysozyme
            "P00703",   # Goose lysozyme
        ]
    },
    "kinase": {
        "description": "Protein kinases — allosteric relay conservation",
        "pdbs": [
            # CDK2
            "1HCK", "1FTJ", "1FTK",
            # ERK2
            "1ERK", "2ERK",
            # p38
            "1P38", "3GI3",
            # Abl
            "1IEP", "2HYY",
            # Src
            "2SRC", "1Y57",
            # PKA
            "1ATP", "1CDK",
            # Aurora
            "1MUO", "2J4Z",
        ],
        "alphafold_uniprots": [
            "P24941",   # CDK2
            "P28482",   # ERK2
            "Q16539",   # p38α
            "P00519",   # Abl1
            "P12931",   # Src
        ]
    },
    "protease": {
        "description": "Proteases across classes — mechanical diversity",
        "pdbs": [
            # HIV protease (aspartyl)
            "1HHP", "1HSG",
            # Trypsin (serine)
            "1S0Q", "1UTN",
            # Thrombin (serine)
            "1PPB", "1UDT",
            # Caspase-3 (cysteine)
            "1NME", "1PAU",
            # MMP-2 (metalloprotease)
            "1CK7",
            # Pepsin (aspartyl)
            "4PEP",
        ],
        "alphafold_uniprots": [
            "P07339",   # Cathepsin D
            "P07858",   # Cathepsin B
        ]
    },
    "motor": {
        "description": "Molecular motors — redistributor archetype across motor types",
        "pdbs": [
            # Myosin II
            "1FMV", "1VOM",
            # Myosin V
            "1W8J", "1W7J",
            # Kinesin
            "1BG2", "3KIN",
            # Dynein (too big for full analysis?)
            # F1-ATPase
            "1BMF",
            # GroEL (chaperonin motor)
            "1XCK",
        ],
        "alphafold_uniprots": [
            "P60709",   # Human beta-actin
            "P35579",   # Human myosin heavy chain
            "P33176",   # Human kinesin
        ]
    },
}


# ═══════════════════════════════════════════════════════════════
# BATCH DOWNLOAD FUNCTIONS
# ═══════════════════════════════════════════════════════════════

def download_pdb_list(pdb_ids, outdir):
    """Download a list of PDB structures."""
    os.makedirs(outdir, exist_ok=True)
    success = 0
    for pdb_id in pdb_ids:
        result = download_pdb_cif(pdb_id, outdir)
        if result:
            success += 1
        time.sleep(0.3)  # Be polite to RCSB
    print(f"\n  Downloaded {success}/{len(pdb_ids)} structures to {outdir}")
    return success


def download_alphafold_list(uniprot_ids, outdir):
    """Download a list of AlphaFold structures."""
    os.makedirs(outdir, exist_ok=True)
    success = 0
    for uid in uniprot_ids:
        result = download_alphafold_cif(uid, outdir)
        if result:
            success += 1
        time.sleep(0.3)
    print(f"\n  Downloaded {success}/{len(uniprot_ids)} AlphaFold models to {outdir}")
    return success


def download_archetype_set(outdir):
    """Download all PDBs needed for cross-species archetype mapping."""
    pdb_ids = get_all_pdb_ids()
    print(f"  Downloading {len(pdb_ids)} PDB structures for archetype mapping...")
    return download_pdb_list(pdb_ids, outdir)


def download_species_survey(survey_name, outdir):
    """Download all structures for a species survey."""
    if survey_name not in SPECIES_SURVEYS:
        print(f"  Unknown survey: {survey_name}")
        print(f"  Available: {', '.join(SPECIES_SURVEYS.keys())}")
        return
    
    survey = SPECIES_SURVEYS[survey_name]
    print(f"\n  {survey['description']}")
    
    pdb_dir = os.path.join(outdir, "pdb")
    af_dir = os.path.join(outdir, "alphafold")
    
    if survey.get("pdbs"):
        print(f"\n  --- PDB crystal structures ---")
        download_pdb_list(survey["pdbs"], pdb_dir)
    
    if survey.get("alphafold_uniprots"):
        print(f"\n  --- AlphaFold predicted structures ---")
        download_alphafold_list(survey["alphafold_uniprots"], af_dir)


def download_all_surveys(outdir):
    """Download everything for all species surveys."""
    for name in SPECIES_SURVEYS:
        print(f"\n{'='*60}")
        print(f"  SURVEY: {name}")
        print(f"{'='*60}")
        download_species_survey(name, os.path.join(outdir, name))


# ═══════════════════════════════════════════════════════════════
# GENERATE DOWNLOAD SCRIPTS (for offline use)
# ═══════════════════════════════════════════════════════════════

def generate_wget_script(outpath="download_structures.sh"):
    """Generate a shell script to download all structures with wget."""
    lines = ["#!/bin/bash", "# Auto-generated PDB/AlphaFold download script", 
             f"# Generated {time.strftime('%Y-%m-%d %H:%M')}",
             "", "set -e", ""]
    
    # All PDB IDs
    all_pdbs = set()
    all_uniprots = set()
    
    for pair in ARCHETYPE_PAIRS.values():
        if pair.get("apo"): all_pdbs.add(pair["apo"])
        if pair.get("holo"): all_pdbs.add(pair["holo"])
    
    for survey in SPECIES_SURVEYS.values():
        for p in survey.get("pdbs", []):
            all_pdbs.add(p)
        for u in survey.get("alphafold_uniprots", []):
            all_uniprots.add(u)
    
    lines.append(f"# {len(all_pdbs)} PDB structures + {len(all_uniprots)} AlphaFold models")
    lines.append("")
    lines.append("mkdir -p structures/pdb structures/alphafold")
    lines.append("")
    
    # PDB downloads
    lines.append("echo '=== Downloading PDB crystal structures ==='")
    for pdb in sorted(all_pdbs):
        lines.append(f'[ -f structures/pdb/{pdb}.cif ] || wget -q -O structures/pdb/{pdb}.cif "https://files.rcsb.org/download/{pdb}.cif" && echo "  {pdb}: OK" || echo "  {pdb}: FAILED"')
    
    lines.append("")
    lines.append("echo '=== Downloading AlphaFold models ==='")
    for uid in sorted(all_uniprots):
        lines.append(f'[ -f structures/alphafold/AF-{uid}-F1-model_v4.cif ] || wget -q -O structures/alphafold/AF-{uid}-F1-model_v4.cif "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.cif" && echo "  AF-{uid}: OK" || echo "  AF-{uid}: FAILED"')
    
    lines.append("")
    lines.append("echo 'Done!'")
    lines.append(f"echo \"Downloaded to structures/\"")
    
    with open(outpath, 'w') as f:
        f.write('\n'.join(lines))
    os.chmod(outpath, 0o755)
    print(f"  Generated {outpath} ({len(all_pdbs)} PDBs + {len(all_uniprots)} AlphaFold models)")
    return outpath


def generate_python_download_script(outpath="download_structures.py"):
    """Generate a standalone Python download script (no dependencies)."""
    
    all_pdbs = set()
    all_uniprots = set()
    
    for pair in ARCHETYPE_PAIRS.values():
        if pair.get("apo"): all_pdbs.add(pair["apo"])
        if pair.get("holo"): all_pdbs.add(pair["holo"])
    
    for survey in SPECIES_SURVEYS.values():
        for p in survey.get("pdbs", []):
            all_pdbs.add(p)
        for u in survey.get("alphafold_uniprots", []):
            all_uniprots.add(u)
    
    script = f'''#!/usr/bin/env python3
"""Download all structures for PhiPsi Atlas cross-species mapping."""
import os, urllib.request, time

PDB_IDS = {json.dumps(sorted(all_pdbs))}

UNIPROT_IDS = {json.dumps(sorted(all_uniprots))}

def dl(url, path):
    if os.path.exists(path): return True
    try:
        urllib.request.urlretrieve(url, path)
        return os.path.getsize(path) > 100
    except: return False

os.makedirs("structures/pdb", exist_ok=True)
os.makedirs("structures/alphafold", exist_ok=True)

print(f"Downloading {{len(PDB_IDS)}} PDB + {{len(UNIPROT_IDS)}} AlphaFold structures...")

ok = 0
for pdb in PDB_IDS:
    p = f"structures/pdb/{{pdb}}.cif"
    if dl(f"https://files.rcsb.org/download/{{pdb}}.cif", p):
        ok += 1; print(f"  {{pdb}}: OK")
    else:
        print(f"  {{pdb}}: FAILED")
    time.sleep(0.3)
print(f"PDB: {{ok}}/{{len(PDB_IDS)}}")

ok = 0
for uid in UNIPROT_IDS:
    p = f"structures/alphafold/AF-{{uid}}-F1-model_v4.cif"
    if dl(f"https://alphafold.ebi.ac.uk/files/AF-{{uid}}-F1-model_v4.cif", p):
        ok += 1; print(f"  AF-{{uid}}: OK")
    else:
        print(f"  AF-{{uid}}: FAILED")
    time.sleep(0.3)
print(f"AlphaFold: {{ok}}/{{len(UNIPROT_IDS)}}")
print("Done!")
'''
    
    with open(outpath, 'w') as f:
        f.write(script)
    os.chmod(outpath, 0o755)
    print(f"  Generated {outpath}")
    return outpath


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PDB & AlphaFold Structure Downloader")
    parser.add_argument("--mode", required=True,
        choices=["pdb", "alphafold", "family", "archetype-set", 
                 "species-survey", "all-surveys", "generate-scripts"],
        help="Download mode")
    parser.add_argument("--ids", nargs="+", help="PDB or UniProt IDs")
    parser.add_argument("--file", help="File with one ID per line")
    parser.add_argument("--query", help="Search query for family mode")
    parser.add_argument("--protein", help="Protein name for species survey")
    parser.add_argument("--max", type=int, default=50, help="Max results for search")
    parser.add_argument("--outdir", default="structures", help="Output directory")
    
    args = parser.parse_args()
    
    # Load IDs from file if specified
    ids = args.ids or []
    if args.file:
        with open(args.file) as f:
            ids.extend(line.strip() for line in f if line.strip() and not line.startswith('#'))
    
    if args.mode == "pdb":
        download_pdb_list(ids, args.outdir)
    
    elif args.mode == "alphafold":
        download_alphafold_list(ids, args.outdir)
    
    elif args.mode == "family":
        if not args.query:
            print("  --query required for family mode")
            sys.exit(1)
        pdb_ids = search_rcsb(args.query, args.max)
        if pdb_ids:
            download_pdb_list(pdb_ids, args.outdir)
    
    elif args.mode == "archetype-set":
        download_archetype_set(args.outdir)
    
    elif args.mode == "species-survey":
        if args.protein:
            download_species_survey(args.protein, args.outdir)
        else:
            print(f"  Available surveys: {', '.join(SPECIES_SURVEYS.keys())}")
    
    elif args.mode == "all-surveys":
        download_all_surveys(args.outdir)
    
    elif args.mode == "generate-scripts":
        generate_wget_script(os.path.join(args.outdir, "download_structures.sh"))
        generate_python_download_script(os.path.join(args.outdir, "download_structures.py"))
