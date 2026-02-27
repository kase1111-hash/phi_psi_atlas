#!/usr/bin/env python3
"""
Batch DALI structural superposition validation.

Downloads AlphaFold CIF files for 10 geometry-only pairs and submits them
to the DALI server for pairwise structural comparison.

Usage:
    python dali_batch.py download          # Download all CIF files
    python dali_batch.py submit            # Submit all pairs to DALI
    python dali_batch.py submit --pair 2   # Submit just pair #2
    python dali_batch.py status            # Check submission status
    python dali_batch.py parse             # Parse all results
    python dali_batch.py manual            # Enter Z-scores manually
    python dali_batch.py report            # Show final report

Requires: requests (pip install requests)
"""

import json
import os
import re
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
DATA_DIR = SCRIPT_DIR / "data"
CIF_DIR = SCRIPT_DIR / "dali_structures"
RESULTS_DIR = SCRIPT_DIR / "dali_results"

# ═══════════════════════════════════════════════════════════════════════
#  PAIR DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════

PAIRS = [
    {"id": 1,  "query": "P00918", "hit": "P59991", "q_name": "CA2",     "h_name": "KRTAP12-2",  "geom": 0.804, "note": "Enzyme vs keratin protein"},
    {"id": 2,  "query": "P60174", "hit": "P35270", "q_name": "TPI",     "h_name": "SPR",        "geom": 0.798, "note": "Both alpha/beta reductases"},
    {"id": 3,  "query": "P62258", "hit": "O43808", "q_name": "14-3-3e", "h_name": "SLC25A17",   "geom": 0.783, "note": "All-alpha vs TM transporter"},
    {"id": 4,  "query": "P11166", "hit": "P09914", "q_name": "GLUT1",   "h_name": "IFIT1",      "geom": 0.783, "note": "TM transporter vs antiviral"},
    {"id": 5,  "query": "P02751", "hit": "P21333", "q_name": "FN1",     "h_name": "FLNA",       "geom": 0.779, "note": "Both giant beta-repeat proteins"},
    {"id": 6,  "query": "P01857", "hit": "P32969", "q_name": "IgG1",    "h_name": "RPL9",       "geom": 0.778, "note": "Ig fold vs ribosomal protein"},
    {"id": 7,  "query": "P38646", "hit": "O15344", "q_name": "GRP75",   "h_name": "MID1",       "geom": 0.767, "note": "HSP70 vs E3 ubiquitin ligase"},
    {"id": 8,  "query": "P07900", "hit": "O14829", "q_name": "HSP90a",  "h_name": "PPEF1",      "geom": 0.765, "note": "Chaperone vs phosphatase"},
    {"id": 9,  "query": "P04626", "hit": "O95302", "q_name": "HER2",    "h_name": "FKBP9",      "geom": 0.757, "note": "RTK vs FK-binding protein"},
    {"id": 10, "query": "P21860", "hit": "P28799", "q_name": "ErbB3",   "h_name": "GRN",        "geom": 0.757, "note": "RTK vs growth factor"},
]


def alphafold_url(uid):
    """AlphaFold CIF download URL — try v6 first, then v4."""
    return f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v6.cif"


def alphafold_api_url(uid):
    """AlphaFold API URL that returns JSON with the real CIF link."""
    return f"https://alphafold.com/api/prediction/{uid}"


def cif_path(uid):
    """Local path for a CIF file."""
    return CIF_DIR / f"AF-{uid}-F1.cif"


# ═══════════════════════════════════════════════════════════════════════
#  DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════

def download_cif(uid, timeout=30):
    """Download a CIF file for one UniProt ID.
    
    Strategy:
    1. Try direct v6 CIF URL
    2. Try direct v4 CIF URL  
    3. Query API for metadata JSON, extract cifUrl, download that
    """
    outpath = cif_path(uid)
    if outpath.exists() and outpath.stat().st_size > 5000:
        return True  # Already have it and it's not junk

    headers = {'User-Agent': 'Mozilla/5.0 PhiPsiAtlas/1.0'}

    # Strategy 1: Direct v6 URL
    for version in ['v6', 'v4', 'v3']:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_{version}.cif"
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
                if data[:5] == b'data_':  # Valid CIF starts with data_
                    with open(outpath, 'wb') as f:
                        f.write(data)
                    print(f"    {uid}: downloaded ({len(data):,} bytes, {version})")
                    return True
        except Exception:
            pass

    # Strategy 2: Query API, parse JSON, follow cifUrl
    api_url = alphafold_api_url(uid)
    try:
        req = urllib.request.Request(api_url, headers=headers)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            api_data = resp.read().decode('utf-8')
            entries = json.loads(api_data)
            if isinstance(entries, list) and entries:
                cif_url = entries[0].get('cifUrl', '')
                if cif_url:
                    req2 = urllib.request.Request(cif_url, headers=headers)
                    with urllib.request.urlopen(req2, timeout=timeout) as resp2:
                        data = resp2.read()
                        if data[:5] == b'data_':
                            with open(outpath, 'wb') as f:
                                f.write(data)
                            print(f"    {uid}: downloaded via API ({len(data):,} bytes)")
                            return True
                        else:
                            print(f"    {uid}: API cifUrl returned non-CIF data")
    except Exception as e:
        print(f"    {uid}: API fallback failed: {e}")

    print(f"    {uid}: FAILED all download methods")
    return False


def cmd_download():
    """Download AlphaFold CIF files for all pairs."""
    CIF_DIR.mkdir(parents=True, exist_ok=True)

    # Collect unique UIDs
    uids = set()
    for p in PAIRS:
        uids.add(p["query"])
        uids.add(p["hit"])

    print(f"  Downloading {len(uids)} unique structures...")

    success = 0
    for uid in sorted(uids):
        if download_cif(uid):
            success += 1
        time.sleep(0.5)

    missing = [uid for uid in sorted(uids) if not cif_path(uid).exists() or cif_path(uid).stat().st_size < 5000]
    if missing:
        print(f"\n  WARNING: {len(missing)} structures missing or too small: {missing}")
        print(f"  Download manually from https://alphafold.ebi.ac.uk/")
    else:
        print(f"\n  All {len(uids)} structures downloaded successfully ({success} files).")


# ═══════════════════════════════════════════════════════════════════════
#  SUBMIT TO DALI
# ═══════════════════════════════════════════════════════════════════════

def submit_dali_pair(pair):
    """Submit a pair to the DALI server.
    
    Uses the DALI API at ekhidna2.biocenter.helsinki.fi.
    Returns the job URL if successful.
    """
    try:
        import requests
        HAS_REQUESTS = True
    except ImportError:
        HAS_REQUESTS = False

    q_cif = cif_path(pair["query"])
    h_cif = cif_path(pair["hit"])

    # Also check for old v4 filenames
    if not q_cif.exists():
        q_cif_v4 = CIF_DIR / f"AF-{pair['query']}-F1-model_v4.cif"
        if q_cif_v4.exists():
            q_cif = q_cif_v4
    if not h_cif.exists():
        h_cif_v4 = CIF_DIR / f"AF-{pair['hit']}-F1-model_v4.cif"
        if h_cif_v4.exists():
            h_cif = h_cif_v4

    if not q_cif.exists() or not h_cif.exists():
        print(f"    Missing CIF files for pair {pair['id']}. Run: python dali_batch.py download")
        return None

    if not HAS_REQUESTS:
        print(f"    'requests' library not installed. Using manual submission mode.")
        return None

    url = "http://ekhidna2.biocenter.helsinki.fi/cgi-bin/sans/dump.cgi"

    try:
        with open(q_cif, 'rb') as f1, open(h_cif, 'rb') as f2:
            files = {
                'file1': (q_cif.name, f1, 'chemical/x-cif'),
                'file2': (h_cif.name, f2, 'chemical/x-cif'),
            }
            data = {
                'method': 'pairwise',
                'title': f"PhiPsiAtlas_pair{pair['id']}_{pair['q_name']}_vs_{pair['h_name']}",
                'address': '',  # no email
            }
            resp = requests.post(url, files=files, data=data, timeout=60)

            if resp.status_code == 200:
                # Try to extract job URL from response
                match = re.search(r'(http://ekhidna2\.biocenter\.helsinki\.fi/[^\s<>"]+)', resp.text)
                if match:
                    job_url = match.group(1)
                    print(f"    Submitted. Job URL: {job_url}")
                    return job_url
                else:
                    print(f"    Submitted but could not parse job URL.")
                    # Save response for debugging
                    debug_file = RESULTS_DIR / f"pair{pair['id']}_response.html"
                    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
                    with open(debug_file, 'w') as df:
                        df.write(resp.text)
                    print(f"    Response saved to {debug_file}")
                    return "submitted"
            else:
                print(f"    HTTP {resp.status_code}")
                return None

    except Exception as e:
        print(f"    Submission failed: {e}")
        return None


def cmd_submit(pair_num=None):
    """Submit pairs to DALI server."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    status_file = RESULTS_DIR / "submission_status.json"

    # Load existing status
    if status_file.exists():
        with open(status_file) as f:
            status = json.load(f)
    else:
        status = {}

    pairs_to_submit = PAIRS if pair_num is None else [p for p in PAIRS if p["id"] == pair_num]

    for pair in pairs_to_submit:
        pid = str(pair["id"])
        if pid in status and status[pid].get("job_url"):
            print(f"  Pair {pair['id']} ({pair['q_name']} vs {pair['h_name']}): already submitted")
            continue

        print(f"  Pair {pair['id']}: {pair['q_name']} ({pair['query']}) vs "
              f"{pair['h_name']} ({pair['hit']})  [geom={pair['geom']:.3f}]")

        job_url = submit_dali_pair(pair)

        status[pid] = {
            "query": pair["query"],
            "hit": pair["hit"],
            "q_name": pair["q_name"],
            "h_name": pair["h_name"],
            "geom_score": pair["geom"],
            "job_url": job_url,
            "submitted_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "z_score": None,
        }
        time.sleep(2)  # Be polite to the server

    with open(status_file, 'w') as f:
        json.dump(status, f, indent=2)
    print(f"\n  Status saved to {status_file}")

    # If submission failed, show manual instructions
    failed = [p for p in pairs_to_submit if not status.get(str(p["id"]), {}).get("job_url")]
    if failed:
        print_manual_instructions(failed)


def print_manual_instructions(pairs=None):
    """Print manual DALI submission instructions."""
    pairs = pairs or PAIRS

    print(f"\n  ╔══════════════════════════════════════════════════════════════╗")
    print(f"  ║  MANUAL DALI SUBMISSION INSTRUCTIONS                        ║")
    print(f"  ╚══════════════════════════════════════════════════════════════╝")
    print()
    print(f"  1. Go to: http://ekhidna2.biocenter.helsinki.fi/dali/")
    print(f"  2. Click 'Pairwise structural alignment'")
    print(f"  3. For each pair below:")
    print(f"     - Upload Structure 1 (query) and Structure 2 (hit)")
    print(f"     - CIF files are in: {CIF_DIR}/")
    print(f"     - Click Submit")
    print(f"     - Record the Z-score from the results page")
    print(f"  4. Run: python dali_batch.py manual")
    print()
    print(f"  {'#':>3s} {'Query':>8s} {'→':>1s} {'Hit':>10s}  {'File 1 (query)':>35s}  {'File 2 (hit)':>35s}  {'Geom':>5s}")
    print(f"  {'-'*3} {'-'*8} {'-'*1} {'-'*10}  {'-'*35}  {'-'*35}  {'-'*5}")

    for p in pairs:
        f1 = f"AF-{p['query']}-F1-model_v4.cif"
        f2 = f"AF-{p['hit']}-F1-model_v4.cif"
        print(f"  {p['id']:3d} {p['q_name']:>8s} → {p['h_name']:>10s}  {f1:>35s}  {f2:>35s}  {p['geom']:5.3f}")

    print()
    print(f"  After getting Z-scores, run: python dali_batch.py manual")


# ═══════════════════════════════════════════════════════════════════════
#  MANUAL Z-SCORE ENTRY
# ═══════════════════════════════════════════════════════════════════════

def cmd_manual():
    """Enter DALI Z-scores manually."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    results_file = RESULTS_DIR / "dali_zscores.json"

    # Load existing
    if results_file.exists():
        with open(results_file) as f:
            existing = json.load(f)
    else:
        existing = {}

    print(f"\n  Enter DALI Z-scores for each pair.")
    print(f"  (Press Enter to skip, 'q' to quit)")
    print(f"  Z > 2: significant | Z > 8: same superfamily | Z > 20: very close")
    print()

    for p in PAIRS:
        pid = str(p["id"])
        existing_z = existing.get(pid, {}).get("z_score")
        skip_note = f"  [current: Z={existing_z}]" if existing_z is not None else ""

        prompt = (f"  Pair {p['id']:2d}: {p['q_name']:>8s} → {p['h_name']:>10s} "
                  f"(geom={p['geom']:.3f}) {p['note']}{skip_note}\n"
                  f"         Z-score = ")
        try:
            val = input(prompt).strip()
        except (EOFError, KeyboardInterrupt):
            print("\n  Stopped.")
            break

        if val.lower() == 'q':
            break
        if not val:
            continue

        try:
            z = float(val)
            existing[pid] = {
                "pair_id": p["id"],
                "query": p["query"],
                "hit": p["hit"],
                "q_name": p["q_name"],
                "h_name": p["h_name"],
                "geom_score": p["geom"],
                "z_score": z,
                "note": p["note"],
                "significant": z > 2,
                "superfamily": z > 8,
            }
            sig = "***" if z > 8 else "**" if z > 4 else "*" if z > 2 else "ns"
            print(f"         → Z={z:.1f} {sig}")
        except ValueError:
            print(f"         → Invalid, skipped")

    with open(results_file, 'w') as f:
        json.dump(existing, f, indent=2)
    print(f"\n  Saved to {results_file}")

    # Show summary
    cmd_report()


# ═══════════════════════════════════════════════════════════════════════
#  REPORT
# ═══════════════════════════════════════════════════════════════════════

def cmd_report():
    """Show DALI validation results and generate figure."""
    results_file = RESULTS_DIR / "dali_zscores.json"

    if not results_file.exists():
        print("  No results yet. Run: python dali_batch.py manual")
        return

    with open(results_file) as f:
        results = json.load(f)

    if not results:
        print("  No Z-scores recorded yet.")
        return

    print(f"\n  ╔══════════════════════════════════════════════════════════════╗")
    print(f"  ║  DALI VALIDATION RESULTS                                    ║")
    print(f"  ╚══════════════════════════════════════════════════════════════╝")
    print()
    print(f"  {'#':>3s} {'Query':>8s} {'→':>1s} {'Hit':>10s}  {'Geom':>5s}  {'DALI Z':>7s}  {'Sig':>4s}  {'Note'}")
    print(f"  {'-'*70}")

    n_sig = 0
    n_strong = 0
    n_total = 0
    geom_scores = []
    z_scores = []

    for p in PAIRS:
        pid = str(p["id"])
        r = results.get(pid, {})
        z = r.get("z_score")

        if z is not None:
            n_total += 1
            geom_scores.append(p["geom"])
            z_scores.append(z)

            if z > 2:
                n_sig += 1
            if z > 8:
                n_strong += 1

            sig = "***" if z > 8 else "**" if z > 4 else "*" if z > 2 else "ns"
            print(f"  {p['id']:3d} {p['q_name']:>8s} → {p['h_name']:>10s}  {p['geom']:5.3f}  {z:7.1f}  {sig:>4s}  {p['note']}")
        else:
            print(f"  {p['id']:3d} {p['q_name']:>8s} → {p['h_name']:>10s}  {p['geom']:5.3f}  {'---':>7s}  {'':>4s}  {p['note']}")

    if n_total > 0:
        print(f"\n  SUMMARY:")
        print(f"    Pairs tested:           {n_total}")
        print(f"    Significant (Z > 2):    {n_sig}/{n_total} ({n_sig/n_total*100:.0f}%)")
        print(f"    Strong (Z > 8):         {n_strong}/{n_total} ({n_strong/n_total*100:.0f}%)")

        if n_total >= 3:
            import numpy as np
            corr = np.corrcoef(geom_scores, z_scores)[0, 1]
            print(f"    Correlation (geom↔Z):   r = {corr:.3f}")
            mean_z = np.mean(z_scores)
            print(f"    Mean Z-score:           {mean_z:.1f}")

        # Generate figure
        try:
            _plot_dali(geom_scores, z_scores, results)
        except Exception as e:
            print(f"    (Figure generation failed: {e})")

    # Also save to the data dir for the paper pipeline
    paper_results = []
    for p in PAIRS:
        pid = str(p["id"])
        r = results.get(pid, {})
        if r.get("z_score") is not None:
            paper_results.append({
                "query": p["query"],
                "hit": p["hit"],
                "query_gene": p["q_name"],
                "hit_gene": p["h_name"],
                "geom_score": p["geom"],
                "dali_z": r["z_score"],
            })
    
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(DATA_DIR / "dali_results.json", 'w') as f:
        json.dump(paper_results, f, indent=2)


def _plot_dali(geom_scores, z_scores, results):
    """Generate Fig 7: geometric score vs DALI Z-score."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np

    fig, ax = plt.subplots(1, 1, figsize=(6, 5))

    colors = ['#27AE60' if z > 2 else '#E74C3C' for z in z_scores]
    ax.scatter(geom_scores, z_scores, c=colors, s=80, edgecolors='black',
               lw=0.5, zorder=5)

    # Label each point
    for p in PAIRS:
        pid = str(p["id"])
        r = results.get(pid, {})
        z = r.get("z_score")
        if z is not None:
            label = f"{p['q_name']}-{p['h_name']}"
            ax.annotate(label, (p["geom"], z), textcoords='offset points',
                       xytext=(5, 5), fontsize=7)

    ax.axhline(2, color='#E74C3C', ls='--', alpha=0.4, label='Z = 2 (significant)')
    ax.axhline(8, color='#27AE60', ls='--', alpha=0.4, label='Z = 8 (superfamily)')

    ax.set_xlabel('Geometric similarity score', fontsize=11)
    ax.set_ylabel('DALI Z-score', fontsize=11)
    ax.set_title('Geometry-only hits: structural validation', fontsize=12)
    ax.legend(fontsize=9, loc='upper left')

    # Correlation annotation
    if len(geom_scores) >= 3:
        corr = np.corrcoef(geom_scores, z_scores)[0, 1]
        ax.text(0.95, 0.05, f'r = {corr:.2f}', transform=ax.transAxes,
                fontsize=10, ha='right', va='bottom',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    fig.tight_layout()

    fig_dir = SCRIPT_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / 'fig7_dali_validation.svg', bbox_inches='tight')
    fig.savefig(fig_dir / 'fig7_dali_validation.png', dpi=300, bbox_inches='tight')
    print(f"    Saved: {fig_dir}/fig7_dali_validation.svg/png")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
#  STATUS CHECK
# ═══════════════════════════════════════════════════════════════════════

def cmd_status():
    """Check status of submitted DALI jobs."""
    status_file = RESULTS_DIR / "submission_status.json"
    if not status_file.exists():
        print("  No submissions yet. Run: python dali_batch.py submit")
        return

    with open(status_file) as f:
        status = json.load(f)

    print(f"\n  DALI SUBMISSION STATUS:")
    for pid, s in sorted(status.items(), key=lambda x: int(x[0])):
        url = s.get("job_url", "not submitted")
        z = s.get("z_score", "pending")
        print(f"  Pair {pid}: {s['q_name']} → {s['h_name']}  "
              f"z={z}  url={url}")


# ═══════════════════════════════════════════════════════════════════════
#  PARSE DALI OUTPUT FILES
# ═══════════════════════════════════════════════════════════════════════

def cmd_parse():
    """Parse DALI output files from dali_results/ directory.
    
    Looks for .txt files matching pair IDs and extracts Z-scores.
    DALI output format typically has a summary line like:
    
    # Structural alignment of mol1A vs mol2A
    # Z-score = 4.2  ...
    """
    results_file = RESULTS_DIR / "dali_zscores.json"

    if results_file.exists():
        with open(results_file) as f:
            existing = json.load(f)
    else:
        existing = {}

    found = 0
    for txt_file in sorted(RESULTS_DIR.glob("*.txt")):
        content = txt_file.read_text()

        # Try to extract Z-score
        # Common DALI patterns:
        #   Z-score = 4.2
        #   Z  rmsd  lali  ...
        #   4.2  2.1  120  ...
        z = None

        # Pattern 1: explicit "Z-score" or "Z ="
        m = re.search(r'Z[-_\s]*score\s*[=:]\s*([\d.]+)', content, re.IGNORECASE)
        if m:
            z = float(m.group(1))

        # Pattern 2: DALI summary table header followed by data
        if z is None:
            m = re.search(r'^\s*Z\s+rmsd\s+lali.*?\n\s*([\d.]+)', content, re.MULTILINE)
            if m:
                z = float(m.group(1))

        # Pattern 3: "No significant structural similarity" → Z < 2
        if z is None and 'no significant' in content.lower():
            z = 0.0

        if z is not None:
            found += 1
            # Try to match to a pair
            for p in PAIRS:
                if p["query"] in content or p["hit"] in content or \
                   p["q_name"] in txt_file.name or f"pair{p['id']}" in txt_file.name:
                    pid = str(p["id"])
                    existing[pid] = {
                        "pair_id": p["id"],
                        "query": p["query"],
                        "hit": p["hit"],
                        "q_name": p["q_name"],
                        "h_name": p["h_name"],
                        "geom_score": p["geom"],
                        "z_score": z,
                        "note": p["note"],
                        "significant": z > 2,
                        "superfamily": z > 8,
                        "source_file": txt_file.name,
                    }
                    print(f"  Pair {p['id']}: {p['q_name']}→{p['h_name']} Z={z:.1f} (from {txt_file.name})")
                    break

    if found:
        with open(results_file, 'w') as f:
            json.dump(existing, f, indent=2)
        print(f"\n  Parsed {found} Z-scores. Saved to {results_file}")
        cmd_report()
    else:
        print(f"  No DALI output files found in {RESULTS_DIR}/")
        print(f"  Save DALI results as .txt files in that directory, then re-run parse.")


# ═══════════════════════════════════════════════════════════════════════
#  MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    args = sys.argv[1:] or ['download']
    cmd = args[0]

    if cmd == 'download':
        cmd_download()

    elif cmd == 'submit':
        pair_num = None
        for i, a in enumerate(args):
            if a == '--pair' and i + 1 < len(args):
                pair_num = int(args[i + 1])
        cmd_submit(pair_num)

    elif cmd == 'status':
        cmd_status()

    elif cmd == 'parse':
        cmd_parse()

    elif cmd == 'manual':
        cmd_manual()

    elif cmd == 'report':
        cmd_report()

    elif cmd == 'instructions':
        print_manual_instructions()

    else:
        print(f"  Unknown command: {cmd}")
        print(f"  Commands: download, submit, status, parse, manual, report, instructions")


if __name__ == "__main__":
    main()
