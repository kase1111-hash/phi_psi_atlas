#!/usr/bin/env python3
"""
redundancy_clustering.py
========================
Full pipeline to test whether the Gini baseline (0.835 ± 0.005) survives
sequence redundancy removal at 80% and 50% identity.

Steps:
  1. Extract sequences from AlphaFold CIF files
  2. Install CD-HIT if needed
  3. Cluster at 80% and 50% identity
  4. Map cluster representatives back to census
  5. Recompute clade-level Gini statistics

Usage:
    python redundancy_clustering.py

Expects:
    G:\Mechanical Taxonomy\alphafold_cache\  (with species subdirectories)
    G:\Mechanical Taxonomy\all_species_census.csv

Outputs:
    G:\Mechanical Taxonomy\all_sequences.fasta
    G:\Mechanical Taxonomy\clusters_80.clstr
    G:\Mechanical Taxonomy\clusters_50.clstr
    G:\Mechanical Taxonomy\redundancy_report.txt
"""

import os
import sys
import csv
import subprocess
import re
import statistics
from pathlib import Path
from collections import defaultdict

BASE_DIR = r"G:\Mechanical Taxonomy"
CACHE_DIR = os.path.join(BASE_DIR, "alphafold_cache")
CENSUS_FILE = os.path.join(BASE_DIR, "all_species_census.csv")
FASTA_FILE = os.path.join(BASE_DIR, "all_sequences.fasta")
REPORT_FILE = os.path.join(BASE_DIR, "redundancy_report.txt")


# ─────────────────────────────────────────────────────────
# STEP 1: Extract sequences from CIF files
# ─────────────────────────────────────────────────────────

def extract_sequence_from_cif(cif_path):
    """Extract amino acid sequence from an AlphaFold mmCIF file.
    
    Looks for _entity_poly.pdbx_seq_one_letter_code or falls back
    to reading ATOM records for CA atoms.
    """
    seq = None
    try:
        with open(cif_path, 'r', encoding='utf-8', errors='replace') as f:
            in_seq = False
            seq_lines = []
            for line in f:
                # Method 1: _entity_poly.pdbx_seq_one_letter_code
                if '_entity_poly.pdbx_seq_one_letter_code' in line:
                    # Could be single-line or multi-line
                    parts = line.strip().split(None, 1)
                    if len(parts) > 1 and parts[1] not in (';', ''):
                        # Single line: value after tag
                        val = parts[1].strip().strip("'\"")
                        if val and val != '?':
                            seq = val.replace('\n', '').replace(' ', '')
                            break
                    else:
                        in_seq = True
                        continue
                
                if in_seq:
                    stripped = line.strip()
                    if stripped == ';':
                        if seq_lines:
                            # End of multi-line value
                            seq = ''.join(seq_lines).replace('\n', '').replace(' ', '')
                            break
                        else:
                            # Start of multi-line value
                            continue
                    else:
                        seq_lines.append(stripped)
    except Exception:
        pass
    
    # Clean sequence: remove non-letter characters
    if seq:
        seq = re.sub(r'[^A-Z]', '', seq.upper())
        if len(seq) >= 30:  # minimum length filter
            return seq
    
    # Method 2: Fall back to ATOM CA records
    try:
        residues = {}
        with open(cif_path, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('ATOM') and ' CA ' in line:
                    parts = line.split()
                    # mmCIF ATOM format varies; try to extract residue
                    try:
                        # Typical: ATOM id type_symbol label_atom_id ... label_comp_id label_seq_id ...
                        for i, p in enumerate(parts):
                            if p in ('ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY',
                                     'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
                                     'THR','TRP','TYR','VAL'):
                                seq_id = int(parts[i+2]) if i+2 < len(parts) else len(residues)+1
                                residues[seq_id] = p
                                break
                    except (ValueError, IndexError):
                        pass
        
        if len(residues) >= 30:
            aa_map = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C',
                      'GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
                      'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P',
                      'SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
            sorted_res = sorted(residues.items())
            seq = ''.join(aa_map.get(r, 'X') for _, r in sorted_res)
            return seq
    except Exception:
        pass
    
    return None


def extract_all_sequences():
    """Extract sequences from all CIF files in the cache."""
    print("STEP 1: Extracting sequences from CIF files...")
    
    # Find all species directories
    species_dirs = []
    for item in os.listdir(CACHE_DIR):
        path = os.path.join(CACHE_DIR, item)
        if os.path.isdir(path):
            species_dirs.append((item, path))
    
    if not species_dirs:
        # Maybe CIF files are directly in cache_dir
        species_dirs = [("all", CACHE_DIR)]
    
    print(f"  Found {len(species_dirs)} species directories")
    
    sequences = {}  # uid -> (species, sequence)
    total_cif = 0
    total_extracted = 0
    species_counts = {}
    
    for species_name, species_path in sorted(species_dirs):
        cif_files = [f for f in os.listdir(species_path) if f.endswith('.cif')]
        total_cif += len(cif_files)
        extracted = 0
        
        for fname in cif_files:
            uid = fname.replace('AF-', '').split('-F1')[0].split('-F')[0]
            seq = extract_sequence_from_cif(os.path.join(species_path, fname))
            if seq:
                sequences[uid] = (species_name, seq)
                extracted += 1
        
        total_extracted += extracted
        species_counts[species_name] = (len(cif_files), extracted)
        print(f"  {species_name}: {extracted}/{len(cif_files)} sequences extracted")
    
    print(f"\n  Total: {total_extracted}/{total_cif} sequences extracted")
    
    # Write FASTA
    print(f"  Writing FASTA to {FASTA_FILE}...")
    with open(FASTA_FILE, 'w') as f:
        for uid, (species, seq) in sequences.items():
            f.write(f">{uid}|{species}\n")
            # Write in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    
    print(f"  Done: {len(sequences)} sequences in FASTA file")
    return sequences, species_counts


# ─────────────────────────────────────────────────────────
# STEP 2: Install CD-HIT
# ─────────────────────────────────────────────────────────

def install_cdhit():
    """Try to install CD-HIT via conda or direct download."""
    print("\nSTEP 2: Installing CD-HIT...")
    
    # Check if already available
    try:
        result = subprocess.run(['cd-hit', '-h'], capture_output=True, text=True)
        if result.returncode == 0 or 'CD-HIT' in result.stdout or 'CD-HIT' in result.stderr:
            print("  CD-HIT already installed!")
            return True
    except FileNotFoundError:
        pass
    
    # Try conda
    try:
        print("  Trying: conda install -y -c bioconda cd-hit")
        result = subprocess.run(
            ['conda', 'install', '-y', '-c', 'bioconda', 'cd-hit'],
            capture_output=True, text=True, timeout=300
        )
        if result.returncode == 0:
            print("  CD-HIT installed via conda!")
            return True
        else:
            print(f"  conda failed: {result.stderr[:200]}")
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print("  conda not available or timed out")
    
    # Try pip (some systems have a wrapper)
    try:
        print("  Trying: pip install cd-hit")
        result = subprocess.run(
            [sys.executable, '-m', 'pip', 'install', 'cd-hit'],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            print("  CD-HIT installed via pip!")
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    # Try winget / choco on Windows
    if sys.platform == 'win32':
        print("  CD-HIT is not natively available on Windows.")
        print("  Options:")
        print("    1. Install WSL and run: sudo apt install cd-hit")
        print("    2. Download from https://github.com/weizhongli/cdhit/releases")
        print("    3. Use the built-in Python fallback (slower but works)")
        print("  Falling back to Python implementation...")
        return False
    
    print("  Could not install CD-HIT. Using Python fallback.")
    return False


# ─────────────────────────────────────────────────────────
# STEP 2b: Python fallback clustering (no CD-HIT needed)
# ─────────────────────────────────────────────────────────

def simple_sequence_identity(seq1, seq2):
    """Quick identity estimate using k-mer overlap (not alignment).
    
    This is an approximation: shared 3-mers / total 3-mers.
    Overestimates identity for short sequences, but sufficient
    for 80% and 50% clustering thresholds.
    """
    k = 3
    if len(seq1) < k or len(seq2) < k:
        return 0.0
    
    kmers1 = set(seq1[i:i+k] for i in range(len(seq1)-k+1))
    kmers2 = set(seq2[i:i+k] for i in range(len(seq2)-k+1))
    
    if not kmers1 or not kmers2:
        return 0.0
    
    shared = len(kmers1 & kmers2)
    total = min(len(kmers1), len(kmers2))
    
    return shared / total if total > 0 else 0.0


def greedy_cluster(sequences, threshold=0.80, max_proteins=None):
    """Greedy incremental clustering (same algorithm as CD-HIT).
    
    For each sequence, compare to existing cluster representatives.
    If identity >= threshold, add to that cluster. Otherwise, start new cluster.
    Sequences sorted by length (longest first) as in CD-HIT.
    
    For large datasets, we subsample for speed.
    """
    items = list(sequences.items())
    items.sort(key=lambda x: -len(x[1][1]))  # sort by sequence length, longest first
    
    if max_proteins and len(items) > max_proteins:
        # For very large datasets, sample uniformly
        import random
        random.seed(42)
        step = len(items) // max_proteins
        items = items[::step][:max_proteins]
        print(f"    Subsampled to {len(items)} proteins for clustering")
    
    clusters = {}  # representative_uid -> [member_uids]
    rep_seqs = {}  # representative_uid -> sequence
    
    for i, (uid, (species, seq)) in enumerate(items):
        if (i + 1) % 1000 == 0:
            print(f"    Processed {i+1}/{len(items)}, {len(clusters)} clusters")
        
        assigned = False
        
        # Compare to existing representatives (check length ratio first for speed)
        for rep_uid, rep_seq in rep_seqs.items():
            # CD-HIT length filter: shorter/longer must be >= threshold
            len_ratio = min(len(seq), len(rep_seq)) / max(len(seq), len(rep_seq))
            if len_ratio < threshold:
                continue
            
            identity = simple_sequence_identity(seq, rep_seq)
            if identity >= threshold:
                clusters[rep_uid].append(uid)
                assigned = True
                break
        
        if not assigned:
            clusters[uid] = [uid]
            rep_seqs[uid] = seq
    
    return clusters


# ─────────────────────────────────────────────────────────
# STEP 3: Run clustering
# ─────────────────────────────────────────────────────────

def run_cdhit(fasta_file, threshold, output_prefix):
    """Run CD-HIT clustering."""
    out_file = output_prefix + ".fasta"
    clstr_file = output_prefix + ".fasta.clstr"
    
    word_size = 5 if threshold >= 0.7 else (4 if threshold >= 0.6 else 3)
    
    cmd = [
        'cd-hit',
        '-i', fasta_file,
        '-o', out_file,
        '-c', str(threshold),
        '-n', str(word_size),
        '-M', '4000',  # 4GB memory
        '-T', '8',     # 8 threads
        '-d', '0',
    ]
    
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    
    if result.returncode != 0:
        print(f"  CD-HIT failed: {result.stderr[:300]}")
        return None
    
    # Parse cluster file
    clusters = {}
    current_rep = None
    with open(clstr_file, 'r') as f:
        for line in f:
            if line.startswith('>Cluster'):
                current_rep = None
            elif '*' in line:
                # Representative sequence
                uid = line.split('>')[1].split('|')[0].split('.')[0]
                current_rep = uid
                clusters[uid] = [uid]
            elif current_rep:
                uid = line.split('>')[1].split('|')[0].split('.')[0]
                clusters[current_rep].append(uid)
    
    return clusters


# ─────────────────────────────────────────────────────────
# STEP 4: Load census and recompute stats
# ─────────────────────────────────────────────────────────

def load_census():
    """Load all_species_census.csv into a uid->data dict."""
    census = {}
    with open(CENSUS_FILE, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for r in reader:
            uid = r['label'].replace('AF-', '').split('-F1')[0].split('-F')[0]
            try:
                census[uid] = {
                    'gini': float(r['gini']),
                    'n_res': int(r['n_valid']),
                    'frac_a': float(r['frac_alpha']),
                    'frac_b': float(r['frac_beta']),
                    'kappa2_mean': float(r['kappa2_mean']),
                }
            except:
                pass
    return census


def compute_stats(uid_set, census, species_map):
    """Compute Gini statistics for a set of UIDs, stratified by species."""
    species_ginis = defaultdict(list)
    all_ginis = []
    
    for uid in uid_set:
        if uid in census:
            g = census[uid]['gini']
            all_ginis.append(g)
            if uid in species_map:
                species_ginis[species_map[uid]].append(g)
    
    return all_ginis, species_ginis


# ─────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────

def main():
    report_lines = []
    def report(msg):
        print(msg)
        report_lines.append(msg)
    
    report("=" * 70)
    report("  REDUNDANCY CLUSTERING ANALYSIS")
    report("  Testing Gini baseline robustness to sequence redundancy")
    report("=" * 70)
    
    # Step 1: Extract sequences
    if os.path.exists(FASTA_FILE):
        report(f"\n  FASTA file exists: {FASTA_FILE}")
        report("  Counting sequences...")
        n_seq = 0
        sequences = {}
        species_counts = {}
        with open(FASTA_FILE, 'r') as f:
            uid = None
            seq_lines = []
            for line in f:
                if line.startswith('>'):
                    if uid and seq_lines:
                        seq = ''.join(seq_lines)
                        sp = uid.split('|')[1] if '|' in uid else 'unknown'
                        sequences[uid.split('|')[0]] = (sp, seq)
                    parts = line[1:].strip().split('|')
                    uid = line[1:].strip()
                    seq_lines = []
                    n_seq += 1
                else:
                    seq_lines.append(line.strip())
            if uid and seq_lines:
                seq = ''.join(seq_lines)
                sp = uid.split('|')[1] if '|' in uid else 'unknown'
                sequences[uid.split('|')[0]] = (sp, seq)
        report(f"  Found {n_seq} sequences")
    else:
        sequences, species_counts = extract_all_sequences()
    
    if not sequences:
        report("ERROR: No sequences extracted!")
        return
    
    # Build species map
    species_map = {uid: sp for uid, (sp, _) in sequences.items()}
    
    # Load census
    report("\nLoading census data...")
    census = load_census()
    report(f"  Census entries: {len(census)}")
    
    # Baseline stats (all proteins)
    all_ginis, sp_ginis = compute_stats(set(census.keys()), census, species_map)
    report(f"\n  BASELINE (all {len(all_ginis)} proteins):")
    report(f"    Mean Gini: {statistics.mean(all_ginis):.4f}")
    report(f"    Median:    {statistics.median(all_ginis):.4f}")
    report(f"    SD:        {statistics.stdev(all_ginis):.4f}")
    
    # Step 2: Try CD-HIT, fall back to Python
    has_cdhit = install_cdhit()
    
    for threshold, label in [(0.80, "80%"), (0.50, "50%")]:
        report(f"\n{'='*70}")
        report(f"  CLUSTERING AT {label} IDENTITY")
        report(f"{'='*70}")
        
        prefix = os.path.join(BASE_DIR, f"clusters_{int(threshold*100)}")
        
        if has_cdhit:
            clusters = run_cdhit(FASTA_FILE, threshold, prefix)
        else:
            report(f"\n  Using Python greedy clustering (approximate)...")
            # For the full dataset, limit to manageable size
            max_n = 30000 if threshold >= 0.7 else 20000
            clusters = greedy_cluster(sequences, threshold=threshold, max_proteins=max_n)
        
        if not clusters:
            report("  Clustering failed!")
            continue
        
        # Get representative UIDs
        rep_uids = set(clusters.keys())
        total_members = sum(len(v) for v in clusters.values())
        
        report(f"\n  Results:")
        report(f"    Input proteins:        {total_members}")
        report(f"    Clusters (non-redundant): {len(clusters)}")
        report(f"    Reduction:             {100*(1 - len(clusters)/total_members):.1f}%")
        report(f"    Largest cluster:       {max(len(v) for v in clusters.values())} members")
        
        # Compute Gini stats on representatives only
        rep_ginis, rep_sp_ginis = compute_stats(rep_uids, census, species_map)
        
        if rep_ginis:
            report(f"\n  GINI BASELINE (non-redundant at {label}):")
            report(f"    n proteins:  {len(rep_ginis)}")
            report(f"    Mean Gini:   {statistics.mean(rep_ginis):.4f}")
            report(f"    Median:      {statistics.median(rep_ginis):.4f}")
            report(f"    SD:          {statistics.stdev(rep_ginis):.4f}")
            
            # Singularity rate
            sing_rate = 100 * len([g for g in rep_ginis if g >= 0.92]) / len(rep_ginis)
            report(f"    Sing. rate:  {sing_rate:.1f}%")
            
            # Per-species stats (if available)
            if rep_sp_ginis:
                species_means = {}
                report(f"\n  Per-species Gini (non-redundant at {label}):")
                report(f"  {'Species':>25s} {'n':>6s} {'Mean Gini':>10s} {'SD':>8s}")
                for sp in sorted(rep_sp_ginis.keys()):
                    vals = rep_sp_ginis[sp]
                    if len(vals) >= 50:
                        m = statistics.mean(vals)
                        s = statistics.stdev(vals) if len(vals) > 1 else 0
                        species_means[sp] = m
                        report(f"  {sp:>25s} {len(vals):>6d} {m:>10.4f} {s:>8.4f}")
                
                if len(species_means) >= 3:
                    cv = statistics.stdev(species_means.values()) / statistics.mean(species_means.values())
                    report(f"\n  Cross-species CV of mean Gini: {100*cv:.2f}%")
                    report(f"  (Original baseline CV: 0.54%)")
    
    # Write report
    with open(REPORT_FILE, 'w') as f:
        f.write('\n'.join(report_lines))
    report(f"\nReport written to: {REPORT_FILE}")


if __name__ == '__main__':
    main()
