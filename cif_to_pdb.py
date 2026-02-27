#!/usr/bin/env python3
"""
Convert AlphaFold mmCIF files to PDB format for DALI submission.

DALI often rejects mmCIF files from AlphaFold. This script converts them
to standard PDB format that DALI accepts.

Usage:
    python cif_to_pdb.py                    # Convert all + rebuild pair folders
    python cif_to_pdb.py dali_structures/AF-P00918-F1.cif   # Convert one file

Requires: No external dependencies (pure Python mmCIF parser).
"""

import re
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent
CIF_DIR = SCRIPT_DIR / "dali_structures"
PDB_DIR = SCRIPT_DIR / "dali_pdb"
PAIR_DIR = SCRIPT_DIR / "dali_pairs"

PAIRS = [
    (1,  "P00918", "P59991", "CA2",     "KRTAP12-2",  0.804),
    (2,  "P60174", "P35270", "TPI",     "SPR",        0.798),
    (3,  "P62258", "O43808", "14-3-3e", "SLC25A17",   0.783),
    (4,  "P11166", "P09914", "GLUT1",   "IFIT1",      0.783),
    (5,  "P02751", "P21333", "FN1",     "FLNA",       0.779),
    (6,  "P01857", "P32969", "IgG1",    "RPL9",       0.778),
    (7,  "P38646", "O15344", "GRP75",   "MID1",       0.767),
    (8,  "P07900", "O14829", "HSP90a",  "PPEF1",      0.765),
    (9,  "P04626", "O95302", "HER2",    "FKBP9",      0.757),
    (10, "P21860", "P28799", "ErbB3",   "GRN",        0.757),
]


def parse_mmcif_atoms(cif_path):
    """Parse _atom_site records from mmCIF file.
    
    Returns list of dicts with atom properties.
    """
    with open(cif_path, 'r') as f:
        lines = f.readlines()

    atoms = []
    in_atom_site = False
    columns = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Detect start of _atom_site loop
        if line == 'loop_':
            # Check if next lines are _atom_site.xxx
            j = i + 1
            cols = []
            while j < len(lines) and lines[j].strip().startswith('_atom_site.'):
                cols.append(lines[j].strip())
                j += 1
            if cols:
                columns = [c.replace('_atom_site.', '') for c in cols]
                in_atom_site = True
                i = j
                continue

        if in_atom_site:
            if not line or line.startswith('_') or line.startswith('#') or line.startswith('loop_'):
                in_atom_site = False
                i += 1
                continue

            # Parse data line — handle quoted strings
            parts = _split_cif_line(line)
            if len(parts) >= len(columns):
                atom = {}
                for k, col in enumerate(columns):
                    atom[col] = parts[k] if k < len(parts) else '.'
                atoms.append(atom)

        i += 1

    return atoms


def _split_cif_line(line):
    """Split a CIF data line respecting quoted strings."""
    parts = []
    current = ''
    in_quote = False
    quote_char = None

    for ch in line:
        if in_quote:
            if ch == quote_char:
                in_quote = False
                parts.append(current)
                current = ''
            else:
                current += ch
        elif ch in ('"', "'"):
            in_quote = True
            quote_char = ch
            current = ''
        elif ch in (' ', '\t'):
            if current:
                parts.append(current)
                current = ''
        else:
            current += ch

    if current:
        parts.append(current)

    return parts


def atoms_to_pdb(atoms, title="DALI"):
    """Convert parsed atoms to PDB format string."""
    pdb_lines = []
    pdb_lines.append(f"HEADER    {title:40s}")
    
    serial = 0
    prev_chain = None
    prev_resi = None

    for atom in atoms:
        group = atom.get('group_PDB', 'ATOM')
        if group not in ('ATOM', 'HETATM'):
            continue

        serial += 1
        if serial > 99999:
            serial = 1  # PDB serial wraps

        atom_name = atom.get('label_atom_id', atom.get('auth_atom_id', 'CA'))
        alt_loc = atom.get('label_alt_id', '.')
        if alt_loc == '.':
            alt_loc = ' '
        
        res_name = atom.get('label_comp_id', atom.get('auth_comp_id', 'UNK'))
        chain = atom.get('label_asym_id', atom.get('auth_asym_id', 'A'))
        if chain == '.':
            chain = 'A'
        chain = chain[0]  # PDB only supports single char

        res_seq = atom.get('label_seq_id', atom.get('auth_seq_id', '1'))
        if res_seq == '.':
            res_seq = '0'
        
        icode = atom.get('pdbx_PDB_ins_code', ' ')
        if icode == '?' or icode == '.':
            icode = ' '

        x = float(atom.get('Cartn_x', 0))
        y = float(atom.get('Cartn_y', 0))
        z = float(atom.get('Cartn_z', 0))
        
        occupancy = float(atom.get('occupancy', 1.0))
        b_factor = float(atom.get('B_iso_or_equiv', 0.0))
        
        element = atom.get('type_symbol', atom_name[0])

        # Format atom name (4 chars, left-justified if len<=3, right-shifted if 4)
        if len(atom_name) < 4:
            atom_name_fmt = f" {atom_name:<3s}"
        else:
            atom_name_fmt = f"{atom_name:<4s}"

        # Standard PDB ATOM record format
        # COLUMNS        DATA TYPE       FIELD         DEFINITION
        # 1-6            Record name     "ATOM  "
        # 7-11           Integer         serial        Atom serial number
        # 13-16          Atom            name          Atom name
        # 17             Character       altLoc        Alternate location
        # 18-20          Residue name    resName       Residue name
        # 22             Character       chainID       Chain identifier
        # 23-26          Integer         resSeq        Residue sequence number
        # 27             AChar           iCode         Code for insertion of residues
        # 31-38          Real(8.3)       x             X coordinate
        # 39-46          Real(8.3)       y             Y coordinate
        # 47-54          Real(8.3)       z             Z coordinate
        # 55-60          Real(6.2)       occupancy     Occupancy
        # 61-66          Real(6.2)       tempFactor    Temperature factor
        # 77-78          LString(2)      element       Element symbol

        pdb_line = (f"{group:<6s}{serial:5d} {atom_name_fmt}{alt_loc}"
                    f"{res_name:>3s} {chain}{int(res_seq):4d}{icode}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}"
                    f"{occupancy:6.2f}{b_factor:6.2f}"
                    f"          {element:>2s}  ")

        # TER record between chains
        if prev_chain is not None and chain != prev_chain:
            pdb_lines.append("TER")
        
        pdb_lines.append(pdb_line)
        prev_chain = chain
        prev_resi = res_seq

    pdb_lines.append("TER")
    pdb_lines.append("END")

    return '\n'.join(pdb_lines) + '\n'


def convert_one(cif_path, output_path=None):
    """Convert one CIF file to PDB format."""
    cif_path = Path(cif_path)
    if output_path is None:
        output_path = PDB_DIR / cif_path.with_suffix('.pdb').name
    else:
        output_path = Path(output_path)

    atoms = parse_mmcif_atoms(cif_path)
    if not atoms:
        print(f"  ERROR: No atoms found in {cif_path.name}")
        return False

    title = cif_path.stem
    pdb_text = atoms_to_pdb(atoms, title=title)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(pdb_text)

    n_atoms = sum(1 for a in atoms if a.get('group_PDB') in ('ATOM', 'HETATM'))
    print(f"  {cif_path.name} → {output_path.name}  ({n_atoms} atoms)")
    return True


def find_cif(uid):
    """Find CIF file for a UID."""
    for pattern in [f"AF-{uid}-F1.cif", f"AF-{uid}-F1-model_v6.cif", f"AF-{uid}-F1-model_v4.cif"]:
        p = CIF_DIR / pattern
        if p.exists() and p.stat().st_size > 1000:
            return p
    return None


def convert_all_and_organize():
    """Convert all CIFs to PDB and organize into pair folders."""
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    PAIR_DIR.mkdir(parents=True, exist_ok=True)

    print(f"  Converting CIF → PDB for all DALI pairs...\n")

    import shutil

    for num, q_uid, h_uid, q_gene, h_gene, geom in PAIRS:
        q_cif = find_cif(q_uid)
        h_cif = find_cif(h_uid)

        if not q_cif:
            print(f"  Pair {num}: MISSING CIF for {q_gene} ({q_uid})")
            continue
        if not h_cif:
            print(f"  Pair {num}: MISSING CIF for {h_gene} ({h_uid})")
            continue

        print(f"  ── Pair {num}: {q_gene} vs {h_gene} (geom={geom:.3f}) ──")

        # Convert to PDB
        q_pdb = PDB_DIR / f"{q_gene}_{q_uid}.pdb"
        h_pdb = PDB_DIR / f"{h_gene}_{h_uid}.pdb"

        convert_one(q_cif, q_pdb)
        convert_one(h_cif, h_pdb)

        # Create pair folder with clear names
        pair_folder = PAIR_DIR / f"pair{num:02d}_{q_gene}_vs_{h_gene}"
        pair_folder.mkdir(parents=True, exist_ok=True)

        dst_q = pair_folder / f"structure1_{q_gene}_{q_uid}.pdb"
        dst_h = pair_folder / f"structure2_{h_gene}_{h_uid}.pdb"

        shutil.copy2(q_pdb, dst_q)
        shutil.copy2(h_pdb, dst_h)
        print()

    print(f"  Done! PDB files in: {PDB_DIR}/")
    print(f"  Pair folders in:    {PAIR_DIR}/")
    print()
    print(f"  Each pair folder contains:")
    print(f"    structure1_<GENE>_<UID>.pdb  ← upload as Structure 1 in DALI")
    print(f"    structure2_<GENE>_<UID>.pdb  ← upload as Structure 2 in DALI")
    print()
    print(f"  DALI server: http://ekhidna2.biocenter.helsinki.fi/dali/")
    print(f"  After getting Z-scores: python dali_batch.py manual")


def main():
    args = sys.argv[1:]

    if args and Path(args[0]).exists():
        # Convert single file
        convert_one(args[0])
    else:
        convert_all_and_organize()


if __name__ == "__main__":
    main()
