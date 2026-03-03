#!/usr/bin/env python3
"""Download all structures for PhiPsi Atlas cross-species mapping."""
import os, urllib.request, time

PDB_IDS = ["193L", "1AKE", "1AKI", "1AP4", "1ATP", "1BG2", "1BIN", "1BMF", "1BZP", "1BZR", "1CDK", "1CFD", "1CK7", "1CLL", "1D8U", "1DLY", "1EA5", "1ERK", "1EVE", "1F54", "1F55", "1FDH", "1FMV", "1FSX", "1FTJ", "1FTK", "1G09", "1GD6", "1GDI", "1HCK", "1HEL", "1HHO", "1HHP", "1HSG", "1IDR", "1IEP", "1IG8", "1IOM", "1IWT", "1IXE", "1JWR", "1LYZ", "1LZA", "1MBD", "1MOH", "1MUO", "1NME", "1P0I", "1P0M", "1P38", "1PAU", "1PPB", "1RJV", "1S0Q", "1SPY", "1UDT", "1UTN", "1V4S", "1V4T", "1VOM", "1W7J", "1W8J", "1XCK", "1Y57", "1Z83", "1ZIN", "1ZIO", "2BCX", "2C95", "2CTS", "2DN1", "2DN2", "2ERK", "2HHB", "2HYY", "2J4Z", "2LHB", "2LZM", "2SRC", "2YHX", "3CLN", "3GI3", "3KIN", "4AKE", "4CTS", "4EY4", "4EY7", "4HHB", "4PEP", "6LYZ"]

UNIPROT_IDS = ["P00519", "P00698", "P00699", "P00703", "P01942", "P02062", "P02144", "P02185", "P02189", "P02585", "P06787", "P07339", "P07858", "P0DP23", "P12931", "P24941", "P28482", "P33176", "P35579", "P60709", "P61626", "P62158", "P63098", "P68871", "P69905", "P80946", "Q16539"]

def dl(url, path):
    if os.path.exists(path): return True
    try:
        urllib.request.urlretrieve(url, path)
        return os.path.getsize(path) > 100
    except: return False

os.makedirs("structures/pdb", exist_ok=True)
os.makedirs("structures/alphafold", exist_ok=True)

print(f"Downloading {len(PDB_IDS)} PDB + {len(UNIPROT_IDS)} AlphaFold structures...")

ok = 0
for pdb in PDB_IDS:
    p = f"structures/pdb/{pdb}.cif"
    if dl(f"https://files.rcsb.org/download/{pdb}.cif", p):
        ok += 1; print(f"  {pdb}: OK")
    else:
        print(f"  {pdb}: FAILED")
    time.sleep(0.3)
print(f"PDB: {ok}/{len(PDB_IDS)}")

ok = 0
for uid in UNIPROT_IDS:
    p = f"structures/alphafold/AF-{uid}-F1-model_v4.cif"
    if dl(f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.cif", p):
        ok += 1; print(f"  AF-{uid}: OK")
    else:
        print(f"  AF-{uid}: FAILED")
    time.sleep(0.3)
print(f"AlphaFold: {ok}/{len(UNIPROT_IDS)}")
print("Done!")
