# quick_check.py — run in G:\Mechanical Taxonomy\
import csv
targets = ['MYH6','MYH7','ATP4A','ATP4B','HEPACAM','HSPD1','HLA-B','GBA1','COX11']
with open('mimicry_candidates.csv') as f:
    for r in csv.DictReader(f):
        if r['human_gene'] in targets:
            print(f"{r['human_gene']} | {r['mimic_species']} | {r['source_type']} | score={r['score']} | {r['mimic_uid']}")