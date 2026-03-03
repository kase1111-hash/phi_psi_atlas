# baseline_disease_rate.py
# Queries UniProt for total human reviewed proteins with disease annotations
# Uses UniProt's REST API — single query, no per-protein calls

import urllib.request
import json

# Count ALL reviewed human proteins
url1 = "https://rest.uniprot.org/uniprotkb/search?query=(organism_id:9606)+AND+(reviewed:true)&format=json&size=1"
with urllib.request.urlopen(url1) as r:
    data = json.loads(r.read())
    total = int(r.headers.get('x-total-results', 0))
    print(f"Total reviewed human proteins: {total}")

# Count reviewed human proteins WITH disease annotation
url2 = "https://rest.uniprot.org/uniprotkb/search?query=(organism_id:9606)+AND+(reviewed:true)+AND+(cc_disease:*)&format=json&size=1"
with urllib.request.urlopen(url2) as r:
    data = json.loads(r.read())
    disease = int(r.headers.get('x-total-results', 0))
    print(f"With disease annotation: {disease}")
    print(f"Disease rate: {100*disease/total:.1f}%")
    print(f"Your singularity disease rate: 23.2%")
    print(f"Enrichment: {23.2/(100*disease/total):.2f}x")