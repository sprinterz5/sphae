"""
Check merged GBK against embedding cache, output only proteins not yet embedded.
Cache key: MD5(aa_sequence) — same sequence across different phages = one cache hit.
Args: merged_gbk cache_3di_fasta new_proteins_gbk cache_hits_json
"""
import sys
import os
import json
import hashlib
from copy import deepcopy
from Bio import SeqIO

merged_gbk    = sys.argv[1]
cache_3di     = sys.argv[2]
new_gbk_out   = sys.argv[3]
hits_json_out = sys.argv[4]

# Build set of already-cached MD5 hashes
cached = set()
if os.path.exists(cache_3di) and os.path.getsize(cache_3di) > 0:
    for rec in SeqIO.parse(cache_3di, "fasta"):
        cached.add(rec.id)

records_in   = list(SeqIO.parse(merged_gbk, "genbank"))
records_new  = []
cache_hits   = {}       # full_protein_id ("contig:gene") -> md5
seen_in_run  = set()    # deduplicate within this batch (same seq, different phages)

for rec in records_in:
    new_feats = [f for f in rec.features if f.type != "CDS"]
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        aa      = feat.qualifiers.get("translation", [""])[0]
        gene_id = feat.qualifiers.get("ID", [""])[0]
        full_id = f"{rec.id}:{gene_id}"
        h       = hashlib.md5(aa.encode()).hexdigest()
        if h in cached or h in seen_in_run:
            cache_hits[full_id] = h
        else:
            new_feats.append(feat)
            seen_in_run.add(h)

    if any(f.type == "CDS" for f in new_feats):
        nr = deepcopy(rec)
        nr.features = new_feats
        records_new.append(nr)

os.makedirs(os.path.dirname(new_gbk_out), exist_ok=True)
with open(new_gbk_out, "w") as fh:
    SeqIO.write(records_new, fh, "genbank")
with open(hits_json_out, "w") as fh:
    json.dump(cache_hits, fh)

new_n = sum(1 for r in records_new for f in r.features if f.type == "CDS")
print(f"[filter_predict] cache_hits={len(cache_hits)} new_proteins={new_n}")
