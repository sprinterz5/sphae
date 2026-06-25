"""
Append new phold predict results to the persistent embedding cache.
Cache key: MD5(original_gbk_aa) — matches filter_for_predict.py hash computation.
Phold masks residues in batch_aa.fasta (e.g. G->X), so we hash from the GBK instead.
Args: batch_3di_fasta new_proteins_gbk cache_3di_out cache_aa_out
"""
import sys
import os
import hashlib
from Bio import SeqIO

batch_3di     = sys.argv[1]
new_prot_gbk  = sys.argv[2]
cache_3di     = sys.argv[3]
cache_aa      = sys.argv[4]

if not os.path.exists(batch_3di) or os.path.getsize(batch_3di) == 0:
    print("[update_cache] batch_3di.fasta empty/missing — nothing to cache")
    sys.exit(0)

# Build full_id -> original_aa from the GBK (same source as filter_for_predict.py)
gbk_aa = {}
for rec in SeqIO.parse(new_prot_gbk, "genbank"):
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        gene_id = feat.qualifiers.get("ID", [""])[0]
        full_id = f"{rec.id}:{gene_id}"
        gbk_aa[full_id] = feat.qualifiers.get("translation", [""])[0]

os.makedirs(os.path.dirname(cache_3di) if os.path.dirname(cache_3di) else ".", exist_ok=True)

added = 0
with open(cache_3di, "a") as f3, open(cache_aa, "a") as fa:
    for rec in SeqIO.parse(batch_3di, "fasta"):
        aa = gbk_aa.get(rec.id, "")
        if not aa:
            continue
        h = hashlib.md5(aa.encode()).hexdigest()
        f3.write(f">{h}\n{rec.seq}\n")
        fa.write(f">{h}\n{aa}\n")
        added += 1

print(f"[update_cache] +{added} entries → {cache_3di}")
