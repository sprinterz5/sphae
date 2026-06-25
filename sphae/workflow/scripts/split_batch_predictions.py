"""
Split batch predict outputs into per-sample files, merging with embedding cache.
Proteins in cache_hits use cached 3Di; remaining come from batch predict output.
Args: sample_name gbk_path predict_dir cache_3di_fasta cache_hits_json
"""
import sys
import os
import json
import shutil
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

sample_name     = sys.argv[1]
gbk_path        = sys.argv[2]
predict_dir     = sys.argv[3]
cache_3di_path  = sys.argv[4]
cache_hits_path = sys.argv[5]

# Load embedding cache (md5 → 3Di sequence)
cache_3di = {}
if os.path.exists(cache_3di_path) and os.path.getsize(cache_3di_path) > 0:
    for rec in SeqIO.parse(cache_3di_path, "fasta"):
        cache_3di[rec.id] = str(rec.seq)

# Load per-batch cache hit map (full_protein_id → md5)
cache_hits = {}
if os.path.exists(cache_hits_path):
    with open(cache_hits_path) as fh:
        cache_hits = json.load(fh)

# Load batch predict output (may not exist if all proteins were cached)
batch_3di, batch_aa = {}, {}
batch_3di_path = os.path.join(predict_dir, "batch_3di.fasta")
batch_aa_path  = os.path.join(predict_dir, "batch_aa.fasta")
if os.path.exists(batch_3di_path) and os.path.getsize(batch_3di_path) > 0:
    for rec in SeqIO.parse(batch_3di_path, "fasta"):
        batch_3di[rec.id] = str(rec.seq)
if os.path.exists(batch_aa_path) and os.path.getsize(batch_aa_path) > 0:
    for rec in SeqIO.parse(batch_aa_path, "fasta"):
        batch_aa[rec.id] = str(rec.seq)

# Build per-sample 3di + aa records
recs_3di, recs_aa = [], []
n_cache = n_new = 0

for gbk_rec in SeqIO.parse(gbk_path, "genbank"):
    for feat in gbk_rec.features:
        if feat.type != "CDS":
            continue
        gene_id = feat.qualifiers.get("ID", [""])[0]
        aa_seq  = feat.qualifiers.get("translation", [""])[0]
        full_id = f"{gbk_rec.id}:{gene_id}"

        if full_id in cache_hits:
            h       = cache_hits[full_id]
            seq_3di = cache_3di.get(h, "")
            n_cache += 1
        else:
            seq_3di = batch_3di.get(full_id, "")
            aa_seq  = batch_aa.get(full_id, aa_seq)
            n_new   += 1

        if seq_3di:
            recs_3di.append(SeqRecord(Seq(seq_3di), id=full_id, description=""))
            recs_aa.append(SeqRecord(Seq(aa_seq),   id=full_id, description=""))

dst_3di = os.path.join(predict_dir, f"{sample_name}_3di.fasta")
dst_aa  = os.path.join(predict_dir, f"{sample_name}_aa.fasta")
with open(dst_3di, "w") as fh:
    SeqIO.write(recs_3di, fh, "fasta")
with open(dst_aa, "w") as fh:
    SeqIO.write(recs_aa, fh, "fasta")
print(f"[split_batch] {sample_name}: {len(recs_3di)} seqs (cache={n_cache} new={n_new})")

# Probabilities CSV: copy batch file, or create dummy with prob=1.0 for fully-cached runs
src_csv = os.path.join(predict_dir, "batch_prostT5_3di_mean_probabilities.csv")
dst_csv = os.path.join(predict_dir, f"{sample_name}_prostT5_3di_mean_probabilities.csv")
if os.path.exists(src_csv):
    shutil.copy2(src_csv, dst_csv)
else:
    # All proteins came from cache — write dummy CSV (prob=1.0 passes all phold filters)
    with open(dst_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["id", "mean_probabilities"])
        for rec in recs_3di:
            writer.writerow([rec.id, "1.0"])
print(f"[split_batch] {sample_name} probabilities → {dst_csv}")
