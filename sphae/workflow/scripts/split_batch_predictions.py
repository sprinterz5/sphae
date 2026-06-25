"""
Split batch predict outputs into per-sample files for phold compare.
phold compare expects {sample}_{3di,aa}.fasta and {sample}_prostT5_3di_mean_probabilities.csv.
Args: sample_name gbk_path predict_dir
"""
import sys
import shutil
from Bio import SeqIO

sample_name = sys.argv[1]
gbk_path    = sys.argv[2]
predict_dir = sys.argv[3]

contig_ids = {rec.id for rec in SeqIO.parse(gbk_path, "genbank")}

for tag in ("3di", "aa"):
    src = f"{predict_dir}/batch_{tag}.fasta"
    dst = f"{predict_dir}/{sample_name}_{tag}.fasta"
    recs = [r for r in SeqIO.parse(src, "fasta") if r.id.split(":")[0] in contig_ids]
    with open(dst, "w") as fh:
        SeqIO.write(recs, fh, "fasta")
    print(f"[split_batch] {sample_name} {tag}: {len(recs)} seqs → {dst}")

# phold compare also reads {prefix}_prostT5_3di_mean_probabilities.csv
# For per-sample it should ideally be filtered too, but phold only uses it for
# confidence scores; a copy of the batch file with renamed prefix works.
src_csv = f"{predict_dir}/batch_prostT5_3di_mean_probabilities.csv"
dst_csv = f"{predict_dir}/{sample_name}_prostT5_3di_mean_probabilities.csv"
shutil.copy2(src_csv, dst_csv)
print(f"[split_batch] {sample_name} probabilities CSV → {dst_csv}")
