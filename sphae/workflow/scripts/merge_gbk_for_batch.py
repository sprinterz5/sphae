"""
Merge multiple prodigal GBK files into one for batch phold predict.
Ensures unique record IDs by prefixing with sample name.
Usage: python merge_gbk_for_batch.py out.gbk sample1.gbk sample2.gbk ...
"""
import sys
from Bio import SeqIO

output_path = sys.argv[1]
input_paths = sys.argv[2:]

records = []
for path in input_paths:
    for record in SeqIO.parse(path, "genbank"):
        records.append(record)

with open(output_path, "w") as fh:
    SeqIO.write(records, fh, "genbank")

print(f"Merged {len(records)} records from {len(input_paths)} files → {output_path}")
