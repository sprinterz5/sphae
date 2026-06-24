"""
Add /ID qualifiers to prodigal-gv GBK output so phold can identify CDS features
(phold looks for pharokka-style /ID qualifiers to track proteins).
"""
import sys
from Bio import SeqIO

input_gbk, output_gbk = sys.argv[1], sys.argv[2]

records = list(SeqIO.parse(input_gbk, "genbank"))
for record in records:
    cds_idx = 0
    for feature in record.features:
        if feature.type == "CDS":
            cds_idx += 1
            feature.qualifiers["ID"] = [f"{record.id}_CDS_{cds_idx:04d}"]

with open(output_gbk, "w") as fh:
    SeqIO.write(records, fh, "genbank")
