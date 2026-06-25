"""
ORF prediction via pyrodigal-gv Python API.
Outputs phold-compatible GBK (with /ID qualifiers) + FAA.
Replaces: prodigal-gv CLI + prodigal_to_phold_gbk.py
"""
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

try:
    import pyrodigal_gv
    finder = pyrodigal_gv.ViralGeneFinder(meta=True)
except ImportError:
    import pyrodigal
    finder = pyrodigal.GeneFinder(meta=True)

fasta_in, gbk_out, faa_out = sys.argv[1], sys.argv[2], sys.argv[3]

records_in = list(SeqIO.parse(fasta_in, "fasta"))
records_out = []
all_proteins = []

for record in records_in:
    genes = finder.find_genes(bytes(record.seq))

    out_record = SeqRecord(
        record.seq,
        id=record.id,
        name=record.id[:16],
        description=record.description.replace(":", " "),
        annotations={"molecule_type": "DNA"},
    )

    for i, gene in enumerate(genes, 1):
        start  = gene.begin - 1
        end    = gene.end
        strand = 1 if gene.strand == 1 else -1
        gene_id = f"{record.id}_CDS_{i:04d}"
        aa_seq  = gene.translate().rstrip("*")

        feature = SeqFeature(
            FeatureLocation(start, end, strand=strand),
            type="CDS",
            qualifiers={
                "ID":          [gene_id],
                "locus_tag":   [gene_id],
                "product":     ["hypothetical protein"],
                "translation": [aa_seq],
                "phrog":       ["No_PHROGs"],
                "function":    ["unknown function"],
            },
        )
        out_record.features.append(feature)
        all_proteins.append((gene_id, aa_seq))

    records_out.append(out_record)

with open(gbk_out, "w") as fh:
    SeqIO.write(records_out, fh, "genbank")

with open(faa_out, "w") as fh:
    for gene_id, aa_seq in all_proteins:
        fh.write(f">{gene_id}\n{aa_seq}\n")

print(f"[run_prodigal_gv] {len(all_proteins)} ORFs → {gbk_out}, {faa_out}")
