from pathlib import Path
import shutil
import pandas as pd
from Bio import SeqIO
import glob, os


def compute_gc(seq):
    seq = str(seq).upper()
    gc = seq.count('G') + seq.count('C')
    return round(gc / len(seq) * 100, 2) if len(seq) > 0 else 0


def copy_files(input_files, params):
    shutil.copy(input_files['genome'], params['genomes'])
    shutil.copy(input_files['gbk'], params['gbks'])

    os.makedirs(params['plots'], exist_ok=True)

    plot_files = (glob.glob(os.path.join(input_files['plot'], "*.svg")) +
                  glob.glob(os.path.join(input_files['plot'], "*.png")))

    if not plot_files:
        print("Warning: No plot files found")
    else:
        for f in plot_files:
            shutil.copy(f, params['plots'])


def is_protein(seq):
    return any(c not in "ACGTUNacgtun" for c in seq if c.isalpha())


def generate_summary(input_files, output_summary, params):
    copy_files(input_files, params)

    records = list(SeqIO.parse(input_files['genome'], "fasta"))

    with open(output_summary, 'w') as summary:
        summary.write(f"Sample: {params['sample']}\n")

        if not records:
            summary.write("Empty genome file\n")
            return

        is_prot = is_protein(str(records[0].seq))
        summary.write(f"Sequence type: {'protein' if is_prot else 'nucleotide'}\n")

        if not is_prot:
            summary.write(f"Number of contigs: {len(records)}\n")
            total_len = sum(len(r.seq) for r in records)
            gc = compute_gc(records[0].seq)
            summary.write(f"Genome length: {total_len}\n")
            summary.write(f"GC percent: {gc}\n")

        # --- CDS count and hypothetical proteins from phynteny GBK ---
        cds_count = 0
        hypo_count = 0
        for record in SeqIO.parse(input_files['gbk'], "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    cds_count += 1
                    product = feature.qualifiers.get("product", [""])[0].lower()
                    if "hypothetical" in product:
                        hypo_count += 1
        summary.write(f"Number of CDS: {cds_count}\n")
        hypo_pct = round(hypo_count / cds_count * 100, 1) if cds_count > 0 else 0
        summary.write(f"Total number of CDS annotated as 'hypothetical protein': {hypo_count} ({hypo_pct}%)\n")

        # --- CheckV quality ---
        checkv_path = input_files.get('checkv', '')
        if checkv_path and os.path.exists(checkv_path) and os.path.getsize(checkv_path) > 0:
            try:
                cv = pd.read_csv(checkv_path, sep='\t')
                if not cv.empty:
                    row = cv.iloc[0]
                    summary.write(f"CheckV quality: {row.get('checkv_quality', 'N/A')}\n")
                    summary.write(f"Completeness: {row.get('completeness', 'N/A')}%\n")
                    summary.write(f"Contamination: {row.get('contamination', 'N/A')}%\n")
            except Exception as e:
                summary.write(f"CheckV: parse error ({e})\n")
        else:
            summary.write("CheckV: not available\n")

        # --- Read GBK once for keyword scans ---
        gbk_text = open(input_files['gbk']).read().lower()

        # --- Integrases ---
        found_integrase = False
        for record in SeqIO.parse(input_files['gbk'], 'genbank'):
            for feature in record.features:
                if feature.type == "CDS" and 'product' in feature.qualifiers:
                    product = feature.qualifiers['product'][0].lower()
                    if 'integra' in product:
                        found_integrase = True
                        gene_id = feature.qualifiers.get('locus_tag',
                                  feature.qualifiers.get('ID', ['unknown']))[0]
                        summary.write(f"\t{gene_id}: product=\"{product}\"\n")

        if not found_integrase:
            summary.write("No Integrases\n")
            if 'integra' in gbk_text:
                summary.write("\t...but possible low-confidence hits detected\n")

        # --- Other mobile elements ---
        summary.write("Recombinases found in genome\n" if 'recombinase' in gbk_text else "No recombinase\n")
        summary.write("Transposases found in genome\n" if 'transposase' in gbk_text else "No transposase\n")
        summary.write("Lysogeny markers found\n" if any(x in gbk_text for x in ('lysogen', 'repressor', 'cl protein')) else "No lysogeny markers\n")

        # --- AMR (phold CARD) ---
        try:
            card_lines = open(input_files['card']).readlines()
            summary.write("AMR genes found\n" if len(card_lines) > 0 else "No AMR genes found\n")
        except Exception:
            summary.write("AMR: no data\n")

        # --- Virulence (phold VFDB) ---
        try:
            vf_lines = open(input_files['vfdb_phold']).readlines()
            summary.write("Virulence genes found\n" if len(vf_lines) > 0 else "No virulence factor genes\n")
        except Exception:
            summary.write("Virulence: no data\n")

        # --- Anti-CRISPR (phold ACR) ---
        try:
            acr_lines = open(input_files['acr']).readlines()
            summary.write("Anti-CRISPR genes found\n" if len(acr_lines) > 0 else "No anti-CRISPR genes\n")
        except Exception:
            summary.write("Anti-CRISPR: no data\n")

        # --- Defense systems (phold) ---
        try:
            def_lines = open(input_files['defense']).readlines()
            summary.write("Defense genes found\n" if len(def_lines) > 0 else "No Defense genes found\n")
        except Exception:
            summary.write("Defense: no data\n")


# --- Snakemake integration ---
input_files = {
    'genome':    snakemake.input.genome,
    'gbk':       snakemake.input.gbk,
    'plot':      snakemake.input.plot,
    'checkv':    snakemake.input.checkv,
    'acr':       snakemake.input.acr,
    'card':      snakemake.input.card,
    'defense':   snakemake.input.defense,
    'vfdb_phold':snakemake.input.vfdb_phold,
}

output_summary = snakemake.output.summary

params = {
    'sample': snakemake.params.sample,
    'genomes': snakemake.params.genomes,
    'gbks':    snakemake.params.gbks,
    'plots':   snakemake.params.plots,
}

generate_summary(input_files, output_summary, params)
