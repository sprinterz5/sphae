"""
Download phage genomes from NCBI as GenBank (not FASTA) so the /host
qualifier is available for Host Evidence scoring.

Usage:
  python ncbi_downloader.py --query "Salmonella phage" --outdir genomes/ --max 20
  python ncbi_downloader.py --accessions NC_001416.1,NC_002371.2 --outdir genomes/

Writes, per accession, into outdir:
  {accession}.fasta        — genome sequence (feeds sphae annotate --genome)
  {accession}.gbk           — full GenBank record (kept for provenance)
Writes once, for the whole batch:
  host_metadata.json        — {accession: {organism, host, lab_host, strain, host_evidence_tier}}

NCBI eutils rate limit without an API key is 3 req/sec; requests are
batched (comma-joined ids) and throttled accordingly.
"""
import sys
import os
import re
import json
import time
import argparse
import urllib.request
import urllib.parse

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH_SIZE = 20
REQUEST_DELAY = 0.4  # ~2.5 req/sec, safely under the 3 req/sec no-key limit


def eutils_get(endpoint, params):
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read().decode()


def esearch_ids(query, max_results):
    xml = eutils_get("esearch.fcgi", {
        "db": "nuccore", "term": query, "retmax": max_results, "retmode": "xml",
    })
    return re.findall(r"<Id>(\d+)</Id>", xml)


def efetch_genbank(ids):
    return eutils_get("efetch.fcgi", {
        "db": "nuccore", "id": ",".join(ids), "rettype": "gb", "retmode": "text",
    })


def split_genbank_records(text):
    """A multi-record efetch response is '//'-terminated GenBank flat files back to back."""
    records = [r.strip() for r in text.split("//\n") if r.strip()]
    return records


def parse_accession(record_text):
    m = re.search(r"^VERSION\s+(\S+)", record_text, re.M)
    return m.group(1) if m else None


def parse_source_qualifiers(record_text):
    m = re.search(r"^     source\s+.*?(?=\n     \S|\nORIGIN)", record_text, re.M | re.S)
    block = m.group(0) if m else ""
    fields = {}
    for key in ("organism", "host", "lab_host", "strain"):
        qm = re.search(rf'/{key}="([^"]+)"', block)
        if qm:
            # qualifier values can wrap across lines in the flat file; collapse whitespace
            fields[key] = re.sub(r"\s+", " ", qm.group(1)).strip()
    return fields


def parse_origin_sequence(record_text):
    m = re.search(r"^ORIGIN.*?\n(.*?)(?=^//|\Z)", record_text, re.M | re.S)
    if not m:
        return ""
    bases = re.findall(r"[acgtnACGTN]+", m.group(1))
    return "".join(bases).upper()


def host_evidence_tier(fields):
    """0.50 if NCBI /host qualifier present, 0.30 if only inferable from organism name, 0.10 unknown."""
    if fields.get("host") or fields.get("lab_host"):
        return 0.50, "NCBI /host or /lab_host qualifier present"
    organism = fields.get("organism", "")
    if organism and len(organism.split()) >= 2:
        # e.g. "Salmonella phage P22" -> genus "Salmonella" inferred from title
        return 0.30, f"host genus inferred from organism name '{organism}'"
    return 0.10, "no host information found"


def download(accessions=None, query=None, max_results=20, outdir="genomes"):
    os.makedirs(outdir, exist_ok=True)

    if query and not accessions:
        ids = esearch_ids(query, max_results)
        if not ids:
            print(f"[ncbi_downloader] no results for query: {query}")
            return
    elif accessions:
        ids = accessions
    else:
        raise ValueError("must provide --query or --accessions")

    metadata = {}
    for i in range(0, len(ids), BATCH_SIZE):
        batch = ids[i:i + BATCH_SIZE]
        text = efetch_genbank(batch)
        for record in split_genbank_records(text):
            accession = parse_accession(record)
            if not accession:
                continue
            fields = parse_source_qualifiers(record)
            seq = parse_origin_sequence(record)
            tier, reason = host_evidence_tier(fields)

            fasta_path = os.path.join(outdir, f"{accession}.fasta")
            gbk_path = os.path.join(outdir, f"{accession}.gbk")
            organism = fields.get("organism", accession)
            with open(fasta_path, "w") as fh:
                fh.write(f">{accession} {organism}\n")
                for j in range(0, len(seq), 70):
                    fh.write(seq[j:j + 70] + "\n")
            with open(gbk_path, "w") as fh:
                fh.write(record + "\n//\n")

            metadata[accession] = {
                **fields,
                "host_evidence_score": tier,
                "host_evidence_reason": reason,
            }
            print(f"[ncbi_downloader] {accession}: {organism} "
                  f"(host={fields.get('host', 'N/A')}, tier={tier})")
        time.sleep(REQUEST_DELAY)

    with open(os.path.join(outdir, "host_metadata.json"), "w") as fh:
        json.dump(metadata, fh, indent=2)
    print(f"[ncbi_downloader] {len(metadata)} genomes -> {outdir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", help="NCBI nuccore search term, e.g. 'Salmonella phage'")
    ap.add_argument("--accessions", help="comma-separated accession list")
    ap.add_argument("--max", type=int, default=20, help="max results for --query")
    ap.add_argument("--outdir", default="genomes")
    args = ap.parse_args()

    accessions = args.accessions.split(",") if args.accessions else None
    download(accessions=accessions, query=args.query, max_results=args.max, outdir=args.outdir)


if __name__ == "__main__":
    main()
