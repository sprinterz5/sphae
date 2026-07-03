"""
Host Match — search the persistent NCBI host index for phage candidates
against a query bacterium. This is a separate, on-demand operation from
scoring: it answers "which phages match THIS bacterium", not a fixed
per-phage property, so it's a query/filter over pre-computed metadata
rather than something baked into phage_quality_score.

The index (~/.sphae/library/host_index.json by default) is built
incrementally by ncbi_downloader.py every time genomes are downloaded —
searching it does no network calls and no recomputation, just a filter
over data already gathered at annotation time.

Match tiers (highest first):
  direct_host   — query matches the NCBI /host or /lab_host qualifier
                  (phage was actually isolated on/near this bacterium)
  organism_name — query matches only within the phage's own organism
                  name (e.g. "Salmonella phage P22" for a "Salmonella" query) —
                  weaker signal, inferred rather than observed

This is substring/genus text matching against metadata sphae already
downloaded — NOT experimental host-range prediction (no CRISPR spacer
alignment, no receptor-binding protein analysis). Absence from these
results means "no NCBI-recorded host match", not "confirmed non-host".

Usage: python host_match.py "<bacterium query>" [--library-dir ~/.sphae/library] [--index-file <path>]
"""
import sys
import os
import json
import argparse

DEFAULT_LIBRARY_DIR = os.path.join(os.path.expanduser("~"), ".sphae", "library")


def load_index(library_dir=None, index_file=None):
    if index_file:
        path = index_file
    else:
        path = os.path.join(library_dir or DEFAULT_LIBRARY_DIR, "host_index.json")
    if not os.path.exists(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def search(index, query):
    q = query.strip().lower()
    matches = []
    for accession, fields in index.items():
        host = (fields.get("host") or "").lower()
        lab_host = (fields.get("lab_host") or "").lower()
        organism = (fields.get("organism") or "").lower()

        if q in host or q in lab_host:
            match_type = "direct_host"
        elif q in organism:
            match_type = "organism_name"
        else:
            continue

        matches.append({
            "accession": accession,
            "organism": fields.get("organism"),
            "host": fields.get("host"),
            "lab_host": fields.get("lab_host"),
            "match_type": match_type,
            "host_evidence_score": fields.get("host_evidence_score"),
        })

    rank = {"direct_host": 0, "organism_name": 1}
    matches.sort(key=lambda m: (rank[m["match_type"]], m["accession"]))
    return matches


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("query", help="bacterium name or genus to search for, e.g. 'Pseudomonas aeruginosa'")
    ap.add_argument("--library-dir", default=DEFAULT_LIBRARY_DIR)
    ap.add_argument("--index-file", default=None, help="search a specific host_metadata.json/host_index.json instead of the library")
    ap.add_argument("--out", default=None, help="optional path to write matches as JSON")
    args = ap.parse_args()

    index = load_index(library_dir=args.library_dir, index_file=args.index_file)
    if not index:
        print(f"[host_match] no host index found at {args.index_file or os.path.join(args.library_dir, 'host_index.json')} "
              f"— run ncbi_downloader.py first")
        return

    matches = search(index, args.query)

    print(f"[host_match] query='{args.query}' — {len(matches)} match(es) out of {len(index)} indexed genomes")
    print(f"{'accession':<18} {'match_type':<14} {'host':<30} {'organism':<40}")
    for m in matches:
        print(f"{m['accession']:<18} {m['match_type']:<14} {(m['host'] or '-'):<30} {(m['organism'] or '-'):<40}")

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(matches, fh, indent=2)
        print(f"[host_match] wrote {len(matches)} matches -> {args.out}")


if __name__ == "__main__":
    main()
