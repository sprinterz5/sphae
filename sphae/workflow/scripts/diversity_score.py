"""
Cocktail Diversity Value — how redundant is each phage relative to
every other phage sphae has ever scored, not just the current batch.
Computed via canonical k-mer Jaccard similarity on the genome FASTA
(no alignment, no external tool).

Diversity is inherently a *relative* measure (unlike the other scoring
components, which are intrinsic properties of a single genome), so it
needs a persistent, cross-run library — comparing only within the
current --output directory would silently miss a near-duplicate that
was scored in a separate batch/output dir. Per-genome canonical k-mer
sets are cached under {library_dir}/genome_kmers/{sample}.kmers and
reused/extended on every call, independent of --output.

For each sample: diversity_score = 1 - (similarity to its closest
relative across the whole library + current batch). A phage nearly
identical to one already scored contributes little (score near 0); a
phage unlike anything seen so far adds real diversity (score near 1).
A phage with nothing to compare against yet (empty library, batch of
one) scores 1.0 — that's an "unknown", not a confirmed "unique".

O(n^2) pairwise comparisons against the whole library — fine for tens
to low hundreds of phages; would need MinHash/sketching to scale to
thousands.

Usage: python diversity_score.py <final_annotate_dir> [--k 21] [--out diversity.json]
       [--library-dir ~/.sphae/library] [--no-library]
Expects {final_annotate_dir}/{sample}/{sample}_genome.fasta (sphae's own layout).
"""
import sys
import os
import glob
import json
import argparse

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")
DEFAULT_LIBRARY_DIR = os.path.join(os.path.expanduser("~"), ".sphae", "library")


def read_fasta_seq(path):
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq_parts.append(line.strip())
    return "".join(seq_parts).upper()


def canonical_kmers(seq, k):
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if "N" in kmer:
            continue
        rc = kmer.translate(COMPLEMENT)[::-1]
        kmers.add(min(kmer, rc))
    return kmers


def jaccard(a, b):
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def kmer_dir(library_dir):
    return os.path.join(library_dir, "genome_kmers")


def load_library(library_dir):
    samples = {}
    d = kmer_dir(library_dir)
    if os.path.isdir(d):
        for path in glob.glob(os.path.join(d, "*.kmers")):
            sample = os.path.basename(path)[:-len(".kmers")]
            with open(path) as fh:
                samples[sample] = set(line.strip() for line in fh if line.strip())
    return samples


def save_to_library(library_dir, sample, kmers):
    d = kmer_dir(library_dir)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{sample}.kmers"), "w") as fh:
        fh.write("\n".join(sorted(kmers)))
        fh.write("\n")


def score_batch(final_annotate_dir, k=21, library_dir=DEFAULT_LIBRARY_DIR):
    genome_paths = glob.glob(os.path.join(final_annotate_dir, "*", "*_genome.fasta"))
    batch_samples = {}
    for path in genome_paths:
        sample = os.path.basename(path).replace("_genome.fasta", "")
        seq = read_fasta_seq(path)
        batch_samples[sample] = canonical_kmers(seq, k)

    library_samples = load_library(library_dir) if library_dir else {}
    # current batch wins on name collision (e.g. a sample re-annotated with a new genome)
    combined = {**library_samples, **batch_samples}

    results = {}
    for name in batch_samples:
        others = [o for o in combined if o != name]
        best_sim, best_match = 0.0, None
        for other in others:
            sim = jaccard(combined[name], combined[other])
            if best_match is None or sim > best_sim:
                best_sim, best_match = sim, other
        diversity = 1.0 if not others else round(1.0 - best_sim, 4)
        results[name] = {
            "diversity_score": diversity,
            "most_similar_to": best_match,
            "similarity_to_most_similar": round(best_sim, 4) if best_match is not None else None,
            "sole_sample_in_batch": not others,
            "compared_against_library_size": len(others),
        }

    if library_dir:
        for name, kmers in batch_samples.items():
            save_to_library(library_dir, name, kmers)

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("final_annotate_dir")
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--out", default=None, help="default: <final_annotate_dir>/diversity.json")
    ap.add_argument("--library-dir", default=DEFAULT_LIBRARY_DIR,
                     help="persistent cross-run k-mer library (default: ~/.sphae/library)")
    ap.add_argument("--no-library", action="store_true",
                     help="compare only within this batch, don't read/write the persistent library")
    args = ap.parse_args()

    library_dir = None if args.no_library else args.library_dir
    results = score_batch(args.final_annotate_dir, k=args.k, library_dir=library_dir)
    out_path = args.out or os.path.join(args.final_annotate_dir, "diversity.json")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)

    print(f"{'sample':<15} {'diversity':>10} {'closest_match':<20} {'similarity':>10} {'vs_n':>5}")
    for name, r in sorted(results.items()):
        match = r["most_similar_to"] or "-"
        sim = f"{r['similarity_to_most_similar']:.4f}" if r["similarity_to_most_similar"] is not None else "-"
        print(f"{name:<15} {r['diversity_score']:>10.4f} {match:<20} {sim:>10} {r['compared_against_library_size']:>5}")
    print(f"[diversity_score] {len(results)} samples -> {out_path} (library: {library_dir or 'disabled'})")


if __name__ == "__main__":
    main()
