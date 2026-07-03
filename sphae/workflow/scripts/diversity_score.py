"""
Cocktail Diversity Value — how redundant is each phage relative to the
others in the current batch. Computed via canonical k-mer Jaccard
similarity on the genome FASTA (no alignment, no external tool).

For each sample: diversity_score = 1 - (similarity to its closest
relative in the batch). A phage nearly identical to one already in the
cocktail contributes little (score near 0); a phage unlike anything
else in the batch adds real diversity (score near 1). A lone phage in
a batch of one has nothing to be redundant with, so it scores 1.0.

O(n^2) pairwise comparisons — fine for tens to low hundreds of phages
per batch; would need MinHash/sketching to scale to thousands.

Usage: python diversity_score.py <final_annotate_dir> [--k 21] [--out diversity.json]
Expects {final_annotate_dir}/{sample}/{sample}_genome.fasta (sphae's own layout).
"""
import sys
import os
import glob
import json
import argparse

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


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


def score_batch(final_annotate_dir, k=21):
    genome_paths = glob.glob(os.path.join(final_annotate_dir, "*", "*_genome.fasta"))
    samples = {}
    for path in genome_paths:
        sample = os.path.basename(path).replace("_genome.fasta", "")
        seq = read_fasta_seq(path)
        samples[sample] = canonical_kmers(seq, k)

    names = list(samples)
    results = {}
    for name in names:
        others = [o for o in names if o != name]
        best_sim, best_match = 0.0, None
        for other in others:
            sim = jaccard(samples[name], samples[other])
            if best_match is None or sim > best_sim:
                best_sim, best_match = sim, other
        diversity = 1.0 if not others else round(1.0 - best_sim, 4)
        results[name] = {
            "diversity_score": diversity,
            "most_similar_to": best_match,
            "similarity_to_most_similar": round(best_sim, 4) if best_match is not None else None,
            "sole_sample_in_batch": not others,
        }
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("final_annotate_dir")
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--out", default=None, help="default: <final_annotate_dir>/diversity.json")
    args = ap.parse_args()

    results = score_batch(args.final_annotate_dir, k=args.k)
    out_path = args.out or os.path.join(args.final_annotate_dir, "diversity.json")
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)

    print(f"{'sample':<15} {'diversity':>10} {'closest_match':<20} {'similarity':>10}")
    for name, r in sorted(results.items()):
        match = r["most_similar_to"] or "-"
        sim = f"{r['similarity_to_most_similar']:.4f}" if r["similarity_to_most_similar"] is not None else "-"
        print(f"{name:<15} {r['diversity_score']:>10.4f} {match:<20} {sim:>10}")
    print(f"[diversity_score] {len(results)} samples -> {out_path}")


if __name__ == "__main__":
    main()
