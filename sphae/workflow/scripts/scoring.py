"""
Phage prioritisation scoring — parses sphae final-annotate outputs
(summary.txt + summary.functions) and computes Safety Gate, Genomic
Safety Score, Annotation Confidence, and Genome Quality (if checkv ran).

Formula (weights per component, out of components currently available):
  Safety Gate x (0.50 Genomic Safety + 0.15 Genome Quality
                 + 0.10 Annotation Confidence
                 + 0.10 Host Evidence + 0.10 Cocktail Diversity
                 + 0.05 Clinical/Literature)

Host Evidence, Cocktail Diversity, Clinical/Literature are not computed
here (need NCBI metadata / multi-phage batch / literature DB) and are
reported as null. total_partial is the gate-weighted sum of components
that ARE available, with weight_used shown so it's not mistaken for a
final score.

Args: final_annotate_dir [sample ...]  (default: all subdirs)
"""
import sys
import os
import re
import json
import glob

CORE_PROTEIN_KEYWORDS = {
    "terminase":   ["terminase"],
    "portal":      ["portal"],
    "capsid_head": ["capsid", "major head", "head protein"],
    "tail":        ["tail"],
    "lysis":       ["holin", "lysin", "endolysin", "spanin"],
}


def parse_summary_txt(path):
    data = {
        "cds_count": None,
        "hypothetical_pct": None,
        "genome_length": None,
        "gc_percent": None,
        "n_contigs": None,
        "integrase": False,
        "recombinase": False,
        "transposase": False,
        "lysogeny": False,
        "amr": False,
        "virulence": False,
        "anti_crispr": False,
        "defense": False,
        "completeness": None,
        "contamination": None,
        "checkv_quality": None,
    }
    with open(path) as fh:
        text = fh.read()

    m = re.search(r"Number of contigs:\s*(\d+)", text)
    if m: data["n_contigs"] = int(m.group(1))
    m = re.search(r"Genome length:\s*(\d+)", text)
    if m: data["genome_length"] = int(m.group(1))
    m = re.search(r"GC percent:\s*([\d.]+)", text)
    if m: data["gc_percent"] = float(m.group(1))
    m = re.search(r"Number of CDS:\s*(\d+)", text)
    if m: data["cds_count"] = int(m.group(1))
    m = re.search(r"hypothetical protein':\s*\d+\s*\(([\d.]+)%\)", text)
    if m: data["hypothetical_pct"] = float(m.group(1))
    m = re.search(r"Completeness:\s*([\d.]+)%", text)
    if m: data["completeness"] = float(m.group(1))
    m = re.search(r"Contamination:\s*([\d.]+)%", text)
    if m: data["contamination"] = float(m.group(1))
    m = re.search(r"CheckV quality:\s*(\S+)", text)
    if m: data["checkv_quality"] = m.group(1)

    data["integrase"]    = bool(re.search(r"_CDS_\d+:.*integrase", text, re.I))
    data["recombinase"]  = "No recombinase" not in text
    data["transposase"]  = "No transposase" not in text
    data["lysogeny"]     = "Lysogeny markers found" in text
    data["amr"]           = "No AMR genes found" not in text
    data["virulence"]     = "No virulence factor genes" not in text
    data["anti_crispr"]   = "No anti-CRISPR genes" not in text
    data["defense"]       = "No Defense genes found" not in text

    return data


def parse_core_proteins(functions_path):
    found = {k: False for k in CORE_PROTEIN_KEYWORDS}
    if not os.path.exists(functions_path) or os.path.getsize(functions_path) == 0:
        return found
    with open(functions_path) as fh:
        text = fh.read().lower()
    for category, keywords in CORE_PROTEIN_KEYWORDS.items():
        if any(kw in text for kw in keywords):
            found[category] = True
    return found


def safety_gate(d):
    if d["amr"] or d["virulence"]:
        return 0.0, "AMR or virulence factor detected"
    if d["integrase"] or d["lysogeny"]:
        return 0.5, "integrase/lysogeny marker — manual review"
    return 1.0, "no major red flags"


def genomic_safety_score(d):
    score = 1.0
    reasons = []
    if d["amr"]:
        score -= 0.40; reasons.append("AMR gene (-0.40)")
    if d["virulence"]:
        score -= 0.40; reasons.append("virulence factor (-0.40)")
    if d["integrase"]:
        score -= 0.30; reasons.append("integrase (-0.30)")
    if d["recombinase"] or d["transposase"]:
        score -= 0.15; reasons.append("recombinase/transposase (-0.15)")
    if d["lysogeny"]:
        score -= 0.25; reasons.append("lysogeny/repressor marker (-0.25)")
    return max(0.0, min(1.0, score)), reasons


def genome_quality_score(d):
    if d["completeness"] is None and d["checkv_quality"] is None:
        return None, ["checkv data not available"]
    score = 1.0
    reasons = []
    if d["completeness"] is None:
        score -= 0.10; reasons.append("completeness unknown (-0.10)")
    elif d["completeness"] < 50:
        score -= 0.35; reasons.append("completeness < 50% (-0.35)")
    elif d["completeness"] < 90:
        score -= 0.15; reasons.append("completeness < 90% (-0.15)")
    if d["contamination"] is not None:
        if d["contamination"] > 10:
            score -= 0.35; reasons.append("contamination > 10% (-0.35)")
        elif d["contamination"] > 5:
            score -= 0.20; reasons.append("contamination > 5% (-0.20)")
    if d["n_contigs"] and d["n_contigs"] > 1:
        score -= 0.20; reasons.append("multiple fragmented contigs (-0.20)")
    if d["genome_length"] and d["genome_length"] < 5000:
        score -= 0.15; reasons.append("genome suspiciously small (-0.15)")
    return max(0.0, min(1.0, score)), reasons


def annotation_confidence_score(d, core_proteins):
    score = 1.0
    reasons = []
    hyp = d["hypothetical_pct"]
    if hyp is not None:
        if hyp > 85:
            score -= 0.35; reasons.append(f"hypothetical {hyp}% > 85% (-0.35)")
        elif hyp > 75:
            score -= 0.20; reasons.append(f"hypothetical {hyp}% > 75% (-0.20)")
        elif hyp > 60:
            score -= 0.10; reasons.append(f"hypothetical {hyp}% > 60% (-0.10)")
    n_core_found = sum(core_proteins.values())
    if n_core_found == 0:
        score -= 0.15; reasons.append("no core phage proteins detected (-0.15)")
    elif n_core_found < 3:
        score -= 0.05; reasons.append(f"only {n_core_found}/5 core protein categories found (-0.05)")
    return max(0.0, min(1.0, score)), reasons


def score_sample(sample_dir, sample_name, host_metadata=None):
    summary_txt = os.path.join(sample_dir, f"{sample_name}_summary.txt")
    functions_f = os.path.join(sample_dir, f"{sample_name}_summary.functions")

    if not os.path.exists(summary_txt):
        return {"sample": sample_name, "error": f"missing {summary_txt}"}

    d = parse_summary_txt(summary_txt)
    core_proteins = parse_core_proteins(functions_f)

    gate, gate_reason = safety_gate(d)
    safety, safety_reasons = genomic_safety_score(d)
    quality, quality_reasons = genome_quality_score(d)
    annot, annot_reasons = annotation_confidence_score(d, core_proteins)

    host_entry = (host_metadata or {}).get(sample_name)
    host_score = host_entry["host_evidence_score"] if host_entry else None
    host_reason = host_entry["host_evidence_reason"] if host_entry else "no host_metadata.json entry for this sample"

    weights = {
        "genomic_safety": 0.50, "genome_quality": 0.15,
        "annotation_confidence": 0.10, "host_evidence": 0.10,
    }
    components = {
        "genomic_safety": safety,
        "genome_quality": quality,
        "annotation_confidence": annot,
        "host_evidence": host_score,
    }
    weight_used = sum(w for k, w in weights.items() if components[k] is not None)
    weighted_sum = sum(weights[k] * components[k] for k in weights if components[k] is not None)
    partial_score = gate * weighted_sum if weight_used > 0 else 0.0

    return {
        "sample": sample_name,
        "safety_gate": {"value": gate, "reason": gate_reason},
        "genomic_safety_score": {"value": safety, "reasons": safety_reasons},
        "genome_quality_score": {"value": quality, "reasons": quality_reasons},
        "annotation_confidence_score": {"value": annot, "reasons": annot_reasons},
        "host_evidence_score": {"value": host_score, "reason": host_reason},
        "cocktail_diversity_score": None,
        "clinical_literature_score": None,
        "core_proteins_detected": core_proteins,
        "raw_fields": d,
        "partial_score": round(partial_score, 4),
        "weight_used": weight_used,
        "weight_total": 0.50 + 0.15 + 0.10 + 0.10 + 0.10 + 0.05,
        "note": "partial_score uses only available components (weight_used/weight_total shown); "
                "cocktail_diversity and clinical_literature not yet computed",
    }


def main():
    final_annotate_dir = sys.argv[1]
    args = sys.argv[2:]

    host_metadata = None
    if "--host-metadata" in args:
        i = args.index("--host-metadata")
        with open(args[i + 1]) as fh:
            host_metadata = json.load(fh)
        del args[i:i + 2]

    samples = args
    if not samples:
        samples = sorted(
            d for d in os.listdir(final_annotate_dir)
            if os.path.isdir(os.path.join(final_annotate_dir, d))
        )

    results = []
    for sample in samples:
        sample_dir = os.path.join(final_annotate_dir, sample)
        result = score_sample(sample_dir, sample, host_metadata=host_metadata)
        results.append(result)
        out_json = os.path.join(sample_dir, f"{sample}_score.json")
        with open(out_json, "w") as fh:
            json.dump(result, fh, indent=2)

    print(f"{'sample':<15} {'gate':>5} {'safety':>7} {'quality':>8} {'annot':>6} {'host':>6} {'partial':>8} {'weight_used':>12}")
    for r in results:
        if "error" in r:
            print(f"{r['sample']:<15} ERROR: {r['error']}")
            continue
        q = r["genome_quality_score"]["value"]
        q_str = f"{q:.2f}" if q is not None else "N/A"
        h = r["host_evidence_score"]["value"]
        h_str = f"{h:.2f}" if h is not None else "N/A"
        print(
            f"{r['sample']:<15} "
            f"{r['safety_gate']['value']:>5.2f} "
            f"{r['genomic_safety_score']['value']:>7.2f} "
            f"{q_str:>8} "
            f"{r['annotation_confidence_score']['value']:>6.2f} "
            f"{h_str:>6} "
            f"{r['partial_score']:>8.4f} "
            f"{r['weight_used']:>10.2f}/{r['weight_total']:.2f}"
        )


if __name__ == "__main__":
    main()
