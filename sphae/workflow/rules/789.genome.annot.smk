import re
import shutil
import sys
from pathlib import Path

"""
PATTERNS
"""
PATTERN_LONG = "{sample}.fasta"

"""
RESOLVER FUNCTION
"""
def resolve_input(wc):
    genome_dir = config['args'].get('genome')

    genome = os.path.join(genome_dir, f"{wc.sample}.fasta") if genome_dir else None

    if genome and Path(genome).exists():
        return genome
    raise ValueError(f"No input found for {wc.sample}")


rule prodigal_annotate_genome:
    """ORF prediction via pyrodigal-gv Python API; outputs phold-compatible GBK"""
    input:
        fasta=resolve_input,
    params:
        script=os.path.join(dir_script, "run_prodigal_gv.py"),
    output:
        gbk=os.path.join(dir_annot, "{sample}-prodigal", "{sample}.gbk"),
        faa=os.path.join(dir_annot, "{sample}-prodigal", "{sample}.faa"),
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime'],
    log:
        os.path.join(dir_log, "prodigal.{sample}.log")
    shell:
        """
        mkdir -p $(dirname {output.gbk})
        python {params.script} {input.fasta} {output.gbk} {output.faa} 2> {log}
        """


rule checkv_run_genome:
    """CheckV genome quality assessment; falls back to empty TSV if not installed"""
    input:
        fasta=resolve_input,
    output:
        quality=os.path.join(dir_annot, "{sample}-checkv", "quality_summary.tsv"),
    params:
        outdir=os.path.join(dir_annot, "{sample}-checkv"),
        db=config['args']['checkv_db'],
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime'],
    log:
        os.path.join(dir_log, "checkv.{sample}.log")
    shell:
        """
        mkdir -p {params.outdir}
        if command -v checkv &> /dev/null; then
            checkv end_to_end {input.fasta} {params.outdir} -t {threads} -d {params.db} 2> {log}
        else
            echo "checkv not found, writing empty quality summary" > {log}
            printf "contig_id\tcheckv_quality\tcompleteness\tcontamination\n" > {output.quality}
        fi
        """


rule phold_predict_batch:
    """
    Run ProstT5 embeddings for ALL samples in ONE GPU pass.
    Model loads once; proteins from all phages are embedded together.
    phold compare (per-sample) reuses the shared predictions_dir.
    """
    input:
        gbks=expand(
            os.path.join(dir_annot, "{sample}-prodigal", "{sample}.gbk"),
            sample=samples_names
        ),
    params:
        merged=BATCH_PREDICT_GBK,
        outdir=BATCH_PREDICT_DIR,
        prefix="batch",
        db=config['args']['phold_db'],
        cpu=PHOLD_CPU_FLAG,
        batch_size=config['params'].get('phold_batch_size', 32),
        script=os.path.join(dir_script, "merge_gbk_for_batch.py"),
    output:
        sentinel=BATCH_PREDICT_SENTINEL,
    conda:
        os.path.join(dir_env, "phold.yaml")
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb=config['resources']['bigjob']['mem_mb'],
        runtime=config['resources']['bigjob']['runtime'],
        gpu=PHOLD_GPU_RESOURCE,
    log:
        os.path.join(dir_log, "phold_predict_batch.log")
    shell:
        """
        mkdir -p $(dirname {params.merged})
        python {params.script} {params.merged} {input.gbks}
        if [[ -s {params.merged} ]] ; then
            TRANSFORMERS_OFFLINE=1 HF_DATASETS_OFFLINE=1 \
            phold predict \
                -i {params.merged} \
                -o {params.outdir} \
                -p {params.prefix} \
                -t {threads} \
                {params.cpu} \
                -d {params.db} \
                -f \
                --batch_size {params.batch_size} \
                2> {log}
        fi
        touch {output.sentinel}
        """


rule phold_split_predictions:
    """Split batch 3Di/AA fastas into per-sample files that phold compare expects."""
    input:
        sentinel=BATCH_PREDICT_SENTINEL,
        gbk=os.path.join(dir_annot, "{sample}-prodigal", "{sample}.gbk"),
    output:
        fasta_3di=os.path.join(BATCH_PREDICT_DIR, "{sample}_3di.fasta"),
        fasta_aa=os.path.join(BATCH_PREDICT_DIR, "{sample}_aa.fasta"),
        prob_csv=os.path.join(BATCH_PREDICT_DIR, "{sample}_prostT5_3di_mean_probabilities.csv"),
    params:
        script=os.path.join(dir_script, "split_batch_predictions.py"),
        predict_dir=BATCH_PREDICT_DIR,
    conda:
        os.path.join(dir_env, "pharokka.yaml")
    threads: 1
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=30,
    log:
        os.path.join(dir_log, "phold_split.{sample}.log")
    shell:
        """
        python {params.script} {wildcards.sample} {input.gbk} {params.predict_dir} 2> {log}
        """


rule phold_compare_genome:
    """
    Per-sample Foldseek structural search using the shared batch predictions dir.
    Runs in parallel across samples (CPU-bound, no GPU needed here).
    """
    input:
        gbk=os.path.join(dir_annot, "{sample}-prodigal", "{sample}.gbk"),
        fasta_3di=os.path.join(BATCH_PREDICT_DIR, "{sample}_3di.fasta"),
        prob_csv=os.path.join(BATCH_PREDICT_DIR, "{sample}_prostT5_3di_mean_probabilities.csv"),
    params:
        predict=BATCH_PREDICT_DIR,
        o=os.path.join(dir_annot, "{sample}-phold"),
        prefix="{sample}",
        db=config['args']['phold_db'],
        foldseek_gpu=PHOLD_FOLDSEEK_GPU_FLAG,
    output:
        gbk=os.path.join(dir_annot, "{sample}-phold", "{sample}.gbk"),
        acr=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"),
    threads:
        config['resources']['smalljob']['threads']
    conda:
        os.path.join(dir_env, "phold.yaml")
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime'],
    log:
        os.path.join(dir_log, "phold_compare.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            if phold compare \
                    -i {input.gbk} \
                    --predictions_dir {params.predict} \
                    -p {params.prefix} \
                    -o {params.o} \
                    -t {threads} \
                    -d {params.db} \
                    -f \
                    {params.foldseek_gpu} \
                    2> {log} ; then
                echo "phold compare OK" >> {log}
            else
                echo "phold compare with foldseek_gpu failed, retrying without GPU" >> {log}
                phold compare \
                    -i {input.gbk} \
                    --predictions_dir {params.predict} \
                    -p {params.prefix} \
                    -o {params.o} \
                    -t {threads} \
                    -d {params.db} \
                    -f \
                    2>> {log}
            fi
        else
            touch {output.gbk}
            touch {output.acr}
            touch {output.card}
            touch {output.defense}
            touch {output.vfdb}
        fi
        """


rule phynteny_run_genome:
    input:
        gbk=os.path.join(dir_annot, "{sample}-phold", "{sample}.gbk")
    params:
        odir=os.path.join(dir_annot, "{sample}-phynteny"),
        model=config['args']['phynteny_db'],
    output:
        pkl=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk")
    conda:
        os.path.join(dir_env, "phynteny.yaml")
    threads:
        config['resources']['smalljob']['threads']
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime'],
        gpu=PHOLD_GPU_RESOURCE,
    log:
        os.path.join(dir_log, "phynteny.{sample}.log")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            phynteny_transformer {input.gbk} -o {params.odir} \
                -m {params.model} -f \
                2> {log}
            touch {output.pkl}
        else
            touch {output.pkl}
        fi
        """


rule phynteny_plotter_genomes:
    input:
        gbk=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk"),
        fasta=resolve_input
    params:
        gff3=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gff3"),
        prefix="{sample}",
        output=os.path.join(dir_annot, "{sample}-phynteny", "plots")
    output:
        plot=directory(os.path.join(dir_annot, "{sample}-phynteny", "plots"))
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime'],
    conda:
        os.path.join(dir_env, "phold.yaml")
    shell:
        """
        if [[ -s {input.gbk} ]] ; then
            genbank_to -g {input.gbk} --gff3 {params.gff3}
            phold plot -i {input.gbk} -f -p {params.prefix} -o {params.output}
        fi
        """


rule summarize_annotations_genome:
    input:
        phold=os.path.join(dir_annot, "{sample}-phold", "{sample}.gbk"),
        pkl=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk"),
    output:
        phold_func=os.path.join(dir_annot, "{sample}-phold", "{sample}_phold.functions"),
        pkl_func=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.functions"),
    resources:
        mem_mb=config['resources']['smalljob']['mem_mb'],
        runtime=config['resources']['smalljob']['runtime'],
    conda:
        os.path.join(dir_env, "phold.yaml")
    shell:
        """
        if [[ -s {input.phold} ]] ; then
            genbank_to -g {input.phold} -f {output.phold_func}
        else
            touch {output.phold_func}
        fi

        if [[ -s {input.pkl} ]] ; then
            genbank_to -g {input.pkl} -f {output.pkl_func}
        else
            touch {output.pkl_func}
        fi
        """


rule annotate_summary_genome:
    input:
        phold_func=os.path.join(dir_annot, "{sample}-phold", "{sample}_phold.functions"),
        pkl_func=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.functions"),
    output:
        summary_gbk=os.path.join(dir_final, "{sample}", "{sample}_summary.functions")
    params:
        tmp=os.path.join(dir_annot, "{sample}-phynteny", "temp")
    localrule: True
    script:
        os.path.join(dir_script, "summary_annot_functions.py")


rule summarize:
    input:
        genome=resolve_input,
        gbk=os.path.join(dir_annot, "{sample}-phynteny", "phynteny.gbk"),
        plot=os.path.join(dir_annot, "{sample}-phynteny", "plots"),
        checkv=os.path.join(dir_annot, "{sample}-checkv", "quality_summary.tsv"),
        acr=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        summary=os.path.join(dir_final, "{sample}", "{sample}_summary.txt")
    params:
        genomes=os.path.join(dir_final, "{sample}", "{sample}_genome.fasta"),
        gbks=os.path.join(dir_final, "{sample}", "{sample}.gbk"),
        plots=directory(os.path.join(dir_final, "{sample}", "{sample}_phynteny")),
        outdir=os.path.join(dir_final),
        sample="{sample}",
    localrule: True
    script:
        os.path.join(dir_script, 'summary-annot.py')


rule accessory_files_genome:
    input:
        summary=os.path.join(dir_final, "{sample}", "{sample}_summary.txt"),
        acr=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "acr_cds_predictions.tsv"),
        card=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "card_cds_predictions.tsv"),
        defense=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "defensefinder_cds_predictions.tsv"),
        vfdb_phold=os.path.join(dir_annot, "{sample}-phold", "sub_db_tophits", "vfdb_cds_predictions.tsv"),
    output:
        amr=os.path.join(dir_final, "{sample}", "{sample}_phold_amr.tsv"),
        vfdb=os.path.join(dir_final, "{sample}", "{sample}_phold_vfdb.tsv"),
        acr=os.path.join(dir_final, "{sample}", "{sample}_phold_acr.tsv"),
        defense=os.path.join(dir_final, "{sample}", "{sample}_phold_defense.tsv"),
    params:
        sample="{sample}"
    localrule: True
    shell:
        """
        cp {input.acr} {output.acr}
        cp {input.card} {output.amr}
        cp {input.defense} {output.defense}
        cp {input.vfdb_phold} {output.vfdb}
        """
