# -----------------------------------------------------
# Virus taxonomy
# -----------------------------------------------------
import pandas as pd
import os


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples_df"], sep="\t")


# get current working directory so absolute paths can be used for input/output files
results = config["results"]


# load report
report: "report/workflow.rst"


# load resources folder path
resources = config["resources"]


# -----------------------------------------------------
# Virus taxonomy rules
# -----------------------------------------------------
localrules: virus_taxonomy_analysis

# -----------------------------------------------------
# 01 geNomad's marker based taxonomy
# -----------------------------------------------------
# run vibrant to identify viral contigs
rule genomad_taxonomy:
    input:
        genomad=resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
        contigs=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        results
        + "09_VIRUS_TAXONOMY/01_genomad/votu_representatives_annotate/votu_representatives_taxonomy.tsv",
    params:
        out_dir=results + "09_VIRUS_TAXONOMY/01_genomad/",
        genomad_dir=resources + "genomad/genomad_db",
    conda:
        "../envs/genomad:1.0.1.yml"
    threads: config["virus_identification"]["genomad_threads"]
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/genomad.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    shell:
        """
        # run vibrant for all virus types
        genomad end-to-end \
        --threads {threads} \
        --verbose \
        --min-score 0.0 \
        --cleanup \
        --splits {threads} \
        --enable-score-calibration \
        {input.contigs} \
        {params.out_dir} \
        {params.genomad_dir}
        """

# -----------------------------------------------------
# 02 MMseqs consensus
# -----------------------------------------------------
rule build_mmseqs:
    output:
        resources + "imgvr_6/virus_tax_db/virus_tax_db"
    params:
        mmseqs_dir=resources + "resources/imgvr_6",
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/build_mmseqs.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="1000",
    shell:
        """
        # download database
        cd {params.mmseqs_dir}
        wget https://zenodo.org/record/6574914/files/virus_tax_db.tar.zst?download=1
        mv 'virus_tax_db.tar.zst?download=1' virus_tax_db.tar.zst
        zstd -d virus_tax_db.tar.zst
        tar -xvf virus_tax_db.tar
        """

rule mmseqs2:
    input:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
        db=resources + "imgvr_6/virus_tax_db/virus_tax_db"
    output:
        results + "09_VIRUS_TAXONOMY/02_mmseqs2/taxonomy_lca.tsv",
    params:
        out_dir=results + "09_VIRUS_TAXONOMY/02_mmseqs2/taxonomy",
        extra_args=config["virus_taxonomy"]["mmseqs_arguments"]
    # conda:
    #     "../envs/mmseqs2:14.7e284--pl5321hf1761c0_0.yml"
    container:
        "docker://quay.io/biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/mmseqs2.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="150000",
    threads: config["virus_taxonomy"]["mmseqs_threads"]
    shell:
        """
        mkdir {params.out_dir}
        cd {params.out_dir}

        mmseqs easy-taxonomy \
        {input.viruses} \
        {input.db} \
        {params.out_dir} \
        tmp \
        --threads {threads} \
        -e 1e-5 \
        -s 6 \
        --blacklist "" \
        --tax-lineage 1 \
        {params.extra_args}
        """


# -----------------------------------------------------
# Analyze taxonomy results
# -----------------------------------------------------
# visualize virus host results
rule virus_taxonomy_analysis:
    input:
        genomad=results + "09_VIRUS_TAXONOMY/01_genomad/votu_representatives_annotate/votu_representatives_taxonomy.tsv",
        mmseqs=results + "09_VIRUS_TAXONOMY/02_mmseqs2/taxonomy_lca.tsv",
    output:
        svg=report(
            results + "09_VIRUS_TAXONOMY/virus_taxonomy_analysis.svg",
            category="Step 09: Virus taxonomy",
        ),
        html=results + "09_VIRUS_TAXONOMY/virus_taxonomy_analysis.html",
        report=results + "09_VIRUS_TAXONOMY/virus_taxonomy_report.tsv",
    params:
        min_genomad_agreement=config["virus_taxonomy"]["genomad_min_agreement"],
        min_mmseqs_agreement=config["virus_taxonomy"]["mmseqs_min_agreement"],
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/virus_taxonomy_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/09_virus_taxonomy_analysis.py"