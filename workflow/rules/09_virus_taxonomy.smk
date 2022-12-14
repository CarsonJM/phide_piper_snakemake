# -----------------------------------------------------
# Virus Taxonomy Module (if include_taxonomy_module = True)
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
localrules:
    virus_taxonomy_analysis

# -----------------------------------------------------
# 01 geNomad
# -----------------------------------------------------
# run genomad to identify virus taxonomy
rule genomad_taxonomy:
    message:
        "Running geNomad to predict virus taxonomy"
    input:
        genomad=resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
        contigs=results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_representatives.fasta",
    output:
        taxonomy=results
        + "09_VIRUS_TAXONOMY/01_genomad/votu_representatives_annotate/votu_representatives_taxonomy.tsv",
        lifestyle=results
        + "09_VIRUS_TAXONOMY/01_genomad/votu_representatives_find_proviruses/votu_representatives_provirus.tsv",
    params:
        out_dir=results + "09_VIRUS_TAXONOMY/01_genomad/",
        genomad_dir=resources + "genomad/genomad_db",
    conda:
        "../envs/genomad:1.3.0--pyhdfd78af_0.yml"
    threads: config["virus_taxonomy"]["genomad_threads"]
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/genomad_taxonomy.tsv"
    resources:
        runtime=config["virus_taxonomy"]["genomad_runtime"],
        mem_mb=config["virus_taxonomy"]["genomad_memory"],
    shell:
        """
        # run genomad on viruses to predict taxonomy
        genomad end-to-end \
        --threads {threads} \
        --verbose \
        --min-score 0.0 \
        --cleanup \
        --splits {threads} \
        {input.contigs} \
        {params.out_dir} \
        {params.genomad_dir}
        """


# -----------------------------------------------------
# Analyze taxonomy results
# -----------------------------------------------------
# visualize virus taxonomy results
rule virus_taxonomy_analysis:
    message:
        "Visualizing virus taxonomy results from geNomad"
    input:
        genomad=results + "09_VIRUS_TAXONOMY/01_genomad/votu_representatives_annotate/votu_representatives_taxonomy.tsv",
    output:
        svg=report(
            results + "09_VIRUS_TAXONOMY/virus_taxonomy_analysis.svg",
            category="Step 09: Virus taxonomy",
        ),
        report=results + "09_VIRUS_TAXONOMY/virus_taxonomy_report.tsv",
    params:
        min_genomad_agreement=config["virus_taxonomy"]["genomad_min_agreement"],
        min_genomad_genes=config["virus_taxonomy"]["genomad_min_genes"],
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/virus_taxonomy_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/09_virus_taxonomy_analysis.py"