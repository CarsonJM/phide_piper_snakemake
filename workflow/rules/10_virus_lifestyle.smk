# -----------------------------------------------------
# Virus lifestyle
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples_df"], sep="\t")
samples = samples_df['sample']


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -----------------------------------------------------
# Virus lifestyle rules
# -----------------------------------------------------
localrules: extract_hq_for_bacphlip

# -----------------------------------------------------
# 01 BACPHLIP
# -----------------------------------------------------
# build bacphlip
rule build_bacphlip:
    message:
        "Building BACPHLIP to determine virus lifestyles"
    output:
        resources + "bacphlip/bacphlip_build_complete",
    params:
        bac_dir=resources + "bacphlip",
    benchmark:
        "benchmark/10_VIRUS_LIFESTYLE/build_bacphlip.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    conda:
        "../envs/bacphlip:0.9.6--py_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/bacphlip:0.9.6--py_0"
    shell:
        """
        pip install bacphlip

        touch {output}
        """

# predict virus lifestyle using bacphlip
rule extract_hq_for_bacphlip:
    message:
        "Extracting HQ viruses to run BACPHLIP on"
    input:
        checkv=results + "05_VIRUS_QUALITY/virus_quality_report.tsv",
        clusters=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_clusters.tsv",
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        results + "10_VIRUS_LIFESTYLE/01_bacphlip/hq_viruses.fna",
    params:
        min_completeness=config["virus_lifestyle"]["min_completeness"]
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/10_VIRUS_LIFESTYLE/bacpextract_hq_for_bacphliphlip.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script: "../scripts/10_extract_hq_viruses.py"


# predict virus lifestyle using bacphlip
rule bacphlip:
    message:
        "Running BACPHLIP to determine virus lifestyles"
    input:
        bacphlip=resources + "bacphlip/bacphlip_build_complete",
        viruses=results + "10_VIRUS_LIFESTYLE/01_bacphlip/hq_viruses.fna",
    output:
        bac_final=results + "10_VIRUS_LIFESTYLE/01_bacphlip/hq_viruses.fna.bacphlip",
        hmm_final=results + "10_VIRUS_LIFESTYLE/01_bacphlip/hq_viruses.fna.hmmsearch.tsv",
    params:
        ind_final=results + "10_VIRUS_LIFESTYLE/01_bacphlip/",
    conda:
        "../envs/bacphlip:0.9.6--py_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/bacphlip:0.9.6--py_0"
    benchmark:
        "benchmark/10_VIRUS_LIFESTYLE/bacphlip.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="10000",
    threads: config["virus_lifestyle"]["bacphlip_threads"]
    shell:
        """
        # run bacphlip
        bacphlip -i {input.viruses} \
        --multi_fasta
        """


# -----------------------------------------------------
# Analyze virus lifestyles
# -----------------------------------------------------
# visualize virus host results
rule virus_lifestyle_analysis:
    message:
        "Visualizing virus lifestyle outputs as determined using BACPHLIP"
    input:
        results + "10_VIRUS_LIFESTYLE/01_bacphlip/hq_viruses.fna.bacphlip",
    output:
        svg=report(
            results + "10_VIRUS_LIFESTYLE/virus_lifestyle_analysis.svg",
            category="Step 10: Virus lifestyle",
        ),
        html=results + "10_VIRUS_LIFESTYLE/virus_lifestyle_analysis.html",
    params:
        bacphlip_prob=config["virus_lifestyle"]["bacphlip_confidence"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/10_VIRUS_LIFESTYLE/virus_lifestyle_analysis.tsv"
    resources:
        runtime="00:30:00",
        mem_mb="1000",
    script:
        "../scripts/10_virus_lifestyle_analysis.py"
