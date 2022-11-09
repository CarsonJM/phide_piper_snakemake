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
localrules: combine_genomad_lifestyles,
            add_sample_to_genomad_lifestyle,
            combine_genomad_lifestyles,
            combine_virus_lifestyle_outputs,
            virus_lifestyle_analysis

# -----------------------------------------------------
# 01 BACPHLIP
# -----------------------------------------------------
# build bacphlip
rule build_bacphlip:
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
rule bacphlip:
    input:
        bacphlip=resources + "bacphlip/bacphlip_build_complete",
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        bac_final=results + "10_VIRUS_LIFESTYLE/01_bacphlip/lifestyles.tsv",
        hmm_final=results + "10_VIRUS_LIFESTYLE/01_bacphlip/hmmsearch.tsv",
    params:
        bac_out=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna.bacphlip",
        hmm_out=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna.hmmsearch.tsv",
        ind_out=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna.BACPHLIP_DIR",
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
        partition="compute-hugemem"
    threads: config["virus_lifestyle"]["bacphlip_threads"]
    shell:
        """
        # run bacphlip
        bacphlip -i {input.viruses} \
        --multi_fasta

        # move output files to bacphlip dir
        mv {params.bac_out} {output.bac_final}
        mv {params.hmm_out} {output.hmm_final}
        mv {params.ind_out} {params.ind_final}
        """


# -----------------------------------------------------
# Analyze virus lifestyles
# -----------------------------------------------------
# visualize virus host results
rule virus_lifestyle_analysis:
    input:
        results + "10_VIRUS_LIFESTYLE/01_bacphlip/lifestyles.tsv",
    output:
        svg=report(
            results + "10_VIRUS_LIFESTYLE/virus_lifestyle_analysis.svg",
            category="Step 08: Virus lifestyle",
        ),
        html=results + "10_VIRUS_LIFESTYLE/virus_lifestyle_analysis.html",
    params:
        bacphlip_prob=config["virus_lifestyle"]["bacphlip_confidence"],
    conda:
        "../envs/jupyter.yml"
    resources:
        runtime="00:30:00",
        mem_mb="1000",
    benchmark:
        "benchmark/10_VIRUS_LIFESTYLE/virus_lifestyle_analysis.tsv"
    script:
        "../scripts/10_virus_lifestyle_analysis.py"
