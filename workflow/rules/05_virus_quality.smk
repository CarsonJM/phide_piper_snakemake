# -----------------------------------------------------
# Virus Quality Module (if input_data = "reads" or "contigs" or "vls")
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["samples_df"], sep="\t")
samples = samples_df["sample"]


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "../report/workflow.rst"


# -----------------------------------------------------
# Virus quality rules
# -----------------------------------------------------
localrules:
    symlink_viruses,
    combine_virus_quality_reports,


# -----------------------------------------------------
# 00 Determine inputs for module
# -----------------------------------------------------
# symlink input vls
rule symlink_vls:
    message:
        "Symlinking input virus-like sequences from {wildcards.sample} to new location"
    input:
        lambda wildcards: samples_df[(samples_df["sample"]) == wildcards.sample][
            "vls"
        ].iloc[0],
    output:
        results + "00_INPUT/{sample}_vls.fasta",
    benchmark:
        "benchmark/05_VIRUS_QUALITY/symlink_vls_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # symlink input vls to renamed files
        ln -s {input} {output}
        """


# if reads or contigs are input, use output from 04_virus_identification
if config["input_data"] == "reads" or config["input_data"] == "contigs":
    vls = (
        results
        + "04_VIRUS_IDENTIFICATION/08_combine_outputs/{sample}/combined_viral_contigs.fasta"
    )
# if input is vls, use symlinked vls
elif config["input_data"] == "vls":
    vls = results + "00_INPUT/{sample}_vls.fasta"


# -----------------------------------------------------
# 01 CheckV: Trim Viruses
# -----------------------------------------------------
# download checkv database
rule download_checkv:
    message:
        "Downloading CheckV database"
    output:
        resources + "checkv/checkv-db-v1.4/genome_db/checkv_reps.dmnd",
    params:
        checkv_dir=resources + "checkv/",
    # conda:
    #     "../envs/checkv:1.0.1--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/checkv:1.0.1--pyhdfd78af_0"
    benchmark:
        "benchmark/05_VIRUS_QUALITY/download_checkv.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="10000",
    shell:
        """
        # download checkv database
        checkv download_database {params.checkv_dir}
        """


# run checkv on viral contigs to trim proviruses
rule checkv_trim:
    message:
        "Trimming host regions from {wildcards.sample} viruses"
    input:
        checkv_db=resources + "checkv/checkv-db-v1.4/genome_db/checkv_reps.dmnd",
        virus_contigs=vls,
    output:
        checkv_proviruses=results
        + "05_VIRUS_QUALITY/01_checkv_trim/{sample}/proviruses.fna",
        checkv_viruses=results + "05_VIRUS_QUALITY/01_checkv_trim/{sample}/viruses.fna",
    params:
        checkv_dir=results + "05_VIRUS_QUALITY/01_checkv_trim/{sample}/",
        checkv_db=resources + "checkv/checkv-db-v1.4",
    # conda:
    #     "../envs/checkv:1.0.1--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/checkv:1.0.1--pyhdfd78af_0"
    benchmark:
        "benchmark/05_VIRUS_QUALITY/checkv_trim_{sample}.tsv"
    resources:
        runtime=config["virus_quality"]["checkv_runtime"],
        mem_mb=config["virus_quality"]["checkv_memory"],
    threads: config["virus_quality"]["checkv_threads"]
    shell:
        """
        # run checkv to trim proviruses
        checkv contamination {input.virus_contigs} {params.checkv_dir} \
        -d {params.checkv_db} \
        -t {threads}
        """


# combine trimmed proviruses and input vls
rule combine_trimmed_vls:
    message:
        "Renaming trimmed viruses and combining with untrimmed viruses from {wildcards.sample}"
    input:
        checkv_proviruses=results
        + "05_VIRUS_QUALITY/01_checkv_trim/{sample}/proviruses.fna",
        checkv_viruses=results + "05_VIRUS_QUALITY/01_checkv_trim/{sample}/viruses.fna",
    output:
        results + "05_VIRUS_QUALITY/01_checkv_trim/{sample}/trimmed_vls.fna",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/05_VIRUS_QUALITY/combine_trimmed_viruses_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/05_combine_trimmed_viruses.py"


# if remove_proviruses = True, do not use trimmed viruses as proviruses will be removed
if config["virus_quality"]["remove_proviruses"]:
    checkv_input = vls
# if remove_proviruses = False, use trimmed viruses to quality filter proviruses
else:
    checkv_input = results + "05_VIRUS_QUALITY/01_checkv_trim/{sample}/trimmed_vls.fna"


# -----------------------------------------------------
# 02 CheckV: Virus quality
# -----------------------------------------------------
# determine virus quality using CheckV
rule checkv:
    message:
        "Running CheckV to determine virus quality and contamination for {wildcards.sample}"
    input:
        checkv_db=resources + "checkv/checkv-db-v1.4/genome_db/checkv_reps.dmnd",
        virus_contigs=checkv_input,
    output:
        checkv_results=results
        + "05_VIRUS_QUALITY/02_checkv/{sample}/quality_summary.tsv",
        checkv_viruses=results + "05_VIRUS_QUALITY/02_checkv/{sample}/viruses.fna",
    params:
        checkv_dir=results + "05_VIRUS_QUALITY/02_checkv/{sample}/",
        checkv_db=resources + "checkv/checkv-db-v1.4",
    # conda:
    #     "../envs/checkv:1.0.1--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/checkv:1.0.1--pyhdfd78af_0"
    benchmark:
        "benchmark/05_VIRUS_QUALITY/checkv_{sample}.tsv"
    resources:
        runtime=config["virus_quality"]["checkv_runtime"],
        mem_mb=config["virus_quality"]["checkv_memory"],
    threads: config["virus_quality"]["checkv_threads"]
    shell:
        """
        # run checkv to determine virus quality
        checkv end_to_end {input.virus_contigs} {params.checkv_dir} \
        -d {params.checkv_db} \
        -t {threads}

        # add sample column to each checkv output
        s={wildcards.sample}
        sed -i "s/$/\t$s/" {output.checkv_results}
        sample="sample"
        sed -i "1s/$s/$sample/" {output.checkv_results}
        """


# -----------------------------------------------------
# 03 Quality filter viruses
# -----------------------------------------------------
# remove low quality virus (and retain high-quality untrimmed counterparts)
rule quality_filter_viruses:
    message:
        "Filtering viruses in {wildcards.sample} not meeting config.yaml criteria"
    input:
        checkv_results=results
        + "05_VIRUS_QUALITY/02_checkv/{sample}/quality_summary.tsv",
        checkv_viruses=results + "05_VIRUS_QUALITY/02_checkv/{sample}/viruses.fna",
        untrimmed_viruses=vls,
    output:
        viruses=results
        + "05_VIRUS_QUALITY/03_quality_filter/{sample}/quality_filtered_viruses.fna",
        untrimmed_viruses=results
        + "05_VIRUS_QUALITY/03_quality_filter/{sample}/untrimmed_quality_filtered_viruses.fna",
    params:
        min_completeness=config["virus_quality"]["min_completeness"],
        min_viral_genes=config["virus_quality"]["min_viral_genes"],
        max_bacterial_genes=config["virus_quality"]["max_bacterial_genes"],
        remove_proviruses=config["virus_quality"]["remove_proviruses"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/05_VIRUS_QUALITY/quality_filter_viruses_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/05_quality_filter_viruses.py"


# -----------------------------------------------------
# Analyze virus quality
# -----------------------------------------------------
# combine checkv reports across samples
rule combine_checkv_reports:
    message:
        "Combining CheckV reports across samples"
    input:
        expand(
            results + "05_VIRUS_QUALITY/02_checkv/{sample}/quality_summary.tsv",
            sample=samples,
        ),
    output:
        results + "05_VIRUS_QUALITY/virus_quality_report.tsv",
    benchmark:
        "benchmark/05_VIRUS_QUALITY/combine_checkv_reports.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# analyze checkv results to visualize genome qualities and provirus prevalence
rule virus_quality_analysis:
    message:
        "Visualizing virus quality for complete dataset"
    input:
        results + "05_VIRUS_QUALITY/virus_quality_report.tsv",
    output:
        svg=report(
            results + "05_VIRUS_QUALITY/virus_quality_figure.svg",
            category="Step 05: Virus quality",
        ),
        html=results + "05_VIRUS_QUALITY/virus_quality_figure.html",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/05_VIRUS_QUALITY/virus_quality_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/05_virus_quality_analysis.py"
