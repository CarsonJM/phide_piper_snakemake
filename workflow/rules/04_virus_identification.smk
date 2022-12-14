# -----------------------------------------------------
# Virus Identification Module (If input_data = "reads" or "contigs")
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
# Virus identification rules
# -----------------------------------------------------
localrules:
    symlink_contigs,
    symlink_preprocessed_reads,
    combine_reports_across_samples,
    merge_reports_within_samples,
    merge_viral_contigs_within_samples,
    virus_identification_analysis,


# -----------------------------------------------------
# 00 Determine input for module
# -----------------------------------------------------
# symlink contigs if contigs are input
rule symlink_contigs:
    input:
        lambda wildcards: samples_df[(samples_df["sample"]) == wildcards.sample][
            "contigs"
        ].iloc[0],
    output:
        results + "00_INPUT/{sample}_contigs.fasta",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/symlink_contigs_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # symlink input contigs to renamed files
        ln -s {input} {output}
        """


# filter symlinked contigs based on contig length
rule filter_symlinked_contigs:
    input:
        results + "00_INPUT/{sample}_contigs.fasta",
    output:
        results + "00_INPUT/{sample}_contigs.fasta",
    params:
        min_length=config["read_assembly"]["min_contig_length"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/filter_symlinked_contigs_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/03_contig_length_filter.py"


# if reads are input type, use output from 03_read_assembly
if config["input_data"] == "reads":
    assembly = (
        results + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs.fasta",
    )
# if contigs are input, use symlinked contigs
elif config["input_data"] == "contigs":
    assembly = (results + "00_INPUT/{sample}_contigs.fasta",)


# symlink preprocessed reads for external hits/abundances
rule symlink_preprocessed_reads:
    input:
        R1=lambda wildcards: samples_df[(+samples_df["sample"]) == wildcards.sample][
            "R1"
        ].iloc[0],
        R2=lambda wildcards: samples_df[(+samples_df["sample"]) == wildcards.sample][
            "R2"
        ].iloc[0],
    output:
        R1=results + "00_INPUT/{sample}_proprcessed_1.fastq.gz",
        R2=results + "00_INPUT/{sample}_preprocessed_2.fastq.gz",
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # symlink input reads to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# if reads are input, use preprocessed reads from 01_read_preprocessing for external hits/abundances
if config["input_data"] == "reads":
    R1 = results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz"
    R2 = results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz"
# if contigs, vls, or viruses are input then symlink reads for external hits/abundances
else:
    R1 = results + "00_INPUT/{sample}_proprocessed_1.fastq.gz"
    R2 = results + "00_INPUT/{sample}_preprocessed_2.fastq.gz"


# -----------------------------------------------------
# 01 geNomad
# -----------------------------------------------------
# download genomad database
rule download_genomad:
    message:
        "Downloading geNomad database"
    output:
        resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
    params:
        genomad_dir=resources + "genomad/",
    conda:
        "../envs/genomad:1.3.0--pyhdfd78af_0.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/download_genomad.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="10000",
    shell:
        """
        # change to genomad directory
        cd {params.genomad_dir}

        # download genomad databases
        genomad download-database .
        """


# run genomad to identify viral contigs
rule genomad:
    message:
        "Running geNomad on {wildcards.sample} to identify viral contigs"
    input:
        genomad_db=resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
        contigs=assembly,
    output:
        summary=results
        + "04_VIRUS_IDENTIFICATION/01_genomad/{sample}/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv",
        viruses=results
        + "04_VIRUS_IDENTIFICATION/01_genomad/{sample}/{sample}_contigs_summary/{sample}_contigs_virus.fna",
    params:
        out_dir=results + "04_VIRUS_IDENTIFICATION/01_genomad/{sample}/",
        genomad_dir=resources + "genomad/genomad_db",
        min_score=config["virus_identification"]["genomad_min_score"],
        max_fdr=config["virus_identification"]["genomad_max_fdr"],
    conda:
        "../envs/genomad:1.3.0--pyhdfd78af_0.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/genomad_{sample}.tsv"
    resources:
        runtime=config["virus_identification"]["genomad_runtime"],
        mem_mb=config["virus_identification"]["genomad_memory"],
    threads: config["virus_identification"]["genomad_threads"]
    shell:
        """
        # run genomad
        genomad end-to-end \
        --threads {threads} \
        --verbose \
        --min-score {params.min_score} \
        --cleanup \
        --splits {threads} \
        {input.contigs} \
        {params.out_dir} \
        {params.genomad_dir}
        """


# -----------------------------------------------------
# 02 Identify external hits
# -----------------------------------------------------
# create a mash sketch of virusdb
rule mash_sketch_virusdb:
    message:
        "Creating mash sketch of {input}"
    input:
        config["virus_db"],
    output:
        config["virus_db"] + ".msh",
    # conda:
    #     "../envs/mash:2.3--ha61e061_0.yml"
    container:
        "docker://quay.io/biocontainers/mash:2.3--ha61e061_0"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mash_sketch_virusdb.tsv"
    resources:
        runtime=config["virus_identification"]["mash_runtime"],
        mem_mb=config["virus_identification"]["mash_memory"],
    threads: config["virus_identification"]["mash_threads"]
    shell:
        """
        # create a mash sketch of virusdb
        mash sketch \
        -p {threads} \
        -i {input}
        """


# screen reads to identify external viruses
rule screen_reads_against_virusdb:
    message:
        "Screening reads against {input.sketch} to identify external viruses present in {wildcards.sample}"
    input:
        R1=R1,
        R2=R2,
        sketch=config["virus_db"] + ".msh",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/02_external_hits/{sample}/virusdb_mash_screen.tab",
    params:
        combined=results
        + "04_VIRUS_IDENTIFICATION/02_external_hits/{sample}/combined_reads.fastq",
    # conda:
    #     "../envs/mash:2.3--ha61e061_0.yml"
    container:
        "docker://quay.io/biocontainers/mash:2.3--ha61e061_0"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/screen_reads_against_virusdb_sketch_{sample}.tsv"
    resources:
        runtime=config["virus_identification"]["mash_runtime"],
        mem_mb=config["virus_identification"]["mash_memory"],
    threads: config["virus_identification"]["mash_threads"]
    shell:
        """
        # combine reads
        cat {input.R1} {input.R2} > {params.combined}

        # screen reads against virusdb
        mash screen \
        -p {threads} \
        {input.sketch} \
        {params.combined} > {output}

        # remove combined fastq to save space
        rm {params.combined}
        """


# filter to keep only external hits
rule extract_external_hits:
    message:
        "Extracting external viruses present in {wildcards.sample}"
    input:
        read_screen=results
        + "04_VIRUS_IDENTIFICATION/02_external_hits/{sample}/virusdb_mash_screen.tab",
        virusdb=config["virus_db"],
    output:
        results + "04_VIRUS_IDENTIFICATION/02_external_hits/{sample}/virusdb_hits.fna",
    params:
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_hashes=config["virus_identification"]["min_mash_hashes"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/extract_virusdb_hits_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/04_extract_virusdb_hits.py"


# -----------------------------------------------------
# 07 Combine outputs
# -----------------------------------------------------
# determine input files for detecting virus sequences
if config["virus_identification"]["run_genomad"]:
    genomad = (
        results
        + "04_VIRUS_IDENTIFICATION/01_genomad/{sample}/{sample}_contigs_summary/{sample}_contigs_virus_summary.tsv"
    )
else:
    genomad = pd.DataFrame()

if config["virus_identification"]["run_external"] and config["input_data"] == "reads":
    external = (
        results
        + "04_VIRUS_IDENTIFICATION/02_external_hits/{sample}/virusdb_mash_screen.tab"
    )
else:
    external = pd.DataFrame()


# combine outputs from all tools
rule merge_reports_within_samples:
    message:
        "Merging all virus identification reports within {wildcards.sample}"
    input:
        genomad_results=genomad,
        external_results=external,
    output:
        results
        + "04_VIRUS_IDENTIFICATION/03_combine_outputs/{sample}/combined_report.csv",
    params:
        run_genomad=config["virus_identification"]["run_genomad"],
        run_external=config["virus_identification"]["run_external"],
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_hashes=config["virus_identification"]["min_mash_hashes"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
        assembly="{sample}",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/merge_reports_within_samples_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/04_merge_reports_within_samples.py"


# determine which inputs are required for virus identification
vls_contigs = []
if config["virus_identification"]["run_external"]:
    vls_contigs.append(
        results + "04_VIRUS_IDENTIFICATION/02_external_hits/{sample}/virusdb_hits.fna",
    )

if config["virus_identification"]["run_genomad"]:
    vls_contigs.append(
        results
        + "04_VIRUS_IDENTIFICATION/01_genomad/{sample}/{sample}_contigs_summary/{sample}_contigs_virus.fna",
    )


# combine viral contigs from all tool outputs using thresholds specified in config.yaml
rule merge_viral_contigs_within_samples:
    message:
        "Merging viral contigs meeting config.yaml criteria within {wildcards.sample}"
    input:
        vls_contigs,
    output:
        results
        + "04_VIRUS_IDENTIFICATION/03_combine_outputs/{sample}/combined_viral_contigs.fasta",
    params:
        run_genomad=config["virus_identification"]["run_genomad"],
        run_external=config["virus_identification"]["run_external"],
        assembly="{sample}",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/merge_viral_contigs_within_samples_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        cat {input} > {output}
        """


# -----------------------------------------------------
# Analyze combined virus data
# -----------------------------------------------------
# combine virus reports across samples
rule combine_reports_across_samples:
    message:
        "Combining viral reports across all samples"
    input:
        expand(
            results
            + "04_VIRUS_IDENTIFICATION/03_combine_outputs/{sample}/combined_report.csv",
            sample=samples,
        ),
    output:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.csv",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/combine_reports_across_samples.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# plot virus counts by tool
rule virus_identification_analysis:
    message:
        "Visualizing virus identification"
    input:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.csv",
    output:
        report(
            results + "04_VIRUS_IDENTIFICATION/virus_identification_boxplot.svg",
            category="Step 04: Virus identification",
        ),
    params:
        run_genomad=config["virus_identification"]["run_genomad"],
        genomad_score=config["virus_identification"]["genomad_min_score"],
        genomad_fdr=config["virus_identification"]["genomad_max_fdr"],
        run_external=config["virus_identification"]["run_external"],
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_hashes=config["virus_identification"]["min_mash_hashes"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/virus_identification_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/04_virus_identification_analysis.py"
