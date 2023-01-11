# -----------------------------------------------------
# Virus Identification Module (If input_data = "reads" or "contigs")
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["samples_df"], sep="\t")
groups_samples = samples_df.loc[:, "group"] + "_" + samples_df.loc[:, "sample"]
samples_df["group_sample"] = samples_df["group"] + "_" + samples_df["sample"]

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
    rename_contigs_within_samples,
    virus_identification_analysis,


# -----------------------------------------------------
# 00 Determine input for module
# -----------------------------------------------------
# symlink contigs if contigs are input
rule symlink_contigs:
    input:
        lambda wildcards: samples_df[
            (samples_df["group"] + "_" + samples_df["sample"])
            == wildcards.group_sample
        ]["contigs"].iloc[0],
    output:
        temp(results + "00_INPUT/{group_sample}_contigs_symlink.fasta"),
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/symlink_contigs_{group_sample}.tsv"
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
        results + "00_INPUT/{group_sample}_contigs_symlink.fasta",
    output:
        temp(results + "00_INPUT/{group_sample}_contigs.fasta"),
    params:
        min_length=config["read_assembly"]["min_contig_length"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/filter_symlinked_contigs_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/03_contig_length_filter.py"


# if reads are input type, use output from 03_read_assembly
if config["input_data"] == "reads":
    assembly = (
        results
        + "03_READ_ASSEMBLY/02_contig_filters/{group_sample}/{group_sample}_contigs.fasta",
    )
# if contigs are input, use symlinked contigs
elif config["input_data"] == "contigs":
    assembly = (results + "00_INPUT/{group_sample}_contigs.fasta",)


# symlink preprocessed reads for external hits/abundances
rule symlink_preprocessed_reads:
    message:
        "Symlinking {wildcards.group_sample_replicate} input files to new location"
    input:
        R1=lambda wildcards: samples_df[
            (
                +samples_df["group"]
                + "_"
                + samples_df["sample"]
                + "_"
                + samples_df["replicate"]
            )
            == wildcards.group_sample_replicate
        ]["R1"].iloc[0],
        R2=lambda wildcards: samples_df[
            (
                +samples_df["group"]
                + "_"
                + samples_df["sample"]
                + "_"
                + samples_df["replicate"]
            )
            == wildcards.group_sample_replicate
        ]["R2"].iloc[0],
    output:
        R1=temp(results + "00_INPUT/{group_sample_replicate}.preprocessed_R1.fastq.gz"),
        R2=temp(results + "00_INPUT/{group_sample_replicate}.preprocessed_R2.fastq.gz"),
    benchmark:
        "benchmark/01_READ_PREPROCESSING/symlink_preprocessed_reads_{group_sample_replicate}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # symlink input reads to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# identify replicate reads
group_sample_replicate = samples_df.loc[:, ["sample", "replicate", "group"]]
group_sample_replicate["group_sample"] = (
    group_sample_replicate.loc[:, "group"]
    + "_"
    + group_sample_replicate.loc[:, "sample"]
)
group_sample_replicate_dictionary = group_sample_replicate.set_index(
    "group_sample"
).to_dict()["replicate"]


# merge sample replicates into single file
rule merge_preprocesssed_replicates:
    message:
        "Merging {wildcards.group_sample} replicates into one file"
    input:
        R1=lambda wildcards: expand(
            results + "00_INPUT/{{group_sample}}_{replicate}.preprocessed_R1.fastq.gz",
            replicate=group_sample_replicate_dictionary[wildcards.group_sample],
        ),
        R2=lambda wildcards: expand(
            results + "00_INPUT/{{group_sample}}_{replicate}.preprocessed_R2.fastq.gz",
            replicate=group_sample_replicate_dictionary[wildcards.group_sample],
        ),
    output:
        R1=temp(
            results
            + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R1.fastq.gz"
        ),
        R2=temp(
            results
            + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R2.fastq.gz"
        ),
    benchmark:
        "benchmark/01_READ_PREPROCESSING/merge_preprocessed_replicates_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # symlink replicates into one combined file
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# if reads are input, use preprocessed reads from 01_read_preprocessing for external hits/abundances
if config["input_data"] == "reads":
    R1 = results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz"
    R2 = results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz"
# if contigs, vls, or viruses are input then symlink reads for external hits/abundances
else:
    R1 = results + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R1.fastq.gz"
    R2 = results + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R2.fastq.gz"


# -----------------------------------------------------
# 01 Identify external hits
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


external_input = []
if config["virus_identification"]["external_input"] == "reads":
    external_input.append(R1)
    external_input.append(R2)
elif config["virus_identification"]["external_input"] == "contigs":
    external_input.append(assembly)


# screen reads to identify external viruses
rule screen_reads_against_virusdb:
    message:
        "Screening reads against {input.sketch} to identify external viruses present in {wildcards.group_sample}"
    input:
        external_input=external_input,
        sketch=config["virus_db"] + ".msh",
    output:
        temp(
            results
            + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/virusdb_mash_screen.tab"
        ),
    params:
        combined=results
        + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/combined_reads.fastq",
    # conda:
    #     "../envs/mash:2.3--ha61e061_0.yml"
    container:
        "docker://quay.io/biocontainers/mash:2.3--ha61e061_0"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/screen_reads_against_virusdb_sketch_{group_sample}.tsv"
    resources:
        runtime=config["virus_identification"]["mash_runtime"],
        mem_mb=config["virus_identification"]["mash_memory"],
    threads: config["virus_identification"]["mash_threads"]
    shell:
        """
        # combine reads
        cat {input.external_input} > {params.combined}

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
        "Extracting external viruses present in {wildcards.group_sample}"
    input:
        read_screen=results
        + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/virusdb_mash_screen.tab",
        virusdb=config["virus_db"],
    output:
        temp(
            results
            + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/virusdb_hits.fna"
        ),
    params:
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_hashes=config["virus_identification"]["min_mash_hashes"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/extract_virusdb_hits_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/04_extract_virusdb_hits.py"


# combine external hits with assembly contigs
rule combine_external_with_assembly:
    message:
        "Combining external hits in {wildcards.group_sample} with assembled contigs"
    input:
        contigs=assembly,
        external_hits=results
        + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/virusdb_hits.fna",
    output:
        temp(
            results
            + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/virusdb_hits_w_assembly.fna"
        ),
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/combine_external_with_assembly_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        cat {input} > {output}
        """


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
        "Running geNomad on {wildcards.group_sample} to identify viral contigs"
    input:
        genomad_db=resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
        contigs=results
        + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/virusdb_hits_w_assembly.fna",
    output:
        report=temp(
            results
            + "04_VIRUS_IDENTIFICATION/02_genomad/{group_sample}/virusdb_hits_w_assembly_summary/virusdb_hits_w_assembly_virus_summary.tsv"
        ),
        sequences=temp(
            results
            + "04_VIRUS_IDENTIFICATION/02_genomad/{group_sample}/virusdb_hits_w_assembly_summary/virusdb_hits_w_assembly_virus.fna"
        ),
    params:
        out_dir=results + "04_VIRUS_IDENTIFICATION/02_genomad/{group_sample}/",
        genomad_dir=resources + "genomad/genomad_db",
        min_score=config["virus_identification"]["genomad_min_score"],
        max_fdr=config["virus_identification"]["genomad_max_fdr"],
    conda:
        "../envs/genomad:1.3.0--pyhdfd78af_0.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/genomad_{group_sample}.tsv"
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
# 07 Combine outputs
# -----------------------------------------------------
# combine outputs from all tools
rule merge_reports_within_samples:
    message:
        "Merging all virus identification reports within {wildcards.group_sample}"
    input:
        genomad_results=results
        + "04_VIRUS_IDENTIFICATION/02_genomad/{group_sample}/virusdb_hits_w_assembly_summary/virusdb_hits_w_assembly_virus_summary.tsv",
        external_results=results
        + "04_VIRUS_IDENTIFICATION/01_external_hits/{group_sample}/virusdb_mash_screen.tab",
    output:
        temp(
            results
            + "04_VIRUS_IDENTIFICATION/03_combine_outputs/{group_sample}/combined_report.tsv"
        ),
    params:
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_hashes=config["virus_identification"]["min_mash_hashes"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
        assembly="{group_sample}",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/merge_reports_within_samples_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/04_merge_reports_within_samples.py"


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
            + "04_VIRUS_IDENTIFICATION/03_combine_outputs/{group_sample}/combined_report.tsv",
            group_sample=groups_samples,
        ),
    output:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.tsv",
    params:
        external_dirs=results + "04_VIRUS_IDENTIFICATION/01_external_hits/",
        genomad_dirs=results + "04_VIRUS_IDENTIFICATION/02_genomad/",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/combine_reports_across_samples.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}

        # remove 
        rm -rf {params.external_dirs}
        rm -rf {params.genomad_dirs}
        """


# combine virus reports across samples
rule rename_contigs_within_samples:
    message:
        "Renaming viral contigs for {wildcards.group_sample}"
    input:
        results
        + "04_VIRUS_IDENTIFICATION/02_genomad/{group_sample}/virusdb_hits_w_assembly_summary/virusdb_hits_w_assembly_virus.fna",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/03_combine_outputs/{group_sample}/combined_viral_contigs.fna",
    params:
        assembly="{group_sample}",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/rename_contigs_within_samples_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/04_merge_viral_contigs_within_samples.py"


# plot virus counts by tool
rule virus_identification_analysis:
    message:
        "Visualizing virus identification"
    input:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.tsv",
    output:
        report(
            results + "04_VIRUS_IDENTIFICATION/virus_identification_boxplot.svg",
            category="Step 04: Virus identification",
        ),
    params:
        genomad_score=config["virus_identification"]["genomad_min_score"],
        genomad_fdr=config["virus_identification"]["genomad_max_fdr"],
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
