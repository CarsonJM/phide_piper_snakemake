# -------------------------------------
# Read Preprocessing Module (if input_data = "reads")
# -------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["samples_df"], sep="\t")
groups_samples = samples_df.loc[:, "group"] + "_" + samples_df.loc[:, "sample"]


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "../report/workflow.rst"


# -------------------------------------
# Preprocessing Rules
# -------------------------------------
localrules:
    symlink_input_reads,
    merge_input_replicates,
    kneaddata_analysis,


# -----------------------------------------------------
# 00 Symlink inputs
# -----------------------------------------------------
# symlink input reads to new paths
rule symlink_input_reads:
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
        R1=temp(results + "00_INPUT/{group_sample_replicate}.R1.fastq.gz"),
        R2=temp(results + "00_INPUT/{group_sample_replicate}.R2.fastq.gz"),
    benchmark:
        "benchmark/01_READ_PREPROCESSING/symlink_reads_{group_sample_replicate}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # symlink input reads to specified paths
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 01 Merge Replicates
# -----------------------------------------------------
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
rule merge_input_replicates:
    message:
        "Merging {wildcards.group_sample} replicates into one file"
    input:
        R1=lambda wildcards: expand(
            results + "00_INPUT/{{group_sample}}_{replicate}.R1.fastq.gz",
            replicate=group_sample_replicate_dictionary[wildcards.group_sample],
        ),
        R2=lambda wildcards: expand(
            results + "00_INPUT/{{group_sample}}_{replicate}.R2.fastq.gz",
            replicate=group_sample_replicate_dictionary[wildcards.group_sample],
        ),
    output:
        R1=temp(
            results
            + "01_READ_PREPROCESSING/01_merge_replicates/{group_sample}.R1.fastq.gz"
        ),
        R2=temp(
            results
            + "01_READ_PREPROCESSING/01_merge_replicates/{group_sample}.R2.fastq.gz"
        ),
    benchmark:
        "benchmark/01_READ_PREPROCESSING/merge_replicates_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # symlink replicates into one combined file
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 02 fastp
# -----------------------------------------------------
# trim and deduplicate reads with fastp
rule fastp:
    message:
        "Trimming and deduplicating {wildcards.group_sample} with fastp"
    input:
        R1=results
        + "01_READ_PREPROCESSING/01_merge_replicates/{group_sample}.R1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/01_merge_replicates/{group_sample}.R2.fastq.gz",
    output:
        R1=temp(results + "01_READ_PREPROCESSING/02_fastp/{group_sample}.R1.fastq.gz"),
        R2=temp(results + "01_READ_PREPROCESSING/02_fastp/{group_sample}.R2.fastq.gz"),
        report=temp(
            results + "01_READ_PREPROCESSING/02_fastp/{group_sample}_fastp.json"
        ),
        html=temp(results + "01_READ_PREPROCESSING/02_fastp/{group_sample}_fastp.html"),
    params:
        extra_args=config["read_preprocessing"]["fastp_arguments"],
        html=results + "01_READ_PREPROCESSING/02_fastp/{group_sample}_fastp.html",
    # conda:
    #     "../envs/fastp:0.23.2--h79da9fb_0.yml"
    container:
        "docker://quay.io/biocontainers/fastp:0.23.2--h79da9fb_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/fastp_{group_sample}.tsv"
    resources:
        runtime=config["read_preprocessing"]["fastp_runtime"],
        mem_mb=config["read_preprocessing"]["fastp_memory"],
    threads: config["read_preprocessing"]["fastp_threads"]
    shell:
        """
        # run fastp to trim and remove 
        fastp --in1 {input.R1} \
        --in2 {input.R2} \
        --out1 {output.R1} \
        --out2 {output.R2} \
        --json {output.report} \
        --html {params.html} \
        --thread {threads} \
        {params.extra_args}
        """


# generate report for trimming and deduplication
rule fastp_multiqc:
    message:
        "Generating a MULTIQC report using fastp results"
    input:
        json=expand(
            results + "01_READ_PREPROCESSING/02_fastp/{group_sample}_fastp.json",
            group_sample=groups_samples,
        ),
        html=expand(
            results + "01_READ_PREPROCESSING/02_fastp/{group_sample}_fastp.html",
            group_sample=groups_samples,
        ),
    output:
        report(
            results + "01_READ_PREPROCESSING/fastp_multiqc_report.html",
            category="Step 01: Read preprocessing",
        ),
    params:
        fastp_dir=results + "01_READ_PREPROCESSING/02_fastp/",
        fastp_out=results + "01_READ_PREPROCESSING/02_fastp/multiqc_report.html",
    # conda:
    #     "../envs/multiqc:1.12--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/fastp_multiqc.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        multiqc {params.fastp_dir} \
        -o {params.fastp_dir} -f

        mv {params.fastp_out} {output}
        """


# -----------------------------------------------------
# 03 KneadData
# -----------------------------------------------------
# download and build kneaddata human bowtie2 database
rule download_kneaddata_database:
    message:
        "Downloading human genome bowtie2 database for KneadData"
    output:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
    params:
        kneaddata_db=resources + "kneaddata/",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/download_kneaddata_database.tsv"
    resources:
        runtime="00:30:00",
        mem_mb="10000",
    shell:
        """
        # download human genome reference to desired directory
        kneaddata_database --download human_genome bowtie2 {params.kneaddata_db}
        """


# Remove human reads with kneaddata
rule kneaddata:
    message:
        "Running KneadData on {wildcards.group_sample} to remove human reads"
    input:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
        R1=results + "01_READ_PREPROCESSING/02_fastp/{group_sample}.R1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/02_fastp/{group_sample}.R2.fastq.gz",
    output:
        log=results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}.log",
        R1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz",
    params:
        out_dir=results + "01_READ_PREPROCESSING/03_kneaddata/",
        human_db=resources + "kneaddata/",
        extra_args=config["read_preprocessing"]["kneaddata_arguments"],
        prefix="{group_sample}",
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq",
        hg37_paired1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_hg37dec_v0.1_bowtie2_paired_contam_1.fastq",
        hg37_paired2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_hg37dec_v0.1_bowtie2_paired_contam_2.fastq",
        hg37_unmatched1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_hg37dec_v0.1_bowtie2_unmatched_contam_1.fastq",
        hg37_unmatched2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_hg37dec_v0.1_bowtie2_unmatched_contam_1.fastq",
        unmatched1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_unmatched_1.fastq",
        unmatched2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_unmatched_2.fastq",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/kneaddata_{group_sample}.tsv"
    resources:
        runtime=config["read_preprocessing"]["kneaddata_runtime"],
        mem_mb=config["read_preprocessing"]["kneaddata_memory"],
    threads: config["read_preprocessing"]["kneaddata_threads"]
    shell:
        """
        # run kneaddata to quality filter and remove host reads
        kneaddata --input {input.R1} --input {input.R2} \
        --output {params.out_dir} \
        --output-prefix {params.prefix} \
        --reference-db {params.human_db} \
        --threads {threads} \
        --bypass-trim \
        --bypass-trf \
        --log {output.log} \
        {params.extra_args}


        gzip {params.R1}
        gzip {params.R2}
        gzip {params.unmatched1}
        gzip {params.unmatched2}

        rm {params.hg37_paired1} {params.hg37_paired2} {params.hg37_unmatched1} {params.hg37_unmatched2}
        """


# Count reads before and after host removal
rule kneaddata_read_counts:
    message:
        "Generating read counts before and after human removal"
    input:
        expand(
            results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}.log",
            group_sample=groups_samples,
        ),
    output:
        results + "01_READ_PREPROCESSING/kneaddata_read_counts.tsv",
    params:
        in_dir=results + "01_READ_PREPROCESSING/03_kneaddata/",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/kneaddata_read_counts.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # generate read counts from kneaddata log files
        kneaddata_read_count_table \
        --input {params.in_dir} \
        --output {output}
        """


# -----------------------------------------------------
# 04 Analysis
# -----------------------------------------------------
# generate report for human read removal
rule kneaddata_analysis:
    message:
        "Visualizing read counts before and after human read removal"
    input:
        results + "01_READ_PREPROCESSING/kneaddata_read_counts.tsv",
    output:
        svg=report(
            results + "01_READ_PREPROCESSING/kneaddata_analysis.svg",
            category="Step 01: Read preprocessing",
        ),
        html=results + "01_READ_PREPROCESSING/kneaddata_analysis.html",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/kneaddata_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/01_kneaddata_analysis.py"
