# -------------------------------------
# Read Preprocessing Module
# -------------------------------------
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


# -------------------------------------
# Preprocessing Rules
# -------------------------------------
localrules:
    symlink_reads,
    merge_replicates,


# -----------------------------------------------------
# 00 Symlink inputs
# -----------------------------------------------------
# symlink input paths to new paths
rule symlink_reads:
    input:
        R1=lambda wildcards: samples_df[
            (+samples_df["sample"] + "_" + samples_df["replicate"])
            == wildcards.sample_replicate
        ]["R1"].iloc[0],
        R2=lambda wildcards: samples_df[
            (+samples_df["sample"] + "_" + samples_df["replicate"])
            == wildcards.sample_replicate
        ]["R2"].iloc[0],
    output:
        R1=results + "00_INPUT/{sample_replicate}.R1.fastq.gz",
        R2=results + "00_INPUT/{sample_replicate}.R2.fastq.gz",
    benchmark:
        "benchmark/01_READ_PREPROCESSING/symlink_reads_{sample_replicate}.tsv"
    resources:
        runtime="00:00:10",
        mem_mb="1000",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 01 Merge Replicates
# -----------------------------------------------------
# identify replicates
sample_replicate = samples_df[["sample", "replicate"]]
sample_replicate_dictionary = sample_replicate.set_index("sample").to_dict()[
    "replicate"
]


# merge sample replicates into single file
rule merge_replicates:
    input:
        R1=lambda wildcards: expand(
            results + "00_INPUT/{{sample}}_{replicate}.R1.fastq.gz",
            replicate=sample_replicate_dictionary[wildcards.sample],
        ),
        R2=lambda wildcards: expand(
            results + "00_INPUT/{{sample}}_{replicate}.R2.fastq.gz",
            replicate=sample_replicate_dictionary[wildcards.sample],
        ),
    output:
        R1=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R2.fastq.gz",
    benchmark:
        "benchmark/01_READ_PREPROCESSING/merge_replicates_{sample}.tsv"
    resources:
        runtime="00:00:10",
        mem_mb="1000",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# -----------------------------------------------------
# 02 fastp
# -----------------------------------------------------
# trim and deduplicate reads with fastp
rule fastp:
    input:
        R1=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/01_merge_replicates/{sample}.R2.fastq.gz",
    output:
        R1=results + "01_READ_PREPROCESSING/02_fastp/{sample}.R1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/02_fastp/{sample}.R2.fastq.gz",
        report=results + "01_READ_PREPROCESSING/02_fastp/{sample}_fastp.json",
    params:
        extra_args=config["read_preprocessing"]["fastp_arguments"],
        html=results + "01_READ_PREPROCESSING/02_fastp/{sample}_fastp.html",
    # conda:
    #     "../envs/fastp:0.23.2--h79da9fb_0.yml"
    container:
        "docker://quay.io/biocontainers/fastp:0.23.2--h79da9fb_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/fastp_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["read_preprocessing"]["fastp_threads"]
    shell:
        """
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
    input:
        expand(
            results + "01_READ_PREPROCESSING/02_fastp/{sample}_fastp.json",
            sample=samples,
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
    input:
        resources + "kneaddata/hg37dec_v0.1.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.2.bt2",
        resources + "kneaddata/hg37dec_v0.1.3.bt2",
        resources + "kneaddata/hg37dec_v0.1.4.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.1.bt2",
        resources + "kneaddata/hg37dec_v0.1.rev.2.bt2",
        R1=results + "01_READ_PREPROCESSING/02_fastp/{sample}.R1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/02_fastp/{sample}.R2.fastq.gz",
    output:
        log=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}.log",
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz",
    params:
        out_dir=results + "01_READ_PREPROCESSING/03_kneaddata/",
        human_db=resources + "kneaddata/",
        extra_args=config["read_preprocessing"]["kneaddata_arguments"],
        prefix="{sample}",
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/01_READ_PREPROCESSING/kneaddata_{sample}.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="10000",
        partition="compute-hugemem",
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
        --run-fastqc-end \
        {params.extra_args}


        gzip {params.R1}
        gzip {params.R2}
        """


# Count reads after host removal
rule kneaddata_read_counts:
    input:
        expand(
            results + "01_READ_PREPROCESSING/03_kneaddata/{sample}.log", sample=samples
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
        runtime="00:01:00",
        mem_mb="1000",
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
        mem_mb="1000",
    script:
        "../scripts/01_kneaddata_analysis.py"
