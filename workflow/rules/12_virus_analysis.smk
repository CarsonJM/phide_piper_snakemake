# -----------------------------------------------------
# Virus abundance
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
# Virus abundance rules
# -----------------------------------------------------
localrules:
    symlink_preprocessed_reads,
    read_counts,
    combine_read_counts_across_samples,
    combine_breadth_and_depth_across_samples,


# -----------------------------------------------------
# 00 Determine inputs for module
# -----------------------------------------------------

# read symlink in 07_virus_clustering.smk
if config["input_data"] == "reads":
    R1 = results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz"
    R2 = results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz"
elif (
    config["input_data"] == "contigs"
    or config["input_data"] == "viruses"
    or config["input_data"] == "processed_viruses"
):
    R1 = results + "00_INPUT/{sample}_proprocessed_1.fastq.gz"
    R2 = results + "00_INPUT/{sample}_preprocessed_2.fastq.gz"


# -----------------------------------------------------
# 01 Create genome catalog to align to
# -----------------------------------------------------
# Align reads to virus catalog using bowtie2
rule build_viruses_bowtie2db:
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog.1.bt2",
    params:
        db=results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog",
    conda:
        "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    threads: 8
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/build_viruses_bowtie2db.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="10000",
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input} {params.db} --threads {threads}
        """


# Align reads to virus catalog using bowtie2
rule align_reads_to_viruses:
    input:
        R1=R1,
        R2=R2,
        db=results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog.1.bt2",
    output:
        bam=results + "12_VIRUS_ANALYSIS/01_align_viruses/bam_files/{sample}.bam",
        log=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.log",
    log:
        results + "00_LOGS/12_VIRUS_ANALYSIS/align_reads_to_viruses_{sample}.log",
    params:
        db=results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog",
        sam=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.sam",
    conda:
        "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    threads: 8
    benchmark:
        "benchmark/07_VIRUS_ABUNDANCE/align_reads_to_viruses_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {params.sam} > {output.log} 2>&1


        # convert sam to bam
        samtools view -S -b {params.sam} > {output.bam}
        rm {params.sam}
        """


# -----------------------------------------------------
# 02 Metapop
# -----------------------------------------------------
rule read_counts:
    input:
        results + "12_VIRUS_ANALYSIS/01_align_viruses/bam_files/{sample}.bam",
    output:
        results + "12_VIRUS_ANALYSIS/02_metapop/{sample}_read_counts.tsv",
    conda:
        "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/read_counts_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        echo -e {wildcards.sample}"\t"$(samtools view -c {input}) > {output}
        """


rule combine_read_counts_across_samples:
    input:
        expand(
            results + "12_VIRUS_ANALYSIS/02_metapop/{sample}_read_counts.tsv",
            sample=samples,
        ),
    output:
        results + "12_VIRUS_ANALYSIS/02_metapop/combined_read_counts.tsv",
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/combine_read_counts_across_samples.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        cat {input} > {output}
        """


# determine which viruses are present in the sample
rule metapop:
    input:
        bam=expand(
            results + "12_VIRUS_ANALYSIS/01_align_viruses/bam_files/{sample}.bam",
            sample=samples,
        ),
        read_counts=results + "12_VIRUS_ANALYSIS/02_metapop/combined_read_counts.tsv",
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        expand(
            results
            + "12_VIRUS_ANALYSIS/02_metapop/MetaPop/03.Breadth_and_Depth/{sample}_breadth_and_depth.tsv",
            sample=samples,
        ),
    log:
        results + "00_LOGS/12_VIRUS_ANALYSIS/metapop.log",
    params:
        bam_dir=results + "12_VIRUS_ANALYSIS/01_align_viruses/bam_files/",
        viruses_dir=results + "07_VIRUS_DIVERSITY/01_votu_clusters/fasta/",
        viruses=results
        + "07_VIRUS_DIVERSITY/01_votu_clusters/fasta/votu_representatives.fna",
        out_dir=results + "12_VIRUS_ANALYSIS/02_metapop/",
        min_id=95,
        min_breadth=50,
        min_depth=1,
        extra_args="",
    conda:
        "../envs/metapop:1.0.2.yml"
    benchmark:
        "benchmark/07_VIRUS_ABUNDANCE/metapop.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: 8
    shell:
        """
        mkdir {params.viruses_dir}
        cp {input.viruses} {params.viruses}

        # run metapop to identify viruses present in samples
        metapop --input_samples {params.bam_dir} \
        --norm {input.read_counts} \
        --reference {params.viruses_dir} \
        --output {params.out_dir} \
        --id_min {params.min_id} \
        --min_cov {params.min_breadth} \
        --min_dep {params.min_depth} \
        {params.extra_args}
        """


# -----------------------------------------------------
# 03 Filter viruses
# -----------------------------------------------------
# extract present viruses
rule combine_breadth_and_depth_across_samples:
    input:
        metapop=expand(
            results
            + "12_VIRUS_ANALYSIS/02_metapop/MetaPop/03.Breadth_and_Depth/{sample}_breadth_and_depth.tsv",
            sample=samples,
        ),
    output:
        results + "12_VIRUS_ANALYSIS/03_filter_viruses/combined_breadth_and_depth.tsv",
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/combine_breadth_and_depth_across_samples.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input.metapop} > {output}
        """


rule extract_present_viruses:
    input:
        viruses=results
        + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_viruses.fna",
        metapop=results
        + "12_VIRUS_ANALYSIS/03_filter_viruses/combined_breadth_and_depth.tsv",
    output:
        results + "12_VIRUS_ANALYSIS/03_filter_viruses/present_viruses.fna",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/extract_present_viruses.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/07_extract_present_viruses.py"
