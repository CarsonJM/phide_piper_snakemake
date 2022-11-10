# -----------------------------------------------------
# Virus abundance Module (will always run)
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
# read symlink in 04_virus_identification.smk
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
    message:
        "Building a bowtie2 db of vOTU representative viruses"
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog.1.bt2",
    params:
        db=results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/build_viruses_bowtie2db.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["virus_analysis"]["bowtie2_threads"]
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input} {params.db} --threads {threads}
        """


# Align reads to virus catalog using bowtie2
rule align_reads_to_viruses:
    message:
        "Aligning reads to vOTU database to determine virus abundances"
    input:
        R1=R1,
        R2=R2,
        db=results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog.1.bt2",
    output:
        bam=results + "12_VIRUS_ANALYSIS/01_align_viruses/bam_files/{sample}.bam",
        log=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.log",
    params:
        db=results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog",
        sam=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.sam",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/align_reads_to_viruses_{sample}.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="10000",
        partition="compute-hugemem",
    threads: config["virus_analysis"]["bowtie2_threads"]
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


# generate report for alignments
rule bowtie2_multiqc:
    message:
        "Running MULTIQC to visualize bowtie2 alignment rates"
    input:
        expand(
            results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.log",
            sample=samples,
        ),
    output:
        report(
            results + "12_VIRUS_ANALYSIS/bowtie2_multiqc.html",
            category="Step 12: Virus analysis",
        ),
    params:
        bt2_input=results + "12_VIRUS_ANALYSIS/01_align_viruses/*.log",
        bt2_dir=results + "12_VIRUS_ANALYSIS/01_align_viruses/",
        bt2_out=results + "12_VIRUS_ANALYSIS/01_align_viruses/multiqc_report.html",
    # conda:
    #     "../envs/multiqc:1.12--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/bowtie2_multiqc.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # run multiqc on bowtie2 outputs
        multiqc {params.bt2_input} \
        -o {params.bt2_dir} -f

        mv {params.bt2_out} {output}
        """


# -----------------------------------------------------
# 02 Metapop
# -----------------------------------------------------
# generate read counts using bam files
rule read_counts:
    message:
        "Generating read counts for {sample} FastQ"
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
        # print read counts out to file
        echo -e {wildcards.sample}"\t"$(samtools view -c {input}) > {output}
        """


# combine read counts files across all sampels
rule combine_read_counts_across_samples:
    message:
        "Combining read counts files across all samples"
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
        # combine read counts files
        cat {input} > {output}
        """


# Run metapop to preprocess and analyze alignments
rule metapop:
    message:
        "Running MetaPop to preprocess and analyze virus alignments"
    input:
        bam=expand(
            results + "12_VIRUS_ANALYSIS/01_align_viruses/bam_files/{sample}.bam",
            sample=samples,
        ),
        read_counts=results + "12_VIRUS_ANALYSIS/02_metapop/combined_read_counts.tsv",
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        breadth=expand(
            results
            + "12_VIRUS_ANALYSIS/02_metapop/MetaPop/03.Breadth_and_Depth/{sample}_breadth_and_depth.tsv",
            sample=samples,
        ),
        pdf=report(
            results
            + "12_VIRUS_ANALYSIS/02_metapop/MetaPop/12.Visualizations/preprocessing_summaries.pdf",
            category="Step 12: Virus analysis",
        ),
    params:
        bam_dir=results + "12_VIRUS_ANALYSIS/01_align_viruses/bam_files/",
        viruses_dir=results + "07_VIRUS_DIVERSITY/01_votu_clusters/fasta/",
        viruses=results
        + "07_VIRUS_DIVERSITY/01_votu_clusters/fasta/votu_representatives.fna",
        out_dir=results + "12_VIRUS_ANALYSIS/02_metapop/",
        min_id=config["virus_analysis"]["min_id"],
        min_breadth=config["virus_analysis"]["min_breadth"],
        min_depth=config["virus_analysis"]["min_depth"],
        extra_args=config["virus_analysis"]["metapop_arguments"],
    conda:
        "../envs/metapop:1.0.2.yml"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/metapop.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["virus_analysis"]["metapop_threads"]
    shell:
        """
        rm -rf {params.viruses_dir}
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
        --threads {threads} \
        {params.extra_args}
        """
