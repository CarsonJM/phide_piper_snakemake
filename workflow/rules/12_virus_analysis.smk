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
        sam=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.sam",
        log=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.log",
    params:
        db=results + "12_VIRUS_ANALYSIS/01_align_viruses/virus_catalog",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/align_reads_to_viruses_{sample}.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="10000",
    threads: config["virus_analysis"]["bowtie2_threads"]
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {output.sam} > {output.log} 2>&1
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
# 02 inStrain
# -----------------------------------------------------
# Run instrain to preprocess and analyze alignments
rule make_stb_file:
    message:
        "Making scaffold to bin file for inStrain"
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        results + "12_VIRUS_ANALYSIS/02_instrain/stb_file.tsv",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/make_stb_file.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    script:
        "../scripts/12_generate_stb_file.py"


# Run instrain to preprocess and analyze alignments
rule instrain_profile:
    message:
        "Running inStrain to preprocess alignments and calculate diversity metrics"
    input:
        sam=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.sam",
        fasta=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
        genes=results
        + "07_VIRUS_DIVERSITY/02_vcontact2/votu_representatives_proteins.fna",
        stb=results + "12_VIRUS_ANALYSIS/02_instrain/stb_file.tsv",
    output:
        sorted_bam=results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.sorted.bam",
        instrain=results
        + "12_VIRUS_ANALYSIS/02_instrain/{sample}/output/{sample}_genome_info.tsv",
    params:
        out_dir=results + "12_VIRUS_ANALYSIS/02_instrain/{sample}",
        min_id=config["virus_analysis"]["min_id"],
        min_breadth=config["virus_analysis"]["min_breadth"],
        min_depth=config["virus_analysis"]["min_depth"],
        extra_args=config["virus_analysis"]["instrain_profile_arguments"],
    container:
        "docker://quay.io/biocontainers/instrain:1.5.7--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/instrain_profile_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["virus_analysis"]["instrain_threads"]
    shell:
        """
        # run instrain profile
        inStrain profile \
        {input.sam} \
        {input.fasta} \
        --output {params.out_dir} \
        --processes {threads} \
        --min_read_ani {params.min_id} \
        --min_cov {params.min_depth} \
        --gene_file {input.genes} \
        --stb {input.stb} \
        --skip_plot_generation \
        {params.extra_args}
        """


# Run instrain to preprocess and analyze alignments
rule instrain_compare:
    message:
        "Running inStrain to preprocess alignments and calculate diversity metrics"
    input:
        profiles=expand(
            results
            + "12_VIRUS_ANALYSIS/02_instrain/{sample}/output/{sample}_genome_info.tsv",
            sample=samples,
        ),
        sam=expand(
            results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.sam",
            sample=samples,
        ),
        stb=results + "12_VIRUS_ANALYSIS/02_instrain/stb_file.tsv",
    output:
        results
        + "12_VIRUS_ANALYSIS/02_instrain/compare/output/compare_genomeWide_compare.tsv",
    params:
        bams_dir=results + "12_VIRUS_ANALYSIS/01_align_viruses/",
        profiles=expand(
            results + "12_VIRUS_ANALYSIS/02_instrain/{sample}/",
            sample=samples,
        ),
        out_dir=results + "12_VIRUS_ANALYSIS/02_instrain/compare",
        min_id=config["virus_analysis"]["min_id"],
        min_breadth=config["virus_analysis"]["min_breadth"],
        min_depth=config["virus_analysis"]["min_depth"],
        extra_args=config["virus_analysis"]["instrain_compare_arguments"],
    container:
        "docker://quay.io/biocontainers/instrain:1.5.7--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/instrain_compare.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["virus_analysis"]["instrain_threads"]
    shell:
        """
        # run instrain profile
        inStrain compare \
        --input {params.profiles} \
        --output {params.out_dir} \
        --processes {threads} \
        --stb {input.stb} \
        --breadth {params.min_breadth} \
        --ani_threshold {params.min_id} \
        --coverage_treshold {params.min_breadth} \
        --skip_plot_generation \
        {params.extra_args}
        """


# -----------------------------------------------------
# 03 CoverM
# -----------------------------------------------------
# Run instrain to preprocess and analyze alignments
rule coverm:
    message:
        "Running inStrain to preprocess alignments and calculate diversity metrics"
    input:
        expand(
            results + "12_VIRUS_ANALYSIS/01_align_viruses/{sample}.sorted.bam",
            sample=samples,
        ),
    output:
        results
        + "12_VIRUS_ANALYSIS/03_coverm/"
        + config["virus_analysis"]["coverm_method"],
    params:
        min_id=config["virus_analysis"]["min_id"],
        method=config["virus_analysis"]["coverm_method"],
        extra_args=config["virus_analysis"]["coverm_arguments"],
    container:
        "docker://quay.io/biocontainers/coverm:0.6.1--h8b24846_3"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/coverm.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["virus_analysis"]["coverm_threads"]
    shell:
        """
        # run coverm
        coverm contig \
        --bam-files {input} \
        --min-read-percent-identity {params.min_id} \
        --methods {params.method} \
        --output-file {output} \
        --threads {threads} \
        --verbose \
        {params.extra_args}
        """
