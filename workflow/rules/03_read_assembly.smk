# -------------------------------------
# Read Assembly Module (If input_data = "reads")
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
# Read assembly rules
# -------------------------------------
localrules:
    combine_spades_assemblies,
    combine_quast_across_samples,


# -----------------------------------------------------
# 01 SPAdes
# -----------------------------------------------------
# Assemble reads using metaspades
rule metaspades:
    message:
        "Assembling {wildcards.group_sample} using metaSPAdes"
    input:
        R1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz",
    output:
        results + "03_READ_ASSEMBLY/01_spades/meta/{group_sample}/contigs.fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/01_spades/meta/{group_sample}",
        extra_args=config["read_assembly"]["spades_arguments"],
    # conda:
    #     "../envs/spades:3.15.4--h95f258a_0.yml"
    container:
        "docker://quay.io/biocontainers/spades:3.15.4--h95f258a_0"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/metaspades_{group_sample}.tsv"
    resources:
        runtime=config["read_assembly"]["spades_runtime"],
        mem_mb=config["read_assembly"]["spades_memory"],
    threads: config["read_assembly"]["spades_threads"]
    shell:
        """
        # assemble reads using spades
        spades.py \
        --meta \
        -1 {input.R1} \
        -2 {input.R2} \
        -o {params.output_dir} \
        --threads {threads} \
        {params.extra_args}
        """


# Assemble reads using metaviralspades
rule metaviralspades:
    message:
        "Assembling {wildcards.group_sample} using metaviralSPAdes"
    input:
        R1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz",
    output:
        results + "03_READ_ASSEMBLY/01_spades/metaviral/{group_sample}/contigs.fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/01_spades/metaviral/{group_sample}",
        extra_args=config["read_assembly"]["spades_arguments"],
    # conda:
    #     "../envs/spades:3.15.4--h95f258a_0.yml"
    container:
        "docker://quay.io/biocontainers/spades:3.15.4--h95f258a_0"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/metaviralspades_{group_sample}.tsv"
    resources:
        runtime=config["read_assembly"]["spades_runtime"],
        mem_mb=config["read_assembly"]["spades_memory"],
    threads: config["read_assembly"]["spades_threads"]
    shell:
        """
        # assemble reads using metaviralspades
        spades.py \
        --metaviral \
        -1 {input.R1} \
        -2 {input.R2} \
        -o {params.output_dir} \
        --threads {threads} \
        {params.extra_args}

        touch {output}
        """


# Assemble reads using rnaviralspades
rule rnaviralspades:
    message:
        "Assembling {wildcards.group_sample} using rnaviralSPAdes"
    input:
        R1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz",
    output:
        results + "03_READ_ASSEMBLY/01_spades/rnaviral/{group_sample}/contigs.fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/01_spades/rnaviral/{group_sample}",
        extra_args=config["read_assembly"]["spades_arguments"],
    # conda:
    #     "../envs/spades:3.15.4--h95f258a_0.yml"
    container:
        "docker://quay.io/biocontainers/spades:3.15.4--h95f258a_0"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/rnaviralspades_{group_sample}.tsv"
    resources:
        runtime=config["read_assembly"]["spades_runtime"],
        mem_mb=config["read_assembly"]["spades_memory"],
    threads: config["read_assembly"]["spades_threads"]
    shell:
        """
        # assemble reads using rnaviralspades
        spades.py \
        --rnaviral \
        -1 {input.R1} \
        -2 {input.R2} \
        -o {params.output_dir} \
        --threads {threads} \
        {params.extra_args}

        touch {output}
        """


# combine assemblies if multiple were performed
assemblies = []
if "meta" in config["read_assembly"]["assembly_modes"]:
    assemblies.append(
        results + "03_READ_ASSEMBLY/01_spades/meta/{group_sample}/contigs.fasta"
    )
if "metaviral" in config["read_assembly"]["assembly_modes"]:
    assemblies.append(
        results + "03_READ_ASSEMBLY/01_spades/metaviral/{group_sample}/contigs.fasta"
    )
if "rnaviral" in config["read_assembly"]["assembly_modes"]:
    assemblies.append(
        results + "03_READ_ASSEMBLY/01_spades/rnaviral/{group_sample}/contigs.fasta"
    )


# combine all spades assembly types
rule combine_spades_assemblies:
    message:
        "Combining {wildcards.group_sample} assemblies from different assemblers"
    input:
        assemblies,
    output:
        results + "03_READ_ASSEMBLY/01_spades/{group_sample}_contigs.fasta",
    params:
        assemblies=results + "03_READ_ASSEMBLY/01_spades/*/*",
    benchmark:
        "benchmark/03_READ_ASSEMBLY/combine_spades_assemblies_{group_sample}.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # combine assemblies from different assemblers
        cat {input} > {output}

        # rm -rf {params.assemblies}
        """


# -----------------------------------------------------
# 02 Filter and dereplicate contigs within samples
# -----------------------------------------------------
# filter contigs based on contig length
rule contig_length_filter:
    message:
        "Filtering {wildcards.group_sample} assemblies to only those longer than {params.min_length}"
    input:
        results + "03_READ_ASSEMBLY/01_spades/{group_sample}_contigs.fasta",
    output:
        results
        + "03_READ_ASSEMBLY/02_contig_filters/{group_sample}/{group_sample}_contigs.fasta",
    params:
        min_length=config["read_assembly"]["min_contig_length"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/contig_length_filter_{group_sample}.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    script:
        "../scripts/03_contig_length_filter.py"


# -----------------------------------------------------
# 03 QUAST
# -----------------------------------------------------
# run quast to determine the quality of assemblies
rule quast:
    message:
        "Running QUAST on {wildcards.group_sample} to determine assembly quality"
    input:
        results
        + "03_READ_ASSEMBLY/02_contig_filters/{group_sample}/{group_sample}_contigs.fasta",
    output:
        results + "03_READ_ASSEMBLY/03_quast/{group_sample}/transposed_report.tsv",
    params:
        output_dir=results + "03_READ_ASSEMBLY/03_quast/{group_sample}",
        min_len=config["read_assembly"]["min_contig_length"],
        labels="{group_sample}",
        extra_args=config["read_assembly"]["quast_arguments"],
    # conda:
    #     "../envs/quast:5.0.2--py27pl5321h8eb80aa_6.yml"
    container:
        "docker://quay.io/biocontainers/quast:5.0.2--py27pl5321h8eb80aa_6"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/quast_{group_sample}.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # assembly analysis using quast
        metaquast.py \
        {input} \
        -o {params.output_dir} \
        --threads {threads} \
        --min-contig {params.min_len} \
        --contig-thresholds 0,1000,5000,10000,{params.min_len} \
        --labels {params.labels} \
        {params.extra_args}
        """


# generate report for assemblies
rule quast_multiqc:
    message:
        "Generating MULTIQC report for QUAST analyses"
    input:
        expand(
            results + "03_READ_ASSEMBLY/03_quast/{group_sample}/transposed_report.tsv",
            group_sample=groups_samples,
        ),
    output:
        report(
            results + "03_READ_ASSEMBLY/quast_multiqc_report.html",
            category="Step 03: Read assembly",
        ),
    params:
        quast_input=results + "03_READ_ASSEMBLY/03_quast/*/report.tsv",
        quast_dir=results + "03_READ_ASSEMBLY/03_quast/",
        quast_out=results + "03_READ_ASSEMBLY/03_quast/multiqc_report.html",
    # conda:
    #     "../envs/multiqc:1.12--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/quast_multiqc.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # Generate MULTIQC report from QUAST results
        multiqc {params.quast_dir} \
        -o {params.quast_dir} -f

        # move report to final destination
        mv {params.quast_out} {output}
        """


# combine quast reports
rule combine_quast_across_samples:
    message:
        "Combining QUAST reports"
    input:
        expand(
            results + "03_READ_ASSEMBLY/03_quast/{group_sample}/transposed_report.tsv",
            group_sample=groups_samples,
        ),
    output:
        results + "03_READ_ASSEMBLY/assembly_report.tsv",
    params:
        quast_dirs=results + "03_READ_ASSEMBLY/03_quast/",
    benchmark:
        "benchmark/03_READ_ASSEMBLY/combine_quast.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}

        rm -rf {params.quast_dirs}
        """
