# -------------------------------------
# Read Assembly Module (If input_data = "reads")
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
        "Assembling {wildcards.sample} using metaSPAdes"
    input:
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz",
    output:
        results + "03_READ_ASSEMBLY/01_spades/{sample}_meta/contigs.fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/01_spades/{sample}_meta/",
        extra_args=config["read_assembly"]["spades_arguments"],
    # conda:
    #     "../envs/spades:3.15.4--h95f258a_0.yml"
    container:
        "docker://quay.io/biocontainers/spades:3.15.4--h95f258a_0"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/metaspades_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="100000",
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
        "Assembling {wildcards.sample} using metaviralSPAdes"
    input:
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz",
    output:
        results + "03_READ_ASSEMBLY/01_spades/{sample}_metaviral/contigs.fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/01_spades/{sample}_metaviral/",
        extra_args=config["read_assembly"]["spades_arguments"],
    # conda:
    #     "../envs/spades:3.15.4--h95f258a_0.yml"
    container:
        "docker://quay.io/biocontainers/spades:3.15.4--h95f258a_0"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/metaviralspades_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="100000",
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
        """


# Assemble reads using rnaviralspades
rule rnaviralspades:
    message:
        "Assembling {wildcards.sample} using rnaviralSPAdes"
    input:
        R1=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz",
        R2=results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz",
    output:
        results + "03_READ_ASSEMBLY/01_spades/{sample}_rnaviral/contigs.fasta",
    params:
        output_dir=results + "03_READ_ASSEMBLY/01_spades/{sample}_rnaviral/",
        extra_args=config["read_assembly"]["spades_arguments"],
    # conda:
    #     "../envs/spades:3.15.4--h95f258a_0.yml"
    container:
        "docker://quay.io/biocontainers/spades:3.15.4--h95f258a_0"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/rnaviralspades_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="100000",
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
        """


# combine assemblies if multiple were performed
assemblies = []
if "meta" in config["read_assembly"]["assembly_modes"]:
    assemblies.append(
        results + "03_READ_ASSEMBLY/01_spades/{sample}_meta/contigs.fasta"
    )
if "metaviral" in config["read_assembly"]["assembly_modes"]:
    assemblies.append(
        results + "03_READ_ASSEMBLY/01_spades/{sample}_metaviral/contigs.fasta"
    )
if "rnaviral" in config["read_assembly"]["assembly_modes"]:
    assemblies.append(
        results + "03_READ_ASSEMBLY/01_spades/{sample}_rnaviral/contigs.fasta"
    )


# combine all spades assembly types
rule combine_spades_assemblies:
    message:
        "Combining {wildcards.sample} assemblies from different assemblers"
    input:
        assemblies,
    output:
        results + "03_READ_ASSEMBLY/01_spades/{sample}_contigs.fasta",
    benchmark:
        "benchmark/03_READ_ASSEMBLY/combine_spades_assemblies_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # combine assemblies from different assemblers
        cat {input} > {output}
        """


# -----------------------------------------------------
# 02 Filter and dereplicate contigs within samples
# -----------------------------------------------------
# filter contigs based on contig length
rule contig_length_filter:
    message:
        "Filtering {wildcards.sample} assemblies to only those longer than {params.min_length}"
    input:
        results + "03_READ_ASSEMBLY/01_spades/{sample}_contigs.fasta",
    output:
        results + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs.fasta",
    params:
        min_length=config["read_assembly"]["min_contig_length"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/contig_length_filter_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/03_contig_length_filter.py"


# blast contigs against one another for dereplication
rule blast_contigs_within_samples:
    message:
        "BLASTing {wildcards.sample} contigs against one another to dererplicate dataset"
    input:
        results + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs.fasta",
    output:
        results + "03_READ_ASSEMBLY/02_contig_filters/{sample}/contigs_blast.tsv",
    params:
        blastdb=results + "03_READ_ASSEMBLY/02_contig_filters/{sample}/contigs_blastdb",
        blast_tsv=results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/contigs_blast.tsv",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/blast_contigs_within_samples_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    threads: config["virus_dereplication"]["blast_threads"]
    shell:
        """
        # make a blast db from contigs
        makeblastdb -in {input} -out {params.blastdb} -dbtype nucl

        # all against all blast
        blastn -query {input} -db {params.blastdb} -out {output} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
        """


# download build mgv repo and HMM files
rule build_mgv:
    message:
        "Cloning MGV repository to obtain dereplication script"
    output:
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
        cluster_script=resources + "mgv/ani_cluster/cluster.py",
        amino_acid_script=resources + "mgv/aai_cluster/amino_acid_identity.py",
    params:
        mgv_dir=resources + "mgv",
    benchmark:
        "benchmark/03_READ_ASSEMBLY/build_mgv.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # clone MGV repository
        rm -rf {params.mgv_dir}
        git clone https://github.com/snayfach/MGV.git {params.mgv_dir}
        """


# keep only one species representative from each sample
rule dereplicate_contigs_within_samples:
    message:
        "Dereplicating {wildcards.sample} contigs using BLAST results"
    input:
        viruses=results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs.fasta",
        blast=results + "03_READ_ASSEMBLY/02_contig_filters/{sample}/contigs_blast.tsv",
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
        cluster_script=resources + "mgv/ani_cluster/cluster.py",
    output:
        results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs_clusters.tsv",
    params:
        ani_tsv=results + "03_READ_ASSEMBLY/02_contig_filters/{sample}/contigs_ani.tsv",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/dereplicate_contigs_within_samples_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # calculate ani and af from blast results
        python {input.blastani_script} -i {input.blast} -o {params.ani_tsv}

        # cluster contigs based on 95% ANI and 85% AF
        python {input.cluster_script} --fna {input.viruses} --ani {params.ani_tsv} --out {output} --min_ani 95 --min_qcov 85 --min_tcov 0
        """


# extract dereplicated contigs
rule extract_dereplicated_contigs_within_samples:
    message:
        "Extracting {wildcards.sample} representatives to retain only one of each replicate"
    input:
        clusters=results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs_clusters.tsv",
        viruses=results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs.fasta",
    output:
        results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs_dereplicated.fasta",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/extract_dereplicated_contigs_within_samples_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/03_extract_clustered_viruses.py"


# -----------------------------------------------------
# 03 QUAST
# -----------------------------------------------------
# run quast to determine the quality of assemblies
rule quast:
    message:
        "Running QUAST on {wildcards.sample} to determine assembly quality"
    input:
        results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs_dereplicated.fasta",
    output:
        results + "03_READ_ASSEMBLY/03_quast/{sample}/transposed_report.tsv",
    params:
        output_dir=results + "03_READ_ASSEMBLY/03_quast/{sample}",
        min_len=config["read_assembly"]["min_contig_length"],
        labels="{sample}",
        extra_args=config["read_assembly"]["quast_arguments"],
    # conda:
    #     "../envs/quast:5.0.2--py27pl5321h8eb80aa_6.yml"
    container:
        "docker://quay.io/biocontainers/quast:5.0.2--py27pl5321h8eb80aa_6"
    benchmark:
        "benchmark/03_READ_ASSEMBLY/quast_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
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
            results + "03_READ_ASSEMBLY/03_quast/{sample}/transposed_report.tsv",
            sample=samples,
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
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # Generate MULTIQC report from QUAST results
        multiqc {params.quast_dir} \
        -o {params.quast_dir} -f

        # move report to final destination
        mv {params.quast_out} {output}
        """
