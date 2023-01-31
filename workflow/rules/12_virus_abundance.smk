# -----------------------------------------------------
# Virus abundance Module (will always run)
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["samples_df"], sep="\t")
group_sample_dictionary = samples_df.groupby(["group"])["sample"].apply(list).to_dict()

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
    combine_instrain_profiles,
    identify_integrated_prophages,
    build_propagate,


# -----------------------------------------------------
# 00 Determine inputs for module
# -----------------------------------------------------
if (
    config["input_data"] == "reads"
    or config["input_data"] == "contigs"
    or config["input_data"] == "vls"
):
    viruses = (
        results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/filtered_viruses.fasta",
    )
elif config["input_data"] == "viruses":
    viruses = results + "00_INPUT/{group_sample}_viruses.fasta"


# read symlink in 04_virus_identification.smk
if config["input_data"] == "reads":
    R1 = results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz"
    R2 = results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz"
elif (
    config["input_data"] == "contigs"
    or config["input_data"] == "viruses"
    or config["input_data"] == "processed_viruses"
):
    R1 = results + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R1.fastq.gz"
    R2 = results + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R2.fastq.gz"


# -----------------------------------------------------
# 01 Create genome catalog to align to
# -----------------------------------------------------
# Align reads to virus catalog using bowtie2
rule build_viruses_bowtie2db:
    message:
        "Building a bowtie2 db of vOTU representative viruses"
    input:
        viruses
    output:
        results + "12_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog.1.bt2",
    params:
        db=results + "12_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog",
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/build_viruses_bowtie2db.tsv"
    resources:
        runtime=config["virus_analysis"]["bowtie2_runtime"],
        mem_mb=config["virus_analysis"]["bowtie2_memory"],
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
        db=results + "12_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog.1.bt2",
    output:
        sam=results + "12_VIRUS_ABUNDANCE/01_align_viruses/{group_sample}.sam",
        log=results + "12_VIRUS_ABUNDANCE/01_align_viruses/{group_sample}.log",
    params:
        db=results + "12_VIRUS_ABUNDANCE/01_align_viruses/virus_catalog",
        extra_args=config["virus_analysis"]["bowtie2_arguments"],
    # conda:
    #     "../envs/kneaddata:0.10.0--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/kneaddata:0.10.0--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/align_reads_to_viruses_{group_sample}.tsv"
    resources:
        runtime=config["virus_analysis"]["bowtie2_runtime"],
        mem_mb=config["virus_analysis"]["bowtie2_memory"],
    threads: config["virus_analysis"]["bowtie2_threads"]
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        {params.extra_args} \
        -S {output.sam} > {output.log} 2>&1
        """


# generate report for alignments
rule bowtie2_multiqc:
    message:
        "Running MULTIQC to visualize bowtie2 alignment rates"
    input:
        expand(
            results + "12_VIRUS_ABUNDANCE/01_align_viruses/{group_sample}.log",
            group_sample=groups_samples,
        ),
    output:
        report(
            results + "12_VIRUS_ABUNDANCE/bowtie2_multiqc.html",
            category="Step 12: Virus analysis",
        ),
    params:
        bt2_input=results + "12_VIRUS_ABUNDANCE/01_align_viruses/*.log",
        bt2_dir=results + "12_VIRUS_ABUNDANCE/01_align_viruses/",
        bt2_out=results + "12_VIRUS_ABUNDANCE/01_align_viruses/multiqc_report.html",
    # conda:
    #     "../envs/multiqc:1.12--pyhdfd78af_0.yml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/bowtie2_multiqc.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
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
        viruses
    output:
        results + "12_VIRUS_ABUNDANCE/02_instrain_profile/stb_file.tsv",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/make_stb_file.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/12_generate_stb_file.py"


# prodigal-gv
rule prodigal_gv:
    input:
        viruses,
    output:
        faa=results + "12_VIRUS_ABUNDANCE/02_instrain_profile/instrain_proteins.faa",
        fna=results + "12_VIRUS_ABUNDANCE/02_instrain_profile/instrain_proteins.fna",
    params:
        extra_args=config['virus_diversity']['prodigal_gv_arguments']
    conda:
        "../envs/prodigal_gv.yml"
    benchmark:
        "benchmark/12_VIRUS_ANALYSIS/prodigal_gv.tsv"
    resources:
        runtime=config["virus_diversity"]["prodigal_gv_runtime"],
        mem_mb=config["virus_diversity"]["prodigal_gv_memory"],
    shell:
        """
        # create gene2genome file for vcontact2
        prodigal-gv \
        -p meta \
        -i {input} \
        -d {output.fna} \
        -a {output.faa} \
        {params.extra_args}
        """


# Run instrain to preprocess and analyze alignments
rule instrain_profile:
    message:
        "Running inStrain to preprocess alignments and calculate diversity metrics"
    input:
        sam=results + "12_VIRUS_ABUNDANCE/01_align_viruses/{group_sample}.sam",
        fasta=viruses,
        genes=results + "12_VIRUS_ABUNDANCE/02_instrain_profile/instrain_proteins.fna",
        stb=results + "12_VIRUS_ABUNDANCE/02_instrain_profile/stb_file.tsv",
    output:
        sorted_bam=results
        + "12_VIRUS_ABUNDANCE/01_align_viruses/{group_sample}.sorted.bam",
        genome=results
        + "12_VIRUS_ABUNDANCE/02_instrain_profile/{group_sample}/output/{group_sample}_genome_info.tsv",
        gene=results
        + "12_VIRUS_ABUNDANCE/02_instrain_profile/{group_sample}/output/{group_sample}_gene_info.tsv",
    params:
        out_dir=results + "12_VIRUS_ABUNDANCE/02_instrain_profile/{group_sample}",
        min_id=config["virus_analysis"]["min_id"],
        min_breadth=config["virus_analysis"]["min_breadth"],
        min_depth=config["virus_analysis"]["min_depth"],
        extra_args=config["virus_analysis"]["instrain_profile_arguments"],
    container:
        "docker://quay.io/biocontainers/instrain:1.5.7--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/instrain_profile_{group_sample}.tsv"
    resources:
        runtime=config["virus_analysis"]["instrain_runtime"],
        mem_mb=config["virus_analysis"]["instrain_memory"],
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
        {params.extra_args}

        # add sample column to each instrain output
        s={wildcards.group_sample}
        sed -i "s/$/\t$s/" {output.genome}
        sample="sample"
        sed -i "1s/$s/$sample/" {output.genome}
        """


# combine instrain reports across samples
rule compute_virus_abundances:
    message:
        "Computing virus abundances using inStrain outputs"
    input:
        results
        + "12_VIRUS_ABUNDANCE/02_instrain_profile/{group_sample}/output/{group_sample}_genome_info.tsv",
    output:
        results
        + "12_VIRUS_ABUNDANCE/02_instrain_profile/{group_sample}/output/{group_sample}_genome_info_w_abundance.tsv",
    params:
        min_breadth=config["virus_analysis"]["min_breadth"],
        recover_low_abundance=config["virus_analysis"]["recover_low_abundance"],
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/compute_virus_abundances_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/12_compute_virus_abundance.py"


# combine instrain reports across samples
rule combine_instrain_profiles:
    message:
        "Combining inStrain reports across samples"
    input:
        expand(
            results
            + "12_VIRUS_ABUNDANCE/02_instrain_profile/{group_sample}/output/{group_sample}_genome_info_w_abundance.tsv",
            group_sample=groups_samples,
        ),
    output:
        results
        + "12_VIRUS_ABUNDANCE/02_instrain_profile/combined_genome_info_w_abundance.tsv",
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/combine_instrain_profiles.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/12_combine_instrain_profiles.py"


# identify replicate reads
group_sample_dictionary = samples_df.groupby(["group"])["sample"].apply(list).to_dict()


# Run instrain to preprocess and analyze alignments
rule instrain_compare:
    message:
        "Running inStrain to compare samples within {wildcards.group}"
    input:
        profiles=lambda wildcards: expand(
            results
            + "12_VIRUS_ABUNDANCE/02_instrain_profile/{{group}}_{sample}/output/{{group}}_{sample}_genome_info.tsv",
            sample=group_sample_dictionary[wildcards.group],
        ),
        sam=lambda wildcards: expand(
            results + "12_VIRUS_ABUNDANCE/01_align_viruses/{{group}}_{sample}.sam",
            sample=group_sample_dictionary[wildcards.group],
        ),
        stb=results + "12_VIRUS_ABUNDANCE/02_instrain_profile/stb_file.tsv",
    output:
        results
        + "12_VIRUS_ABUNDANCE/03_instrain_compare/{group}/output/{group}_genomeWide_compare.tsv",
    params:
        bams_dir=results + "12_VIRUS_ABUNDANCE/01_align_viruses/",
        profiles=lambda wildcards, input: [
            file.rpartition("/")[0].rpartition("/")[0] + "/" for file in input.profiles
        ],
        out_dir=results + "12_VIRUS_ABUNDANCE/03_instrain_compare/{group}",
        min_id=config["virus_analysis"]["min_id"],
        min_breadth=config["virus_analysis"]["min_breadth"],
        min_depth=config["virus_analysis"]["min_depth"],
        extra_args=config["virus_analysis"]["instrain_compare_arguments"],
    container:
        "docker://quay.io/biocontainers/instrain:1.5.7--pyhdfd78af_0"
    benchmark:
        "benchmark/12_VIRUS_ABUNDANCE/instrain_compare_{group}.tsv"
    resources:
        runtime=config["virus_analysis"]["instrain_runtime"],
        mem_mb=config["virus_analysis"]["instrain_memory"],
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
        {params.extra_args}
        """
