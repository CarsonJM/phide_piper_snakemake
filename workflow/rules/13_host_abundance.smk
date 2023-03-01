# -----------------------------------------------------
# Host abundance Module
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["samples_df"], sep="\t")


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "../report/workflow.rst"


# -----------------------------------------------------
# Host abundance rules
# -----------------------------------------------------


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
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps_viruses.fasta",
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
    R1 = results + "00_INPUT/{group_sample}_proprocessed_1.fastq.gz"
    R2 = results + "00_INPUT/{group_sample}_preprocessed_2.fastq.gz"


# -----------------------------------------------------
# 01 MetaPhlan
# -----------------------------------------------------
# download metaphlan database
rule download_metaphlan_db:
    message:
        "Downloading MetaPhlan4 database mpa_vJan21_CHOCOPhlAnSGB_202103"
    output:
        bowtie2_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2l",
        spa_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.pkl",
    params:
        mpa_dir=resources + "metaphlan/",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/13_HOST_ABUNDANCE/download_metaphlan_db.tsv"
    resources:
        runtime="4h",
        mem_mb="100GB",
    threads: config["host_abundance"]["metaphlan_threads"]
    shell:
        """
        rm -rf {params.mpa_dir}

        # download metaphlan database
        metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103 \
        --nproc {threads} \
        --bowtie2db {params.mpa_dir}
        """


# run metaphlan
rule metaphlan:
    message:
        "Running MetaPhlan4 on {wildcards.group_sample}"
    input:
        mpa_db=resources + "metaphlan/mpa_vJan21_CHOCOPhlAnSGB_202103.1.bt2l",
        R1=R1,
        R2=R2,
    output:
        mpa=results + "13_HOST_ABUNDANCE/01_metaphlan/{group_sample}_profile.tsv",
        sam=results + "13_HOST_ABUNDANCE/01_metaphlan/{group_sample}.sam.bz2",
    params:
        extra_args=config["host_abundance"]["metaphlan_arguments"],
        mpa_dir=resources + "metaphlan/",
        index="mpa_vJan21_CHOCOPhlAnSGB_202103",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/13_HOST_ABUNDANCE/metaphlan_{group_sample}.tsv"
    resources:
        runtime=config["host_abundance"]["metaphlan_runtime"],
        mem_mb=config["host_abundance"]["metaphlan_memory"],
    threads: config["host_abundance"]["metaphlan_threads"]
    shell:
        """
        # run metaphlan
        metaphlan {input.R1},{input.R2} \
        --input_type fastq \
        --nproc {threads} \
        --index {params.index} \
        --bowtie2db {params.mpa_dir} \
        --no_map \
        --samout {output.sam} \
        --output_file {output.mpa} \
        {params.extra_args}
        """


# run convert SGB to GTDB taxonomy
rule sgb_to_gtdb_taxonomy:
    message:
        "Converting MetaPhlan4 SGB taxonomy to standardized GTDB taxonomy for {wildcards.group_sample}"
    input:
        results + "13_HOST_ABUNDANCE/01_metaphlan/{group_sample}_profile.tsv",
    output:
        results + "13_HOST_ABUNDANCE/01_metaphlan/{group_sample}_profile_gtdb.tsv",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/13_HOST_ABUNDANCE/sgb_to_gtdb_taxonomy_{group_sample}.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # convert sgb to gtdb taxonomy
        sgb_to_gtdb_profile.py \
        --input {input} \
        --output {output}
        """


# combine metaphlan profils across samples
rule merge_metaphlan_profiles:
    message:
        "Combining MetaPhlan4 results across all samples"
    input:
        sgb=expand(
            results + "13_HOST_ABUNDANCE/01_metaphlan/{group_sample}_profile.tsv",
            group_sample=groups_samples,
        ),
        gtdb=expand(
            results + "13_HOST_ABUNDANCE/01_metaphlan/{group_sample}_profile_gtdb.tsv",
            group_sample=groups_samples,
        ),
    output:
        sgb=results + "13_HOST_ABUNDANCE/metaphlan_merged_profiles.tsv",
        gtdb=results + "13_HOST_ABUNDANCE/metaphlan_merged_profiles_gtdb.tsv",
    params:
        in_dir=results + "13_HOST_ABUNDANCE/01_metaphlan/",
    # conda:
    #     "../envs/humann:3.6--pyh7cba7a3_1.yml"
    container:
        "docker://quay.io/biocontainers/humann:3.6--pyh7cba7a3_1"
    benchmark:
        "benchmark/13_HOST_ABUNDANCE/merge_metaphlan_profiles.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # merge metaphlan profiles
        merge_metaphlan_tables.py {input.sgb} > {output.sgb}

        # merge metaphlan gtdb profiles
        merge_metaphlan_tables.py --gtdb_profiles {input.gtdb} > {output.gtdb}
        """


# -----------------------------------------------------
# 01 MetaPhlan
# -----------------------------------------------------
# merge viruses with predicted host
# virus abundance alongside precited host abundance
