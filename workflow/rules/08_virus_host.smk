# -----------------------------------------------------
# Virus Host Module (if include_host_module = True)
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
report: "report/workflow.rst"


# -----------------------------------------------------
# Virus host rules
# -----------------------------------------------------


# -----------------------------------------------------
# 00 Determine inputs for module
# -----------------------------------------------------
# if input_data = "viruses" symlink input viruses
rule symlink_preprocessed_viruses:
    message:
        "Symlinking preprocessed viruses to new paths"
    input:
        lambda wildcards: samples_df[(samples_df["sample"]) == wildcards.sample][
            "viruses"
        ].iloc[0],
    output:
        results + "00_INPUT/{sample}_preprocessed_viruses.fasta",
    benchmark:
        "benchmark/06_VIRUS_HOST/symlink_preprocessed_viruses_{sample}.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # symlink input viruses to renamed files
        ln -s {input} {output}
        """


if (
    config["input_data"] == "reads"
    or config["input_data"] == "contigs"
    or config["input_data"] == "vls"
):
    viruses = (
        results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps_untrimmed.fasta",
    )
elif config["input_data"] == "viruses":
    viruses = results + "00_INPUT/{sample}_viruses.fasta"


# -----------------------------------------------------
# 01 iPHoP
# -----------------------------------------------------
# download iphop database
rule download_iphop:
    message:
        "Downloading iPHoP database"
    output:
        resources + "iphop/iphop_download_complete",
    params:
        iphop_dir=resources + "iphop/",
    # conda:
    #     "../envs/iphop:1.1.0.yaml"
    container:
        "/gscratch/stf/carsonjm/apptainer/iphop-1.1.0.sif"
    benchmark:
        "benchmark/08_VIRUS_HOST/download_iphop.tsv"
    resources:
        runtime="4h",
        mem_mb="10GB",
    shell:
        """
        # download iphop test database
        iphop download --db_dir {params.iphop_dir} \
        --no_prompt \
        --split 

        # create file indicating download is complete
        touch {output}
        """


# verify that iphop download is correct
rule verify_iphop_db:
    message:
        "Verifying iPHoP download"
    input:
        resources + "iphop/iphop_download_complete",
    output:
        resources + "iphop/iphop_download_verified",
    params:
        iphop_dir=resources + "iphop/Sept_2021_pub/",
    # conda:
    #     "../envs/iphop:1.1.0.yaml"
    container:
        "/gscratch/stf/carsonjm/apptainer/iphop-1.1.0.sif"
    benchmark:
        "benchmark/08_VIRUS_HOST/verify_iphop_db.tsv"
    resources:
        runtime="4h",
        mem_mb="10GB",
    shell:
        """
        # verify iphop download
        iphop download --db_dir {params.iphop_dir} \
        --full_verify

        # create file to indicate that download is correct
        touch {output}
        """


# run iphop splinter to
checkpoint iphop_split:
    message:
        "Splitting virus fasta for iPHoP"
    input:
        viruses=viruses,
    output:
        directory(results + "08_VIRUS_HOST/01_iphop/input_fastas/"),
    params:
        out_dir=results + "08_VIRUS_HOST/01_iphop/input_fastas/",
        num_seqs=config["virus_host"]["iphop_split_num_seqs"],
    # conda:
    #     "../envs/iphop:1.1.0.yaml"
    container:
        "/gscratch/stf/carsonjm/apptainer/iphop-1.1.0.sif"
    benchmark:
        "benchmark/08_VIRUS_HOST/iphop_split.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        rm -rf {params.out_dir}
        mkdir {params.out_dir}

        # run iphop predict
        iphop split \
        --input_file {input.viruses} \
        --split_dir {params.out_dir} \
        --n_seq {params.num_seqs}
        """


# run iphop to predict hosts for viral sequences
rule iphop:
    message:
        "Predicting hosts for viral sequences"
    input:
        viruses=results + "08_VIRUS_HOST/01_iphop/input_fastas/batch_{batch}.fna",
        iphop=resources + "iphop/iphop_download_verified",
    output:
        results
        + "08_VIRUS_HOST/01_iphop/{batch}/Host_prediction_to_genus_m"
        + str(config["virus_host"]["iphop_min_score"])
        + ".csv",
    params:
        min_score=config["virus_host"]["iphop_min_score"],
        out_dir=results + "08_VIRUS_HOST/01_iphop/{batch}/",
        db_dir=resources + "iphop/Sept_2021_pub/",
        extra_args=config["virus_host"]["iphop_arguments"],
    # conda:
    #     "../envs/iphop:1.1.0.yaml"
    container:
        "/gscratch/stf/carsonjm/apptainer/iphop-1.1.0.sif"
    threads: config["virus_host"]["iphop_threads"]
    benchmark:
        "benchmark/08_VIRUS_HOST/iphop_{batch}.tsv"
    resources:
        runtime=config["virus_host"]["iphop_runtime"],
        mem_mb=config["virus_host"]["iphop_memory"],
    shell:
        """
        # run iphop predict
        iphop predict \
        --fa_file {input.viruses} \
        --out_dir {params.out_dir} \
        --db_dir {params.db_dir} \
        --num_threads {threads} \
        --min_score {params.min_score} \
        {params.extra_args}
        """


# input function for rule aggregate, return paths to all files produced by the checkpoint 'somestep'
def combine_iphop_input(wildcards):
    checkpoint_output = checkpoints.iphop_split.get(**wildcards).output[0]
    return expand(
        results
        + "08_VIRUS_HOST/01_iphop/{batch}/Host_prediction_to_genus_m"
        + str(config["virus_host"]["iphop_min_score"])
        + ".csv",
        batch=glob_wildcards(
            os.path.join(checkpoint_output, "batch_{batch}.fna")
        ).batch,
    )


# combine checkv reports across samples
rule combine_iphop_reports:
    message:
        "Combining iPHoP reports splits"
    input:
        combine_iphop_input,
    output:
        results + "08_VIRUS_HOST/01_iphop/iphop_report.csv",
    benchmark:
        "benchmark/08_VIRUS_HOST/combine_iphop_reports.csv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# -----------------------------------------------------
# Host analysis
# -----------------------------------------------------
# analyze virus host outputs
rule virus_host_analysis:
    message:
        "Visualizing iPHoP host taxonomy outputs"
    input:
        derep_reps=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_clusters.tsv",
        iphop=results + "08_VIRUS_HOST/01_iphop/iphop_report.csv",
    output:
        svg=report(
            results + "08_VIRUS_HOST/virus_host_taxonomy_figure.svg",
            category="Step 08: Virus hosts",
        ),
        report=results + "08_VIRUS_HOST/virus_host_taxonomy_report.tsv",
    benchmark:
        "benchmark/08_VIRUS_HOST/virus_host_analysis.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/08_virus_host_analysis.py"
