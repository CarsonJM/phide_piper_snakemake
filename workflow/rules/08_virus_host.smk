# -----------------------------------------------------
# Virus host
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
rule symlink_preprocessed_viruses:
    input:
        lambda wildcards: samples_df[(samples_df["sample"]) == wildcards.sample][
            "viruses"
        ].iloc[0],
    output:
        results + "00_INPUT/{sample}_preprocessed_viruses.fasta",
    benchmark:
        "benchmark/06_VIRUS_HOST/symlink_preprocessed_viruses_{sample}.tsv"
    resources:
        runtime="00:00:10",
        mem_mb="1000",
    shell:
        """
        # symlink input paths to renamed files
        ln -s {input} {output}
        """


if (
    config["input_data"] == "reads"
    or config["input_data"] == "contigs"
    or config["input_data"] == "vls"
):
    viruses = (
        results
        + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_untrimmed.fna"
    )
elif config["input_data"] == "viruses":
    viruses = results + "00_INPUT/{sample}_viruses.fasta"


# -----------------------------------------------------
# 01 iPHoP
# -----------------------------------------------------
rule download_iphop:
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
        runtime="04:00:00",
        mem_mb="50000",
    shell:
        """
        # download iphop test database
        iphop download --db_dir {params.iphop_dir} \
        --no_prompt \
        --split 

        touch {output}
        """


rule verify_iphop_db:
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
        "benchmark/08_VIRUS_HOST/download_iphop.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="50000",
    shell:
        """
        # download iphop test database
        iphop download --db_dir {params.iphop_dir} \
        --full_verify

        touch {output}
        """


rule iphop:
    input:
        viruses=viruses,
        iphop=resources + "iphop/iphop_download_verified",
    output:
        results
        + "08_VIRUS_HOST/01_iphop/Host_prediction_to_genus_m"
        + str(config["virus_host"]["iphop_min_score"])
        + ".csv",
    params:
        min_score=config["virus_host"]["iphop_min_score"],
        out_dir=results + "08_VIRUS_HOST/01_iphop/",
        db_dir=resources + "iphop/Sept_2021_pub/",
        extra_args=config["virus_host"]["iphop_arguments"],
    # conda:
    #     "../envs/iphop:1.1.0.yaml"
    container:
        "/gscratch/stf/carsonjm/apptainer/iphop-1.1.0.sif"
    threads: config["virus_host"]["iphop_threads"]
    benchmark:
        "benchmark/08_VIRUS_HOST/iphop.tsv"
    resources:
        runtime="4:00:00",
        mem_mb="100000",
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


# -----------------------------------------------------
# 02 PHIST
# -----------------------------------------------------
# build phist
rule build_phist:
    output:
        resources + "phist/out/predictions.csv",
    params:
        phist_dir=resources + "phist",
    conda:
        "../envs/phist.yml"
    benchmark:
        "benchmark/08_VIRUS_HOST/build_phist.tsv"
    resources:
        runtime="1:00:00",
        mem_mb="1000",
    shell:
        """
        # git clone phist
        rm -rf {params.phist_dir}
        git clone --recurse-submodules https://github.com/refresh-bio/PHIST {params.phist_dir}

        # build phist
        cd {params.phist_dir}
        make
        mkdir ./out

        # test phist
        python3 phist.py ./example/virus ./example/host ./out/common_kmers.csv ./out/predictions.csv
        """


# organize viruses to by input into phist
rule split_viruses_for_phist:
    input:
        viruses,
    output:
        results + "08_VIRUS_HOST/02_phist/virus_fastas/viruses_prepared",
    params:
        fasta_dir=results + "08_VIRUS_HOST/02_phist/virus_fastas/",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/08_VIRUS_HOST/split_viruses_for_phist.tsv"
    resources:
        runtime="1:00:00",
        mem_mb="1000",
    script:
        "../scripts/08_split_viruses_for_phist.py"


# run phist using uhgg
rule phist:
    input:
        phist_build=resources + "phist/out/predictions.csv",
        virus_fastas=results + "08_VIRUS_HOST/02_phist/virus_fastas/viruses_prepared",
    output:
        table=results + "08_VIRUS_HOST/02_phist/common_kmers.csv",
        predictions=results + "08_VIRUS_HOST/02_phist/predictions.csv",
    params:
        virus_dir=results + "08_VIRUS_HOST/02_phist/virus_fastas/",
        bacteria_db_dir=config["bacteria_db"],
        phist_script=resources + "phist/phist.py",
        out_dir=results + "08_VIRUS_HOST/02_phist/",
    conda:
        "../envs/phist.yml"
    threads: config["virus_host"]["phist_threads"]
    benchmark:
        "benchmark/08_VIRUS_HOST/phist.tsv"
    resources:
        runtime="10:00:00",
        mem_mb="100000",
    shell:
        """
        # run phist using uhgg
        python3 {params.phist_script} {params.virus_dir} {params.bacteria_db_dir} {params.out_dir} \
        -t {threads}
        """


# determine host taxonomy from refseq phist
rule phist_host_taxonomy:
    input:
        bacteria_db_metadata=config["bacteria_db_meta"],
        phist=results + "08_VIRUS_HOST/02_phist/common_kmers.csv",
    output:
        taxonomy=results + "08_VIRUS_HOST/02_phist/phist_host_taxonomy.csv",
        report=results + "08_VIRUS_HOST/02_phist/phist_host_report.csv",
    params:
        min_common_kmers=config["virus_host"]["min_phist_common_kmers"],
        min_agreement=config["virus_host"]["min_phist_agreement"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/08_VIRUS_HOST/phist_host_taxonomy.tsv"
    resources:
        runtime="1:00:00",
        mem_mb="1000",
    script:
        "../scripts/08_phist_host_taxonomy.py"


# -----------------------------------------------------
# Host analysis
# -----------------------------------------------------
rule virus_host_analysis:
    input:
        iphop=results
        + "08_VIRUS_HOST/01_iphop/Host_prediction_to_genus_m"
        + str(config["virus_host"]["iphop_min_score"])
        + ".csv",
        phist=results + "08_VIRUS_HOST/02_phist/phist_host_taxonomy.csv",
    output:
        results + "08_VIRUS_HOST/virus_host_taxonomy.csv",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/08_virus_host_analysis.py"
