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
        runtime="00:10:00",
        mem_mb="1000",
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
        + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_untrimmed.fna"
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
        runtime="4:00:00",
        mem_mb="10000",
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
        runtime="4:00:00",
        mem_mb="10000",
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
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
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
        runtime="12:00:00",
        mem_mb="100000",
        partition="compute-hugemem",
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
        "Combining iPHoP reports across samples"
    input:
        combine_iphop_input,
    output:
        results + "08_VIRUS_HOST/01_iphop/iphop_report.tsv",
    benchmark:
        "benchmark/08_VIRUS_HOST/combine_iphop_reports.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# -----------------------------------------------------
# 02 PHIST
# -----------------------------------------------------
# clone phist repository
rule build_phist:
    message:
        "Cloning PHIST repository"
    output:
        resources + "phist/out/predictions.csv",
    params:
        phist_dir=resources + "phist",
    conda:
        "../envs/phist.yml"
    benchmark:
        "benchmark/08_VIRUS_HOST/build_phist.tsv"
    resources:
        runtime="00:30:00",
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


# organize viruses for input into phist
rule split_viruses_for_phist:
    message:
        "Splitting viruses into single files for input into PHIST"
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
        runtime="00:30:00",
        mem_mb="1000",
    script:
        "../scripts/08_split_viruses_for_phist.py"


# run phist
rule phist:
    message:
        "Predicting hosts for viruses using PHIST"
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
        runtime="12:00:00",
        mem_mb="100000",
        partition="compute-hugemem",
    shell:
        """
        # run phist
        python3 {params.phist_script} {params.virus_dir} {params.bacteria_db_dir} {params.out_dir} \
        -t {threads}
        """


# determine host taxonomy from refseq phist
rule phist_host_taxonomy:
    message:
        "Predicting host taxonomy using all host genomes with > {params.min_common_kmers} common kmers and LCA where total kmers have > {params.min_agreement}% agreement"
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
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/08_phist_host_taxonomy.py"


# -----------------------------------------------------
# Host analysis
# -----------------------------------------------------
# analyze virus host outputs
rule virus_host_analysis:
    message:
        "Visualizing iPHoP and PHIST taxonomy outputs"
    input:
        iphop=results + "08_VIRUS_HOST/01_iphop/iphop_report.tsv",
        phist=results + "08_VIRUS_HOST/02_phist/phist_host_taxonomy.csv",
    output:
        report=results + "08_VIRUS_HOST/virus_host_taxonomy.csv",
        svg=report(
            results + "08_VIRUS_HOST/virus_host_taxonomy.svg",
            category="Step 08: Virus hosts",
        ),
        html=results + "08_VIRUS_HOST/virus_host_taxonomy.html",
    benchmark:
        "benchmark/08_VIRUS_HOST/virus_host_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/08_virus_host_analysis.py"
