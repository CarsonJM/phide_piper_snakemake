# -----------------------------------------------------
# Virus function Module (if include_function_module = True)
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples_df"], sep="\t")
samples = samples_df['sample']


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -----------------------------------------------------
# Virus function rules
# -----------------------------------------------------


# -----------------------------------------------------
# 01 DRAMv
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
    viruses = results + "00_INPUT/{sample}_viruses.fasta"


# build dram
rule build_dramv:
    message:
        "Downloading DRAM databases"
    output:
        resources + "dram/DRAM_downloaded",
    params:
        dram_dir=resources + "dram/DRAM_data",
    conda:
        "../envs/dram:1.3.4--pyhdfd78af_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/dram:1.3.4--pyhdfd78af_0"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/build_dramv.tsv"
    resources:
        runtime=config["virus_function"]["dramv_runtime"],
        mem_mb=config["virus_function"]["dramv_memory"],
    threads: config["virus_function"]["dramv_threads"]
    shell:
        """
        # download dram databases except uniref
        DRAM-setup.py prepare_databases \
        --skip_uniref \
        --output_dir {params.dram_dir} \
        --threads {threads} \
        --verbose

        touch {output}
        """


# updated dram's config file with database locations
rule update_dram:
    message:
        "Updating DRAM's config file with database locations"
    input:
        resources + "dram/DRAM_downloaded",
    output:
        resources + "dram/DRAM_updated",
    params:
        dram_dir=resources + "dram/DRAM_data",
    conda:
        "../envs/dram:1.3.4--pyhdfd78af_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/dram:1.3.4--pyhdfd78af_0"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/update_dram.tsv"
    resources:
        runtime=config["virus_function"]["dramv_runtime"],
        mem_mb=config["virus_function"]["dramv_memory"],
    threads: config["virus_function"]["dramv_threads"]
    shell:
        """
        # update descriptions db
        DRAM-setup.py update_description_db

        touch {output}
        """


rule download_virsorter2:
    message:
        "Downloading VirSorter2 database"
    output:
        resources + "virsorter2/Done_all_setup",
    params:
        vs2_dir=resources + "virsorter2/",
    conda:
        "../envs/virsorter:2.2.3--pyhdfd78af_1.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/download_virsorter2.tsv"
    resources:
        runtime="30m",
        mem_mb="10GB",
    threads: config["virus_function"]["virsorter2_threads"]
    shell:
        """
        # download virsorter2 database
        # remove the whole directory specified by -d
        rm -rf db
        # run setup
        virsorter setup -d {params.vs2_dir} -j {threads}
        """

# run virsorter2 to identifiy viral contigs
rule virsorter2_dram:
    message:
        "Obtaining VirSorter2 output that is required to run DRAM-v"
    input:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_representatives.fasta",
        vs2_db=resources + "virsorter2/Done_all_setup",
    output:
        tab=results
        + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/viral-affi-contigs-for-dramv.tab",
        fasta=results + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/final-viral-combined-for-dramv.fa",
    params:
        vs2_db=resources + "virsorter2",
        out_dir=results + "11_VIRUS_FUNCTION/01_virsorter2/",
        extra_args=config["virus_function"]["virsorter2_arguments"],
    conda:
        "../envs/virsorter:2.2.3--pyhdfd78af_1.yml"
    threads: config["virus_function"]["virsorter2_threads"]
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/virsorter2_dram.tsv"
    resources:
        runtime=config["virus_function"]["virsorter2_runtime"],
        mem_mb=config["virus_function"]["virsorter2_memory"],
    threads:
        config["virus_function"]["virsorter2_threads"],
    shell:
        """
        # remove output directory
        rm -rf {params.out_dir}

        # run virsorter2
        virsorter run all --keep-original-seq \
        --db-dir {params.vs2_db} \
        -w {params.out_dir} \
        -i {input.viruses} \
        -j {threads} \
        --prep-for-dramv \
        --min-score 0.0 \
        --rm-tmpdir \
        --viral-gene-enrich-off \
        {params.extra_args}
        """


# run dramv to annotate proteins
rule dramv_annotate:
    message:
        "Annotating proteins using DRAM-v"
    input:
        dram=resources + "dram/DRAM_updated",
        vs2=results + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/viral-affi-contigs-for-dramv.tab",
        viruses=results + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/final-viral-combined-for-dramv.fa",
    output:
        results + "11_VIRUS_FUNCTION/02_dramv/annotations.tsv"
    params:
        out_dir = results + "11_VIRUS_FUNCTION/02_dramv/",
    conda:
        "../envs/dram:1.3.4--pyhdfd78af_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/dram:1.3.4--pyhdfd78af_0"
    resources:
        runtime=config["virus_function"]["dramv_runtime"],
        mem_mb=config["virus_function"]["dramv_memory"],
    threads: config["virus_function"]["dramv_threads"]
    shell:
        """
        rm -rf {params.out_dir}

        # annotate proteins with dramv
        DRAM-v.py annotate \
        -i {input.viruses} \
        -v {input.vs2} \
        --verbose \
        --keep_tmp_dir \
        -o {params.out_dir} \
        --threads {threads}
        """


# run dramv distill
rule dramv_distill:
    message:
        "Distilling DRAM-v annotations"
    input:
        results + "11_VIRUS_FUNCTION/02_dramv/annotations.tsv"
    output:
        report(
            results + "11_VIRUS_FUNCTION/02_dramv/distillate/product.html",
            category="Step 11: Virus Function",
        ),
    params:
        out_dir = results + "11_VIRUS_FUNCTION/02_dramv/distillate",
    conda:
        "../envs/dram:1.3.4--pyhdfd78af_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/dram:1.3.4--pyhdfd78af_0"
    resources:
        runtime="10m",
        mem_mb="10GB",
    threads: config["virus_function"]["dramv_threads"]
    shell:
        """
        rm -rf {params.out_dir}
        
        # run dramv distill
        DRAM-v.py distill \
        -i {input} \
        -o {params.out_dir}
        """