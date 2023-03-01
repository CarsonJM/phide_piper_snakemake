# -----------------------------------------------------
# Virome enrichment Module (if input_data = "reads")
# -----------------------------------------------------
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


# -----------------------------------------------------
# Enrichment rules
# -----------------------------------------------------
localrules:
    combine_viromeqc_across_samples,
    virus_enrichment_analysis,


# -----------------------------------------------------
# 01 ViromeQC
# -----------------------------------------------------
# clone viromeqc and download database
rule build_viromeqc:
    message:
        "Cloning ViromeQC and downloading the database"
    output:
        amphora=resources + "viromeqc/index/amphora_bacteria.dmnd",
        lsu=resources + "viromeqc/index/SILVA_132_LSURef_tax_silva.clean.1.bt2",
        ssu=resources + "viromeqc/index/SILVA_132_SSURef_Nr99_tax_silva.clean.1.bt2",
    params:
        viromeqc_dir=resources + "viromeqc",
        index_dir=resources + "viromeqc/index",
        viromeqc_tmp=resources + "viromeqc/tmp/",
    conda:
        "../envs/viromeqc.yml"
    benchmark:
        "benchmark/02_VIRUS_ENRICHMENT/build_viromeqc.tsv"
    resources:
        runtime="4h",
        mem_mb="10GB",
    shell:
        """
        # clone viromeqc
        rm -rf {params.viromeqc_dir}
        git clone --recurse-submodules https://github.com/SegataLab/viromeqc.git {params.viromeqc_dir}

        cd {params.viromeqc_dir}
        python viromeQC.py --install

        # # install the viromeqc databases
        # mkdir {params.index_dir}
        # cd {params.index_dir}
        # wget -O amphora_markers.zip https://zenodo.org/record/4020594/files/amphora_markers.zip?download=1
        # wget -O SILVA_132_LSURef_tax_silva_clean.zip https://zenodo.org/record/4020594/files/SILVA_132_LSURef_tax_silva_clean.zip?download=1
        # wget -O SILVA_132_SSURef_Nr99_tax_silva.clean.zip https://zenodo.org/record/4020594/files/SILVA_132_SSURef_Nr99_tax_silva.clean.zip?download=1

        # # unzip databases
        # unzip amphora_markers.zip
        # unzip SILVA_132_LSURef_tax_silva_clean.zip
        # unzip SILVA_132_SSURef_Nr99_tax_silva.clean.zip

        # make tmp dir
        mkdir {params.viromeqc_tmp}
        """


# determine virus enrichment with viromeqc
rule viromeqc:
    message:
        "Running ViromeQC on {wildcards.group_sample} to determine virus enrichment"
    input:
        amphora=resources + "viromeqc/index/amphora_bacteria.dmnd",
        lsu=resources + "viromeqc/index/SILVA_132_LSURef_tax_silva.clean.1.bt2",
        ssu=resources + "viromeqc/index/SILVA_132_SSURef_Nr99_tax_silva.clean.1.bt2",
        R1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz",
        R2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz",
    output:
        results + "02_VIRUS_ENRICHMENT/01_viromeqc/{group_sample}_vqc.tsv",
    params:
        viromeqc_script=resources + "viromeqc/viromeQC.py",
        temp=resources + "viromeqc/tmp/{group_sample}",
        extra_arguments=config["virus_enrichment"]["viromeqc_arguments"],
    conda:
        "../envs/viromeqc.yml"
    benchmark:
        "benchmark/02_VIRUS_ENRICHMENT/viromeqc_{group_sample}.tsv"
    resources:
        runtime=config["virus_enrichment"]["viromeqc_runtime"],
        mem_mb=config["virus_enrichment"]["viromeqc_memory"],
    threads: config["virus_enrichment"]["viromeqc_threads"]
    shell:
        """
        # make dir to act as tmp
        mkdir {params.temp}

        # determine virome enrichment using viromeqc
        {params.viromeqc_script} \
        --input {input.R1} {input.R2} \
        --output {output} \
        --bowtie2_threads {threads} \
        --diamond_threads {threads} \
        --tempdir {params.temp} \
        {params.extra_arguments}

        # add sample column to each vqc output
        s={wildcards.group_sample}
        sed -i "s/$/\t$s/" {output}
        sample="sample"
        sed -i "1s/$s/$sample/" {output}

        # remove tmp dir
        rm -rf {params.temp}
        """


# combine viromeqc results across all samples
rule combine_viromeqc_across_samples:
    message:
        "Combining ViromeQC results across samples"
    input:
        expand(
            results + "02_VIRUS_ENRICHMENT/01_viromeqc/{group_sample}_vqc.tsv",
            group_sample=groups_samples,
        ),
    output:
        results + "02_VIRUS_ENRICHMENT/virus_enrichment_report.tsv",
    benchmark:
        "benchmark/02_VIRUS_ENRICHMENT/combine_viromeqc_across_samples.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    shell:
        """
        # combine all viromeqc outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# -----------------------------------------------------
# Analyze viral enrichment
# -----------------------------------------------------
# analyze viromeqc results to visualize virus enrichment
rule virus_enrichment_analysis:
    message:
        "Visualizing ViromeQC results for all samples"
    input:
        results + "02_VIRUS_ENRICHMENT/virus_enrichment_report.tsv",
    output:
        report(
            results + "02_VIRUS_ENRICHMENT/virus_enrichment_figure.svg",
            category="Step 02: Virus enrichment",
        ),
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/02_VIRUS_ENRICHMENT/virus_enrichment_analysis.tsv"
    resources:
        runtime="10m",
        mem_mb="10GB",
    script:
        "../scripts/02_virus_enrichment_analysis.py"
