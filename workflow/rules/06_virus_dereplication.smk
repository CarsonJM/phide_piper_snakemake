# -----------------------------------------------------
# Virus Dereplication Module (if input_data = "reads" or "contigs" or "vls")
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
# Virus clustering rules
# -----------------------------------------------------
localrules: combine_viruses


# -----------------------------------------------------
# 01 Combine viruses across samples
# -----------------------------------------------------
# combine viruses across samples
rule combine_viruses:
    message:
        "Combining viruses, untrimmed viruses, and proteins across samples"
    input:
        viruses=expand(results + "05_VIRUS_QUALITY/03_quality_filter/{sample}/quality_filtered_viruses.fna", sample=samples),
        untrimmed_viruses=expand(results + "05_VIRUS_QUALITY/03_quality_filter/{sample}/untrimmed_quality_filtered_viruses.fna", sample=samples),
    output:
        viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        untrimmed_viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_untrimmed_viruses.fasta",
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/combine_viruses.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # combine viruses
        cat {input.viruses} > {output.viruses}

        # combine untrimmed viruses (for host identification)
        cat {input.untrimmed_viruses} > {output.untrimmed_viruses}
        """


# -----------------------------------------------------
# 02 Dereplicate viruses across samples
# -----------------------------------------------------
# blast viral genomes against one another
rule blast_viruses_v_viruses:
    message:
        "BLASTing viral genomes against one another for dereplication across samples"
    input:
        results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
    output:
        results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_v_viruses_blast.tsv",
    params:
        blastdb=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_v_viruses_blastdb",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/blast_viruses_v_viruses.tsv"
    resources:
        runtime=config["virus_dereplication"]["blast_runtime"],
        mem_mb=config["virus_dereplication"]["blast_memory"],
    threads: config["virus_dereplication"]["blast_threads"]
    shell:
        """
        # make a blast db from viruses
        makeblastdb -in {input} -out {params.blastdb} -dbtype nucl

        # all against all blast
        blastn -query {input} -db {params.blastdb} -out {output} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
        """


# dereplicate viruses
rule dereplicate_viruses:
    message:
        "Dereplicating viruses across samples at {params.min_ani} ANI and {params.min_qcov} AF"
    input:
        viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        blast_tsv=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_v_viruses_blast.tsv",
    output:
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_replicates.tsv",
    params:
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
        cluster_script=resources + "mgv/ani_cluster/cluster.py",
        ani_tsv=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_v_viruses_ani.tsv",
        min_ani=config["virus_dereplication"]["min_ani"],
        min_qcov=config["virus_dereplication"]["min_qcov"],
        min_tcov=config["virus_dereplication"]["min_tcov"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/dereplicate_viruses.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # calculate ani and af from blast results
        python {params.blastani_script} -i {input.blast_tsv} -o {params.ani_tsv}

        # cluster viruses 95% ani and 85% af (recommended)
        python {params.cluster_script} --fna {input.viruses} --ani {params.ani_tsv} --out {output} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov}
        """


# extract replicate representatives
rule extract_virus_replicate_representatives:
    message:
        "Extracting longest member of each replicate cluster as the representative"
    input:
        clusters=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_replicates.tsv",
        viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        untrimmed_viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_untrimmed_viruses.fasta",
    output:
        viruses=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_viruses.fna",
        untrimmed_viruses=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_untrimmed_viruses.fna",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/extract_virus_replicate_representatives.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/06_extract_dereplicated_viruses.py"


# -----------------------------------------------------
# Analyze clustering
# -----------------------------------------------------
# plot virus dereplication results
rule virus_dereplication_anlysis:
    message:
        "Visualizing the virus dereplication results"
    input:
        results
            + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_replicates.tsv",
    output:
        svg=report(
            results + "06_VIRUS_DEREPLICATION/virus_dereplication_figure.svg",
            category="Step 06: Virus dereplication",
        ),
        html=results + "06_VIRUS_DEREPLICATION/virus_dereplication_figure.html",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/virus_dereplication_anlysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/06_virus_dereplication_analysis.py"
