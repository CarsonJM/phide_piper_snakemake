# -----------------------------------------------------
# Virus Dereplication
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
    input:
        viruses=expand(results + "05_VIRUS_QUALITY/03_quality_filter/{sample}/quality_filtered_viruses.fna", sample=samples),
        untrimmed_viruses=expand(results + "05_VIRUS_QUALITY/03_quality_filter/{sample}/untrimmed_quality_filtered_viruses.fna", sample=samples),
        proteins=expand(results + "05_VIRUS_QUALITY/03_quality_filter/{sample}/quality_filtered_proteins.faa", sample=samples)
    output:
        viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        untrimmed_viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_untrimmed_viruses.fasta",
        proteins=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_virus_proteins.faa",
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/combine_viruses.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        cat {input.viruses} > {output.viruses}
        cat {input.untrimmed_viruses} > {output.untrimmed_viruses}
        cat {input.proteins} > {output.proteins}
        """


# -----------------------------------------------------
# 02 Dereplicate viruses across samples
# -----------------------------------------------------
# blast viral genomes against one another
rule blast_viruses_v_viruses:
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
        runtime="12:00:00",
        mem_mb="06000",
        partition="compute-hugemem",
    threads: config["virus_dereplication"]["blast_threads"]
    shell:
        """
        # make a blast db from phage contigs
        makeblastdb -in {input} -out {params.blastdb} -dbtype nucl

        # all against all blast
        blastn -query {input} -db {params.blastdb} -out {output} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
        """


# dereplicate viral genomes
rule dereplicate_viruses:
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
        mem_mb="1000",
    shell:
        """
        # calculate ani and af from blast results
        python {params.blastani_script} -i {input.blast_tsv} -o {params.ani_tsv}

        # cluster phage genomes based on 95% ani and 85% af
        python {params.cluster_script} --fna {input.viruses} --ani {params.ani_tsv} --out {output} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov}
        """


# extract replicate representatives
rule extract_virus_replicate_representatives:
    input:
        clusters=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/viruses_replicates.tsv",
        viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        untrimmed_viruses=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_untrimmed_viruses.fasta",
        proteins=results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_virus_proteins.faa",
    output:
        viruses=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_viruses.fna",
        untrimmed_viruses=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_untrimmed_viruses.fna",
        proteins=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_virus_proteins.faa",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/extract_virus_replicate_representatives.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/06_extract_dereplicated_viruses.py"


# -----------------------------------------------------
# Analyze clustering
# -----------------------------------------------------
# plot virus dereplication results
rule virus_dereplication_anlysis:
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
        mem_mb="1000",
    script:
        "../scripts/06_virus_dereplication_analysis.py"
