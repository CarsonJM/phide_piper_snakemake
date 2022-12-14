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
localrules: 
    combine_viruses,
    virus_dereplication_analysis,



# -----------------------------------------------------
# 01 Combine viruses across samples
# -----------------------------------------------------
# combine viruses across samples
rule combine_viruses:
    message:
        "Combining viruses, untrimmed viruses, and proteins across samples"
    input:
        expand(results + "05_VIRUS_QUALITY/02_quality_filter/{sample}/quality_filtered_viruses.fna", sample=samples),
    output:
        results
        + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/combine_viruses.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # combine viruses
        cat {input} > {output}
        """


# -----------------------------------------------------
# 02 Dereplicate viruses across samples
# -----------------------------------------------------
rule make_dereplicate_blastdb:
    input: 
        results + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
    output: 
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/combined_viruses_blastdb.ndb",
    params:
        db=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/combined_viruses_blastdb",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/make_dereplicate_blastdb.tsv"
    resources:
        runtime=config["virus_dereplication"]["blast_runtime"],
        mem_mb=config["virus_dereplication"]["blast_memory"],
    threads: config["virus_dereplication"]["blast_threads"]
    shell: 
        """
        makeblastdb -dbtype nucl -in {input} -out {params.db}
        """

rule dereplicate_blast:
    input:
        fasta=results + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        blastdb=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/combined_viruses_blastdb.ndb",
    output: 
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_blast.tsv",
    params:
        min_blast_ident = config["virus_dereplication"]["blast_min_id"],
        db=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/combined_viruses_blastdb",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/dereplicate_blast.tsv"
    resources:
        runtime=config["virus_dereplication"]["blast_runtime"],
        mem_mb=config["virus_dereplication"]["blast_memory"],
    threads: config["virus_dereplication"]["blast_threads"]
    shell:
        """
        blastn -task megablast \
        -max_target_seqs 25000 \
        -perc_identity {params.min_blast_ident} \
        -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue qlen slen" \
        -num_threads {threads} \
        -query {input.fasta} \
        -db {params.db} \
        -out {output}
        """

rule dereplicate_anicalc:
    input:
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_blast.tsv",
    output: 
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_ani.tsv",
    params:
        blast_max_evalue = config['virus_dereplication']['blast_max_evalue'],
    conda:
        "../envs/leiden_clustering.yml"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/dereplicate_anicalc.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/06_anicalc.py"

rule dereplicate_aniclust:
    input:
        fasta=results + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        ani=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_ani.tsv",
    output: 
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_clusters.tsv",
    params:
        leiden_resolution = config['virus_dereplication']['leiden_resolution'],
        min_ani = config['virus_dereplication']['min_ani'],
        min_cov = config['virus_dereplication']['min_cov'],
        avg_ani = config['virus_dereplication']['avg_ani'],
        seed = config['virus_dereplication']['random_seed'],
    conda:
        "../envs/leiden_clustering.yml"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/dereplicate_aniclust.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script: 
        "../scripts/06_aniclust.py"


rule get_replicate_representatives:
    input:
        fasta=results + "06_VIRUS_DEREPLICATION/01_combine_viruses/combined_viruses.fasta",
        clusters=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_clusters.tsv",
    output:
        representatives_list=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_cluster_reps.txt",
        representatives_fasta=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps.fasta",
    # conda:
    #     "../envs/seqkit:2.1.0--h9ee0642_0"
    container:
        "docker://quay.io/biocontainers/seqkit:2.1.0--h9ee0642_0"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/get_replicate_representatives.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        awk '{{print $1}}' {input.clusters} > {output.representatives_list} && \
        seqkit grep -f {output.representatives_list} {input.fasta} > {output.representatives_fasta}
        """


# -----------------------------------------------------
# Analyze clustering
# -----------------------------------------------------
# plot virus dereplication results
rule virus_dereplication_anlysis:
    message:
        "Visualizing the virus dereplication results"
    input:
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_clusters.tsv",
    output:
        report(
            results + "06_VIRUS_DEREPLICATION/virus_dereplication_figure.svg",
            category="Step 06: Virus dereplication",
        ),
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/06_VIRUS_DEREPLICATION/virus_dereplication_anlysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/06_virus_dereplication_analysis.py"
