# -----------------------------------------------------
# Virus Diversity Module (if input_data = "reads" or "contigs" or "vls")
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
    virus_diversity_votu_anlysis

# -----------------------------------------------------
# 01 Species-level vOTU clustering
# -----------------------------------------------------
rule make_votu_blastdb:
    input: 
        results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps.fasta",
    output: 
        results + "07_VIRUS_DIVERSITY/01_votu_clustering/dereplicate_reps_blastdb.ndb",
    params:
        db=results + "07_VIRUS_DIVERSITY/01_votu_clustering/dereplicate_reps_blastdb",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/make_votu_blastdb.tsv"
    resources:
        runtime=config["virus_diversity"]["blast_runtime"],
        mem_mb=config["virus_diversity"]["blast_memory"],
    threads: config["virus_diversity"]["blast_threads"]
    shell: 
        """
        makeblastdb -dbtype nucl -in {input} -out {params.db}
        """

rule votu_blast:
    input:
        fasta=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps.fasta",
        blastdb=results + "07_VIRUS_DIVERSITY/01_votu_clustering/dereplicate_reps_blastdb.ndb",
    output: 
        results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_blast.tsv",
    params:
        min_blast_ident = config["virus_diversity"]["blast_min_id"],
        db=results + "07_VIRUS_DIVERSITY/01_votu_clustering/dereplicate_reps_blastdb",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/votu_blast.tsv"
    resources:
        runtime=config["virus_diversity"]["blast_runtime"],
        mem_mb=config["virus_diversity"]["blast_memory"],
    threads: config["virus_diversity"]["blast_threads"]
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

rule votu_anicalc:
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_blast.tsv",
    output: 
        results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_ani.tsv",
    params:
        blast_max_evalue = config['virus_diversity']['blast_max_evalue'],
    conda:
        "../envs/leiden_clustering.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/votu_anicalc.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/06_anicalc.py"

rule votu_aniclust:
    input:
        fasta=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps.fasta",
        ani=results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_ani.tsv",
    output: 
        results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_clusters.tsv",
    params:
        leiden_resolution = config['virus_diversity']['leiden_resolution'],
        min_ani = config['virus_diversity']['min_ani'],
        min_cov = config['virus_diversity']['min_cov'],
        avg_ani = config['virus_diversity']['avg_ani'],
        seed = config['virus_diversity']['random_seed'],
    conda:
        "../envs/leiden_clustering.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/votu_aniclust.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script: 
        "../scripts/06_aniclust.py"


rule get_votu_representatives:
    input:
        fasta=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps.fasta",
        clusters=results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_clusters.tsv",
    output:
        representatives_list=results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_cluster_reps.txt",
        representatives_fasta=results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_representatives.fasta",
    # conda:
    #     "../envs/seqkit:2.1.0--h9ee0642_0"
    container:
        "docker://quay.io/biocontainers/seqkit:2.1.0--h9ee0642_0"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/get_votu_representatives.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        awk '{{print $1}}' {input.clusters} > {output.representatives_list} && \
        seqkit grep -f {output.representatives_list} {input.fasta} > {output.representatives_fasta}
        """


# prodigal-gv
rule prodigal_gv:
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_representatives.fasta",
    output:
        faa=results + "07_VIRUS_DIVERSITY/02_proteins/votu_representatives_proteins.faa",
        fna=results + "07_VIRUS_DIVERSITY/02_proteins/votu_representatives_proteins.fna",
    params:
        extra_args=config['virus_diversity']['prodigal_gv_arguments']
    conda:
        "../envs/prodigal_gv.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/prodigal_gv.tsv"
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

# -----------------------------------------------------
# Analyze clustering
# -----------------------------------------------------
# plot virus dereplication results
rule virus_diversity_votu_anlysis:
    message:
        "Visualizing the virus votu clustering results"
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_clusters.tsv",
    output:
        svg=report(
            results + "07_VIRUS_DIVERSITY/virus_diversity_votu_figure.svg",
            category="Step 07: Virus diversity",
        ),
        report=results + "07_VIRUS_DIVERSITY/votu_diversity_report.tsv",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/votu_anlysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/06_virus_dereplication_analysis.py"