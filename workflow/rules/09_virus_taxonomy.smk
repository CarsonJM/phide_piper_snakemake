# -----------------------------------------------------
# Virus Taxonomy Module (if include_taxonomy_module = True)
# -----------------------------------------------------
import pandas as pd
import os


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples_df"], sep="\t")


# get current working directory so absolute paths can be used for input/output files
results = config["results"]


# load report
report: "report/workflow.rst"


# load resources folder path
resources = config["resources"]


# -----------------------------------------------------
# Virus taxonomy rules
# -----------------------------------------------------

# -----------------------------------------------------
# 01 geNomad's marker based taxonomy
# -----------------------------------------------------
# run genomad to identify virus taxonomy
rule genomad_taxonomy:
    message:
        "Running geNomad to predict virus taxonomy"
    input:
        genomad=resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
        contigs=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
    output:
        results
        + "09_VIRUS_TAXONOMY/01_genomad/votu_representatives_annotate/votu_representatives_taxonomy.tsv",
    params:
        out_dir=results + "09_VIRUS_TAXONOMY/01_genomad/",
        genomad_dir=resources + "genomad/genomad_db",
    conda:
        "../envs/genomad:1.0.1.yml"
    threads: config["virus_taxonomy"]["genomad_threads"]
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/genomad.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    shell:
        """
        # run genomad on viruses to predict taxonomy
        genomad end-to-end \
        --threads {threads} \
        --verbose \
        --min-score 0.0 \
        --cleanup \
        --splits {threads} \
        --enable-score-calibration \
        {input.contigs} \
        {params.out_dir} \
        {params.genomad_dir}
        """

# -----------------------------------------------------
# 02 MMseqs consensus
# -----------------------------------------------------
# download mmseqs database and extract 
rule build_mmseqs:
    message:
        "Downloading and extracting MMSeqs2 NCBI NR database for assigning viral taxonomy"
    output:
        resources + "imgvr_6/virus_tax_db/virus_tax_db"
    params:
        mmseqs_dir=resources + "resources/imgvr_6",
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/build_mmseqs.tsv"
    resources:
        runtime="00:30:00",
        mem_mb="1000",
    shell:
        """
        # download database
        cd {params.mmseqs_dir}
        wget https://zenodo.org/record/6574914/files/virus_tax_db.tar.zst?download=1
        mv 'virus_tax_db.tar.zst?download=1' virus_tax_db.tar.zst
        zstd -d virus_tax_db.tar.zst
        tar -xvf virus_tax_db.tar
        """

# run mmseqs2 against ncbi nr
rule mmseqs2:
    message:
        "Running MMSeqs2 against NCBI NR to predict consensus taxonomy"
    input:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
        db=resources + "imgvr_6/virus_tax_db/virus_tax_db"
    output:
        results + "09_VIRUS_TAXONOMY/02_mmseqs2/taxonomy_lca.tsv",
    params:
        out_dir=results + "09_VIRUS_TAXONOMY/02_mmseqs2/taxonomy",
        extra_args=config["virus_taxonomy"]["mmseqs_arguments"]
    # conda:
    #     "../envs/mmseqs2:14.7e284--pl5321hf1761c0_0.yml"
    container:
        "docker://quay.io/biocontainers/mmseqs2:14.7e284--pl5321hf1761c0_0"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/mmseqs2.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="150000",
    threads: config["virus_taxonomy"]["mmseqs_threads"]
    shell:
        """
        # create output directory
        mkdir {params.out_dir}
        cd {params.out_dir}

        # run mmseqs easy taxonomy to identify lca
        mmseqs easy-taxonomy \
        {input.viruses} \
        {input.db} \
        {params.out_dir} \
        tmp \
        --threads {threads} \
        -e 1e-5 \
        -s 6 \
        --blacklist "" \
        --tax-lineage 1 \
        {params.extra_args}
        """


# -----------------------------------------------------
# 03 External comparison
# -----------------------------------------------------
# compare assembled phages to external database
rule external_mash:
    message:
        "Running MASH against external database to determine novelty of assembled phages"
    input:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
        virus_db=config["virus_db"] + ".msh",
    output:
        sketch=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna.msh",
        mash=results + "09_VIRUS_TAXONOMY/03_external_comparison/votu_reps_v_external_mash.tsv",
    params:
        out_dir=results + "09_VIRUS_TAXONOMY/03_external_comparison/"
    # conda:
    #     "../envs/mash:2.3--ha61e061_0.yml"
    container:
        "docker://quay.io/biocontainers/mash:2.3--ha61e061_0"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/external_mash.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="150000",
    threads: config["virus_comparison"]["mash_threads"]
    shell:
        """
        # create output directory
        rm -rf {params.out_dir}
        mkdir {params.out_dir}

        # mash sketch assembled phage genomes
        mash sketch \
        {input.viruses} \
        -s 1000 \
        -i \
        -p {threads}

        # run mash distance against reference genomes
        mash dist \
        {output.sketch} {input.virus_db} \
        -p {threads} > {output.mash}
        """


# filter highly similar genomes for species comparisons
rule identify_novel_genera_and_families:
    message:
        "Identifying approximate genus/family level hits to reference databases"
    input:
        mash_dist=results + "09_VIRUS_TAXONOMY/03_external_comparison/votu_reps_v_external_mash.tsv",
        virus_db=config["virus_db"],
        meta=config["virus_db_meta"],
    output:
        mash_results=results + "09_VIRUS_TAXONOMY/03_external_comparison/genera_and_family_comparisons.csv",
        hi_sim_viruses=results + "09_VIRUS_TAXONOMY/03_external_comparison/hi_sim_viruses.fasta",
    params:
        out_dir=results + "09_VIRUS_TAXONOMY/03_external_comparison/",
        max_dist_species=config["virus_comparison"]["max_dist_species"],
        max_dist_genus=config["virus_comparison"]["max_dist_genus"],
        max_dist_family=config["virus_comparison"]["max_dist_family"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/identify_novel_genera_and_families.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/09_identify_novel_genera_and_families.py"


# blast viruses against highly similar database genomes and determine ANI and AF
rule blast_viruses_v_external:
    message:
        "BLASTing viruses against external database and calculating ANI and AF"
    input:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
        hi_sim_viruses=results + "09_VIRUS_TAXONOMY/03_external_comparison/hi_sim_viruses.fasta",
    output:
        results + "09_VIRUS_TAXONOMY/03_external_comparison/viruses_v_external_blast.tsv",
    params:
        blastdb=results + "09_VIRUS_TAXONOMY/03_external_comparison/hi_sim_viruses",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/blast_viruses_v_external.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["virus_comparison"]["blast_threads"]
    shell:
        """
        # make a blast db from dereplicated viruses
        makeblastdb -in {input.hi_sim_viruses} -out {params.blastdb} -dbtype nucl

        # blast dereplicated viruses against one another
        blastn -query {input.viruses} -db {params.blastdb} -out {output} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
        """


# blast viruses against highly similar database genomes and determine ANI and AF
rule ani_external:
    message:
        "BLASTing viruses against external database and calculating ANI and AF"
    input:
        results + "09_VIRUS_TAXONOMY/03_external_comparison/viruses_v_external_blast.tsv",
    output:
        results + "09_VIRUS_TAXONOMY/03_external_comparison/viruses_v_external_ani.tsv",
    params:
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/ani_external.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # calculate ani and af from blast results
        python {params.blastani_script} -i {input} -o {output}
        """

# identify top hit for each virus
rule identify_novel_species:
    message:
        "Identifying top external hits and determining if viruses represent novel species"
    input:
        results + "09_VIRUS_TAXONOMY/03_external_comparison/viruses_v_external_ani.tsv",
    output:
        results + "09_VIRUS_TAXONOMY/03_external_comparison/species_comparisons.tsv",
    params:
        min_ani=config["virus_comparison"]["min_ani"],
        min_qcov=config["virus_comparison"]["min_virus_cov"],
        min_tcov=config["virus_comparison"]["min_ref_cov"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/identify_novel_species.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/09_identify_novel_species.py"


# -----------------------------------------------------
# Analyze taxonomy results
# -----------------------------------------------------
# visualize virus taxonomy results
rule virus_taxonomy_analysis:
    message:
        "Visualizing virus taxonomy results from both tools"
    input:
        genomad=results + "09_VIRUS_TAXONOMY/01_genomad/votu_representatives_annotate/votu_representatives_taxonomy.tsv",
        mmseqs=results + "09_VIRUS_TAXONOMY/02_mmseqs2/taxonomy_lca.tsv",
    output:
        svg=report(
            results + "09_VIRUS_TAXONOMY/virus_taxonomy_analysis.svg",
            category="Step 09: Virus taxonomy",
        ),
        html=results + "09_VIRUS_TAXONOMY/virus_taxonomy_analysis.html",
        report=results + "09_VIRUS_TAXONOMY/virus_taxonomy_report.tsv",
    params:
        min_genomad_agreement=config["virus_taxonomy"]["genomad_min_agreement"],
        min_mmseqs_agreement=config["virus_taxonomy"]["mmseqs_min_agreement"],
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/virus_taxonomy_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/09_virus_taxonomy_analysis.py"