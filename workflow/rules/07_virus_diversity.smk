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
    virus_diversity_genus_anlysis

# -----------------------------------------------------
# 01 Cluster viruses into vOTUs
# -----------------------------------------------------
# blast dereplicated viruses against one another
rule blast_derep_viruses_v_derep_viruses:
    message:
        "BLASTing all dereplicated viruses against one another for vOTU clustering"
    input:
        results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_viruses.fna",
    output:
        results + "07_VIRUS_DIVERSITY/01_votu_clusters/viruses_blast.tsv",
    params:
        blastdb=results + "07_VIRUS_DIVERSITY/01_votu_clusters/viruses_blastdb",
    # conda:
    #     "../envs/blast:2.12.0--h3289130_3.yml"
    container:
        "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/blast_viral_contigs_cluster.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    threads: config["virus_diversity"]["blast_threads"]
    shell:
        """
        # make a blast db from dereplicated viruses
        makeblastdb -in {input} -out {params.blastdb} -dbtype nucl

        # blast dereplicated viruses against one another
        blastn -query {input} -db {params.blastdb} -out {output} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
        """


# cluster viruses using blast results
rule cluster_viruses:
    message:
        "Clustering viruses into vOTUs based on {params.min_ani} ANI and {params.min_qcov} AF"
    input:
        fasta=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_viruses.fna",
        blast=results + "07_VIRUS_DIVERSITY/01_votu_clusters/viruses_blast.tsv",
    output:
        results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_clusters.tsv",
    params:
        blastani_script=resources + "mgv/ani_cluster/blastani.py",
        cluster_script=resources + "mgv/ani_cluster/cluster.py",
        ani_tsv=results + "07_VIRUS_DIVERSITY/01_votu_clusters/viruses_ani.tsv",
        min_ani=config["virus_diversity"]["min_ani"],
        min_qcov=config["virus_diversity"]["min_qcov"],
        min_tcov=config["virus_diversity"]["min_tcov"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/cluster_viruses.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # calculate ani and af from blast results
        python {params.blastani_script} -i {input.blast} -o {params.ani_tsv}

        # cluster phage genomes based on 95% ani and 85% af
        python {params.cluster_script} --fna {input.fasta} --ani {params.ani_tsv} --out {output} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov}
        """


# extract votu representatives
rule extract_votu_reps:
    message:
        "Extracting longest member of each vOTU as representative"
    input:
        clusters=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_clusters.tsv",
        viruses=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_viruses.fna",
        untrimmed_viruses=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_untrimmed_viruses.fna",
        proteins=results + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicated_virus_proteins.faa",
    output:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
        untrimmed_viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_untrimmed.fna",
        proteins=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_proteins.faa",
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/extract_votu_reps.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="10000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/06_extract_dereplicated_viruses.py"


# -----------------------------------------------------
# 02 vContact2
# -----------------------------------------------------
# create gene2genome file for input to vcontact2
rule gene2genome:
    message:
        "Creating a file that links genes to genomes for vContact2"
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_proteins.faa",
    output:
        results + "07_VIRUS_DIVERSITY/02_vcontact2/vcontact2_gene2genome.csv",
    container:
        "docker://sonnenburglab/vcontact2"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/gene2genome.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # create gene2genome file for vcontact2
        vcontact2_gene2genome \
        --proteins {input} \
        --output {output} \
        --source-type Prodigal-FAA
        """


# run vcontact2 to determine the number of genus-level clusters
rule vcontact2:
    message:
        "Run vContact2 to group viruses into approximately genus-level clusters"
    input:
        proteins=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_proteins.faa",
        g2g=results + "07_VIRUS_DIVERSITY/02_vcontact2/vcontact2_gene2genome.csv",
    output:
        clusters=results + "07_VIRUS_DIVERSITY/02_vcontact2/output/viral_cluster_overview.csv",
        network=results + "07_VIRUS_DIVERSITY/02_vcontact2/output/c1.ntw",
    params:
        out_dir=results + "07_VIRUS_DIVERSITY/02_vcontact2/output/"
    container:
        "docker://sonnenburglab/vcontact2"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/vcontact2.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="100000",
        partition="compute-hugemem",
    threads: config["virus_diversity"]["vcontact2_threads"]
    shell:
        """
        # clear output directory
        rm -rf {params.out_dir}

        # run vcontact2
        vcontact2 \
        --raw-proteins {input.proteins} \
        --rel-mode 'Diamond' \
        --proteins-fp {input.g2g} \
        --pcs-mode MCL \
        --vcs-mode ClusterONE \
        --c1-bin /opt/conda/bin/cluster_one-1.0.jar \
        --threads {threads} \
        --output-dir {params.out_dir}
        """


# -----------------------------------------------------
# 03 Caudoviricetes phylogeny (if include_phylogeny_module = True)
# -----------------------------------------------------
# download vogdb hmms
rule prepare_vogdb_hmms:
    message:
        "Downloading VOGdb HMMs"
    output:
        resources + "ccp77/vogdb_downloaded",
    params:
        ccp77=resources + "ccp77/vog.hmm.tar.gz",
        ccp77_dir=resources + "ccp77/",
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/prepare_vogdb_hmms.tsv"
    resources:
        runtime="00:30:00",
        mem_mb="10000",
    shell:
        """
        # download vogdb hmms
        wget http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz -O {params.ccp77}

        # extract vogdb hmms
        cd {params.ccp77_dir}
        tar -xvf {params.ccp77}
        touch {output}
        """


# extract ccp77 hmms from vogdb dataset
rule extract_ccp77_hmms:
    message:
        "Extracting 77 CCP HMMs specified in https://doi.org/10.1038/s41564-019-0448-z"
    input:
        resources + "ccp77/vogdb_downloaded",
        ccp77="workflow/resources/CCP77",
    output:
        resources + "ccp77/ccp77/VOG21972.hmm",
    params:
        vogdb_dir=resources + "ccp77/",
        ccp77_dir=resources + "ccp77/ccp77/",
        ccp77_hmms=resources + "ccp77/ccp77/ccp77.hmm",
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/extract_ccp77_hmms.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    conda:
        "../envs/jupyter.yml"
    script:
        "../scripts/07_filter_vogdb.py"


# hmmsearch ccp77 hmms
rule ccp77_hmmsearch:
    message:
        "Searching viral genes against CCP77 HMMs to identify Caudoviricetes marker genes"
    input:
        resources + "ccp77/ccp77/VOG21972.hmm",
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_proteins.faa",
    output:
        ccp77_hmms=resources + "ccp77/ccp77/ccp77.hmm",
        table=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/01_hmmsearch/ccp77.out",
    params:
        ccp77_dir=resources + "ccp77/ccp77/",
    # conda:
    #     "../envs/hmmer:3.3.2--h87f3376_2.yml"
    container:
        "docker://quay.io/biocontainers/hmmer:3.3.2--h87f3376_2"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/ccp77_hmmsearch.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="10000",
        partition="compute-hugemem",
    threads: config["virus_phylogeny"]["hmmsearch_threads"]
    shell:
        """
        # combine ccp77 hmms
        cat {params.ccp77_dir}*.hmm > {output.ccp77_hmms}

        # hmmsearch ccp77
        hmmsearch --noali \
        --cpu {threads} \
        --tblout {output.table} \
        {output.ccp77_hmms} \
        {input.viruses}
        """


# extract the top hits from the hmmsearch
rule extract_top_ccp77_hits:
    message:
        "Extracting top CCP77 HMM hits for each genome"
    input:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_proteins.faa",
        hmmsearch=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/01_hmmsearch/ccp77.out",
    output:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/top_hits_extracted",
    params:
        max_evalue=config["virus_phylogeny"]["max_evalue"],
        out_dir=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/02_top_hits/",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/extract_top_ccp77_hits.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/07_extract_top_hits.py"


# align top ccp77 hmm hits
rule ccp77_hmmalign:
    message:
        "Aligning all CCP77 hits for all viral genomes"
    input:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/top_hits_extracted",
    output:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/top_hits_aligned",
    params:
        in_dir=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/02_top_hits/",
        hmm_dir=resources + "ccp77/ccp77/",
        out_dir=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/03_hmmalign/",
    # conda:
    #     "../envs/hmmer:3.3.2--h87f3376_2.yml"
    container:
        "docker://quay.io/biocontainers/hmmer:3.3.2--h87f3376_2"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/ccp77_hmmalign.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="10000",
        partition="compute-hugemem",
    shell:
        """
        # change cwd
        cd {params.in_dir}

        # make output dir
        mkdir {params.out_dir}

        # align hits to CCP77 for each viral genome
        for file in *.faa
        do
        hmm={params.hmm_dir}"${{file%%.faa*}}".hmm
        hmmalign $hmm $file > {params.out_dir}$file.sth
        esl-reformat clustal {params.out_dir}$file.sth > {params.out_dir}$file.aln
        done

        touch {output}
        """


# trim multiple sequence alignment
rule trim_msa:
    message:
        "Trimming MSA using TrimAL"
    input:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/top_hits_aligned",
    output:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/msa_trimmed",
    params:
        in_dir=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/03_hmmalign/",
        out_dir=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/04_trimal/",
    # conda:
    #     "../envs/trimal:1.4.1--h9f5acd7_6.yml"
    container:
        "docker://quay.io/biocontainers/trimal:1.4.1--h9f5acd7_6"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/trim_msa.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="10000",
        partition="compute-hugemem",
    shell:
        """
        # change cwd
        cd {params.in_dir}
        mkdir {params.out_dir}

        # trim msa using trimal
        for file in *.aln
        do
        trimal -gt 0.5 -in $file -out {params.out_dir}$file
        done

        touch {output}
        """


# concatenate and filter msa
rule concatenate_and_filter_msa:
    message:
        "Concatenating MSA and filtering based on a minimum of {params.min_hits} hits per genome and max of {params.max_percent_gaps} percent gaps"
    input:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/msa_trimmed",
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_proteins.faa",
    output:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/05_concatenated_msa/concat.faa",
    params:
        max_percent_gaps=config["virus_phylogeny"]["max_percent_gaps"],
        min_hits=config["virus_phylogeny"]["min_hits_per_genome"],
        in_dir=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/04_trimal/",
        out_dir=results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/05_concatenated_msa/",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/concatenate_and_filter_msa.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/07_concatenate_and_filter.py"


# create phylogenetic tree using fasttree
rule fasttree:
    message:
        "Creating phylogenetic tree using fasttree"
    input:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/05_concatenated_msa/concat.faa",
    output:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/06_fasttree/fasttree",
    # conda:
    #     "../envs/fasttree:2.1.10--h779adbc_6.yml"
    container:
        "docker://quay.io/biocontainers/fasttree:2.1.10--h779adbc_6"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/fasttree.tsv"
    resources:
        runtime="12:00:00",
        mem_mb="50000",
        partition="compute-hugemem",
    shell:
        """
        fasttree -wag -gamma -mlacc 2 -slownni {input} > {output}
        """


# visualize virus phylogeny tree
rule virus_phylogeny_analysis:
    message:
        "Visualizing phylogenetic tree using biopython"
    input:
        results + "07_VIRUS_DIVERSITY/03_caudoviricetes_phylogeny/06_fasttree/fasttree",
    output:
        results + "07_VIRUS_DIVERSITY/virus_phylogeny.png",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/virus_phylogeny_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/07_virus_phylogeny_analysis.py"


# -----------------------------------------------------
# Analyze clustering
# -----------------------------------------------------
# analyze diversity at votu level
rule virus_diversity_votu_anlysis:
    message:
        "Visualizing vOTU diversity"
    input:
        results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_clusters.tsv",
    output:
        svg=report(
            results + "07_VIRUS_DIVERSITY/virus_diversity_votu_figure.svg",
            category="Step 07: Virus clustering",
        ),
        html=results + "07_VIRUS_DIVERSITY/virus_diversity_votu_figure.html",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/virus_diversity_votu_anlysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/06_virus_dereplication_analysis.py"


# analyze diversity at genus level via a network
rule virus_diversity_genus_anlysis:
    message:
        "Visualizing genus-level diversity and networks"
    input:
        network=results + "07_VIRUS_DIVERSITY/02_vcontact2/output/c1.ntw",
        clusters=results + "07_VIRUS_DIVERSITY/02_vcontact2/output/viral_cluster_overview.csv",
    output:
        network=report(
            results + "07_VIRUS_DIVERSITY/virus_diversity_genus_network.png",
            category="Step 07: Virus clustering",
        ),
        svg=report(
            results + "07_VIRUS_DIVERSITY/virus_diversity_genus_figure.svg",
            category="Step 07: Virus clustering",
        ),
        html=results + "07_VIRUS_DIVERSITY/virus_diversity_genus_figure.html",
    conda:
        "../envs/networkx.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/virus_diversity_genus_anlysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/07_virus_network_analysis.py"

