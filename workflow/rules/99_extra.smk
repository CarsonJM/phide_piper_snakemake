# -----------------------------------------------------
# 02 CripsrIdentify
# -----------------------------------------------------
rule build_crispridentify:
    output:
        resources + "crispridentify/CRISPRidentify.py",
    params:
        ci_dir=resources + "crispridentify/",
    shell:
        """
        # clone crispridentify repo
        git clone https://github.com/BackofenLab/CRISPRidentify.git {params.ci_dir}
        """


rule crispridentify:
    input:
        ci_script=resources + "crispridentify/CRISPRidentify.py",
    output:
        results + "13_SPACER_PREDICTION/02_crispridentify/Complete_spacer_dataset.fasta",
    params:
        bacteria_db_dir=config["bacteria_db_dir"],
        out_dir=results + "13_SPACER_PREDICTION/02_crispridentify/",
        ci_dir=resources + "crispridentify/",
    conda:
        "../envs/crispridentify.yml"
    threads: config["spacer_prediction"]["crispridentify_threads"]
    shell:
        """
        cd {params.ci_dir}

        python {input.ci_script} --input_folder {params.bacteria_db_dir} \
        --result_folder {params.out_dir} --fasta_report True --parallel {threads}
        """


# -----------------------------------------------------
# 03 CripsrDetect
# -----------------------------------------------------
rule build_crisprdetect:
    output:
        resources + "crisprdetect/CRISPRDetect_2.4/CRISPRDetect.pl",
    params:
        cd_dir=resources + "crisprdetect/",
        cd=resources + "crisprdetect/CRISPRDetect_2.4.zip",
    shell:
        """
        rm -rf {params.cd_dir}
        mkdir {params.cd_dir}
        cd {params.cd_dir}
        wget https://github.com/davidchyou/CRISPRDetect_2.4/raw/master/CRISPRDetect_2.4.zip -O {params.cd}
        unzip {params.cd}
        """


rule crisprdetect:
    input:
        script=resources + "crisprdetect/CRISPRDetect_2.4/CRISPRDetect.pl",
        bacteria_db=config["bacteria_db"],
    output:
        results + "13_SPACER_PREDICTION/03_crisprdetect/merged_spacers.gff",
    params:
        bacteria_db_dir=config["bacteria_db"].rpartition("/")[0],
        out_dir=results + "13_SPACER_PREDICTION/03_crisprdetect/",
        cd_dir=resources + "crisprdetect/CRISPRDetect_2.4/",
    conda:
        "../envs/crisprdetect.yml"
    threads: config["spacer_prediction"]["crisprdetect_threads"]
    log:
        results + "00_LOGS/13_crisprdetect.log",
    shell:
        """
        # change to cd_dir
        cd {params.cd_dir}

        touch {output}

        for file in {params.bacteria_db_dir}/*.fasta
        do
        file_suffix=$(sed 's%{params.bacteria_db_dir}/%%g' <<< $file)
        perl {input.script} -f $file -array_quality_score_cutoff 3 -o {params.out_dir}$file_suffix -T {threads} >> {log}
        cat {params.out_dir}"$file_suffix".gff >> {output}
        rm {params.out_dir}$file_suffix*
        done
        """


# -----------------------------------------------------
# 06 Cenote Taker2
# -----------------------------------------------------
# build kraken database using custom virus database
rule build_cenote_taker2:
    output:
        resources + "cenote_taker2/run_cenote-taker2.py",
    params:
        ct2_dir=resources + "cenote_taker2/",
    conda:
        "../envs/cenote_taker2.yml"
    container:
        "/home/carsonjm/apptainer/cenote_taker2.sif"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/build_cenote_taker2.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # clone cenote taker2 repository
        git clone https://github.com/mtisza1/Cenote-Taker2.git {params.ct2_dir}
        """


# download cenote taker2 databases
rule download_cenote_taker2:
    input:
        resources + "cenote_taker2/run_cenote-taker2.py",
    output:
        resources + "cenote_taker2/taxdump/nodes.dmp",
    log:
        results + "00_LOGS/04_VIRUS_IDENTIFICATION/download_cenote_taker2.log",
    params:
        ct2_dir=resources + "cenote_taker2/",
    conda:
        "../envs/cenote_taker2.yml"
    container:
        "/home/carsonjm/apptainer/cenote_taker2.sif"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/download_cenote_taker2.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="1000",
    shell:
        """
        # change to cenote taker2 dir
        cd {params.ct2_dir}

        # download cenote_taker2 databases
        python update_ct2_databases.py \
        --hmm True \
        --protein True \
        --rps True \
        --taxdump True > {log} 2>&1
        """


# run cenote taker2 discovery
rule cenote_taker2:
    input:
        contigs=assembly,
        download=resources + "cenote_taker2/cenote_taker2_download_complete",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/06_cenote_taker2/{sample}/{sample}_CONTIG_SUMMARY.tsv",
    log:
        results + "00_LOGS/04_VIRUS_IDENTIFICATION/cenote_taker2_{sample}.log",
    params:
        script=resources + "cenote_taker2/unlimited_breadsticks.py",
        out_dir=results + "04_VIRUS_IDENTIFICATION/06_cenote_taker2/",
        min_circular_hallmarks=config["virus_identification"][
            "cenote_taker2_min_hallmark_circular"
        ],
        min_linear_hallmarks=config["virus_identification"][
            "cenote_taker2_min_hallmark_linear"
        ],
        extra_args=config["virus_identification"]["cenote_taker2_arguments"],
    container:
        "/home/carsonjm/apptainer/cenote_taker2.sif"
    threads: config["virus_identification"]["cenote_taker2_threads"]
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/cenote_taker2_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="32000",
    shell:
        """
        cd {params.out_dir}

        # run cenote taker2 discovery module
        python {params.script} \
        --contigs {input.contigs} \
        --run_title {wildcards.sample} \
        --mem 32 \
        --cpu {threads} \
        --lin_minimum_hallmark_genes {params.min_linear_hallmarks} \
        --circ_minimum_hallmark_genes {params.min_circular_hallmarks} \
        {params.extra_args} > {log} 2>&1
        """


# -----------------------------------------------------
# 05 vOTU consensus
# -----------------------------------------------------
rule parse_votu_clusters:
    input:
        viruses=results + "07_VIRUS_ABUNDANCE/03_filter_viruses/present_viruses.fna",
        clusters=results
        + "05_VIRUS_CLUSTERING/03_cluster_viruses/viruses_virusdb_votu_clusters.tsv",
        metadata=config["virus_db_meta"],
    output:
        taxonomy=results
        + "08_VIRUS_HOST/05_votu_consensus/votu_consensus_host_taxonomy.csv",
        report=results + "08_VIRUS_HOST/05_votu_consensus/votu_consensus_report.csv",
    params:
        min_agreement=config["virus_host"]["min_votu_agreement"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/08_votu_consensus.py.ipynb"


# -----------------------------------------------------
# 01 CRISPR Spacers
# -----------------------------------------------------
# make blastdb of uhgg spacers
rule make_spacer_blastdb:
    input:
        config["spacer_db"],
    output:
        resources + "crispr_spacers/crispr_spacers_db.ndb",
    params:
        spacers_blastdb=resources + "crispr_spacers/crispr_spacers_db",
    conda:
        "../envs/blast:2.12.0--h3289130_3.yml"
    shell:
        """
        # make blastdb
        makeblastdb \
        -in {input} \
        -out {params.spacers_blastdb} \
        -dbtype nucl
        """


# run crispropendb
rule blast_viruses_spacers:
    input:
        viruses=results + "07_VIRUS_ABUNDANCE/03_filter_viruses/present_viruses.fna",
        spacers_blastdb=resources + "crispr_spacers/crispr_spacers_db.ndb",
    output:
        results + "08_VIRUS_HOST/01_crispr_spacers/viruses_v_spacers_blast.tsv",
    params:
        spacers_blastdb=resources + "crispr_spacers/crispr_spacers_db",
    conda:
        "../envs/blast:2.12.0--h3289130_3.yml"
    threads: config["virus_host"]["blast_threads"]
    shell:
        """
        # blast viruses against uhgg spacers
        blastn \
        -query {input.viruses} \
        -db {params.spacers_blastdb} \
        -out {output} \
        -outfmt '6 std qlen slen' \
        -dust no \
        -word_size 18 \
        -num_threads {threads}
        """


# assign host taxonomy based on crispropendb blast
rule spacer_host_taxonomy:
    input:
        spacer_metadata=config["spacer_db_meta"],
        spacer_blast=results
        + "08_VIRUS_HOST/01_crispr_spacers/viruses_v_spacers_blast.tsv",
    output:
        report=results + "08_VIRUS_HOST/01_crispr_spacers/spacer_host_report.csv",
        taxonomy=results + "08_VIRUS_HOST/01_crispr_spacers/spacer_host_taxonomy.csv",
    params:
        max_mismatch=config["virus_host"]["max_spacer_mismatch"],
        min_spacer_coverage=config["virus_host"]["min_spacer_coverage"],
        min_agreement=config["virus_host"]["min_spacer_agreement"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/08_crispr_spacer_taxonomy.py.ipynb"


# -----------------------------------------------------
# 03 BLAST
# -----------------------------------------------------
# script to convert gff to fasta
# for file in /labdata3/hoffdata/Shared/resources/bacteria_db/bacteria/*.gff
# do
#     csplit $file /##FASTA/
#     rm xx00
#     mv xx01 $file.fasta
# done


rule extract_phist_hits:
    input:
        results + "08_VIRUS_HOST/03_phist/blast_genomes.csv",
    output:
        results + "08_VIRUS_HOST/02_blast/host_genomes_with_hits.fasta",
    shell:
        """
        cat $(cat {input}) > {output}
        """


rule make_phist_hits_blastdb:
    input:
        results + "08_VIRUS_HOST/02_blast/host_genomes_with_hits.fasta",
    output:
        results + "08_VIRUS_HOST/02_blast/host_genomes_with_hits.ndb",
    params:
        spacers_blastdb=results + "08_VIRUS_HOST/02_blast/host_genomes_with_hits",
    conda:
        "../envs/blast.yml"
    shell:
        """
        # make blastdb
        makeblastdb \
        -in {input} \
        -out {params.spacers_blastdb} \
        -dbtype nucl
        """


rule blast_viruses_hosts:
    input:
        viruses=results + "07_VIRUS_ABUNDANCE/03_filter_viruses/present_viruses.fna",
        spacers_blastdb=results + "08_VIRUS_HOST/02_blast/host_genomes_with_hits.ndb",
    output:
        results + "08_VIRUS_HOST/02_blast/host_genomes_blast_hits.tsv",
    params:
        spacers_blastdb=results + "08_VIRUS_HOST/02_blast/host_genomes_with_hits",
    conda:
        "../envs/blast.yml"
    threads: config["virus_host"]["blast_threads"]
    shell:
        """
        # blast viruses against uhgg spacers
        blastn \
        -query {input.viruses} \
        -db {params.spacers_blastdb} \
        -out {output} \
        -outfmt '6 std qlen slen' \
        -dust no \
        -word_size 18 \
        -num_threads {threads}
        """


rule blast_host_taxonomy:
    input:
        blast=results + "08_VIRUS_HOST/02_blast/host_genomes_blast_hits.tsv",
        bacteria_db_metadata=config["bacteria_db_meta"],
    output:
        taxonomy=results + "08_VIRUS_HOST/02_blast/blast_host_taxonomy.csv",
        report=results + "08_VIRUS_HOST/02_blast/blast_host_report.csv",
    params:
        min_length=config["virus_host"]["min_blast_len"],
        min_identity=config["virus_host"]["min_blast_identity"],
        min_agreement=config["virus_host"]["min_blast_agreement"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/08_blast_host_taxonomy.py.ipynb"


# -----------------------------------------------------
# 04 RAFAH
# -----------------------------------------------------
rule build_rafah:
    output:
        h3f=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3f",
        h3i=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3i",
        h3m=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3m",
        h3p=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3p",
        predict_script=resources + "rafah/RaFAH_Predict_Host.R",
        filename=resources + "rafah/MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData",
        rafah=resources + "rafah/RaFAH.pl",
        valid_cols=resources + "rafah/HP_Ranger_Model_3_Valid_Cols.txt",
    params:
        rafah_dir=resources + "rafah/",
    conda:
        "../envs/rafah.yml"
    shell:
        """
        # download rafah files
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/README.md
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_v0.3_hmm_models.tgz
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_v0.3_Ranger_Model.tgz
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH.pl
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/HP_Ranger_Model_3_Valid_Cols.txt
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_Predict_Host.R
        wget -P {params.rafah_dir} https://sourceforge.net/projects/rafah/files/RaFAH_v0.3_Files/RaFAH_Train_Model.R

        # decompress models
        cd {params.rafah_dir}
        tar -zxvf RaFAH_v0.3_hmm_models.tgz
        tar -zxvf RaFAH_v0.3_Ranger_Model.tgz
        """


# rule train_rafah:
#     input:
#     output:
#     log:
#     params:
#     conda:
#     threads:
#     benchmark:
#     resources:
#     shell:


rule rafah:
    input:
        h3f=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3f",
        h3i=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3i",
        h3m=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3m",
        h3p=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm.h3p",
        rafah=resources + "rafah/RaFAH.pl",
        valid_cols=resources + "rafah/HP_Ranger_Model_3_Valid_Cols.txt",
        predict_script=resources + "rafah/RaFAH_Predict_Host.R",
        filename=resources + "rafah/MMSeqs_Clusters_Ranger_Model_1+2+3_Clean.RData",
        viruses=results + "07_VIRUS_ABUNDANCE/03_filter_viruses/present_viruses.fna",
    output:
        results + "08_VIRUS_HOST/04_rafah/RaFAH_Seq_Info_Prediction.tsv",
    params:
        virus_dir=results + "07_VIRUS_ABUNDANCE/03_filter_viruses/",
        hmm=resources + "rafah/HP_Ranger_Model_3_Filtered_0.9_Valids.hmm",
        out_dir=results + "08_VIRUS_HOST/04_rafah",
        extra_args=config["virus_host"]["rafah_arguments"],
        min_score=config["virus_host"]["rafah_min_score"],
    conda:
        "../envs/rafah.yml"
    threads: config["virus_host"]["rafah_threads"]
    shell:
        """
        cd {params.out_dir}

        perl {input.rafah} --predict \
        --genomes_dir {params.virus_dir} \
        --extension .fna \
        --valid_ogs_file {input.valid_cols} \
        --hmmer_db_file_name {params.hmm} \
        --r_script_predict_file_name {input.predict_script} \
        --r_model_file_name {input.filename} \
        --threads {threads} \
        --min_cutoff {params.min_score} \
        {params.extra_args}
        """


rule rafah_taxonomy:
    input:
        rafah=results + "08_VIRUS_HOST/04_rafah/RaFAH_Seq_Info_Prediction.tsv",
    output:
        taxonomy=results + "08_VIRUS_HOST/04_rafah/rafah_host_taxonomy.csv",
        report=results + "08_VIRUS_HOST/04_rafah/rafah_consensus_report.csv",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/08_rafah_taxonomy.py.ipynb"

# -----------------------------------------------------
# 02 Cluster viruses into genera
# -----------------------------------------------------
rule vgenus_make_diamond_db:
    input:
        results + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_proteins.faa",
    output:
        results + "07_VIRUS_DIVERSITY/02_genera_clustering/virus_proteins.dmnd"
    log:
        results + "00_LOGS/07_VIRUS_DIVERSITY/vgenus_make_diamond_db.log"
    conda:
        "../envs/diamond:2.0.15--hb97b32f_1.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/vgenus_make_diamond_db.tsv"
    threads:
        config["virus_diversity"]["diamond_threads"]
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        diamond makedb \
        --in {input} \
        --db {output} \
        --threads {threads}
        """


rule vgenus_all_v_all_diamond:
    input:
        proteins=results + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_proteins.faa",
        db=results + "07_VIRUS_DIVERSITY/02_genera_clustering/virus_proteins.dmnd"
    output:
        results + "07_VIRUS_DIVERSITY/02_genera_clustering/virus_v_virus_blastp.tsv"
    log:
        results + "00_LOGS/07_VIRUS_DIVERSITY/vgenus_all_v_all_diamond.log"
    params:
        out_dir=results + "07_VIRUS_DIVERSITY/02_genera_clustering/vcontact2/",
    conda:
        "../envs/diamond:2.0.15--hb97b32f_1.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/vgenus_all_v_all_diamond.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="50000",
    threads:
        config["virus_diversity"]["diamond_threads"]
    shell:
        """
        diamond blastp \
        --query {input.proteins} \
        --db {input.db} \
        --out {output} \
        --outfmt 6 \
        --evalue 1e-5 \
        --max-target-seqs 10000 \
        --query-cover 50 \
        --subject-cover 50
        """


rule vgenus_compute_and_filter_aai:
    input:
        proteins=results + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_proteins.faa",
        amino_acid_script=resources + "mgv/aai_cluster/amino_acid_identity.py",
        blastp_tsv=results + "07_VIRUS_DIVERSITY/02_genera_clustering/virus_v_virus_blastp.tsv"
    output:
        results + "07_VIRUS_DIVERSITY/02_genera_clustering/genus_edges.tsv",
    log:
        results + "00_LOGS/07_VIRUS_DIVERSITY/vgenus_compute_and_filter_aai.log"
    params:
        aai_tsv=results + "07_VIRUS_DIVERSITY/02_genera_clustering/aai.tsv",
        filter_aai_script=resources + "mgv/aai_cluster/filter_aai.py",
        min_percent_shared=config["virus_diversity"]["min_percent_shared"],
        min_num_shared=config["virus_diversity"]["min_num_shared"],
        min_aai=config["virus_diversity"]["min_aai"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/vgenus_compute_and_filter_aai.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="50000",
    shell:
        """
        python {input.amino_acid_script} \
        --in_faa {input.proteins} \
        --in_blast {input.blastp_tsv} \
        --out_tsv {params.aai_tsv}

        python {params.filter_aai_script} \
        --in_aai {params.aai_tsv} \
        --min_percent_shared {params.min_percent_shared} \
        --min_num_shared {params.min_num_shared} \
        --min_aai {params.min_aai} \
        --out_tsv {output}
        """


rule vgenus_mcl_clustering:
    input:
        results + "07_VIRUS_DIVERSITY/02_genera_clustering/genus_edges.tsv",
    output:
        results + "07_VIRUS_DIVERSITY/02_genera_clustering/genus_clusters.tsv",
    log:
        results + "00_LOGS/07_VIRUS_DIVERSITY/vgenus_mcl_clustering.log"
    conda:
        "../envs/mcl:14.137--pl5321hec16e2b_8.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/vgenus_mcl_clustering.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="50000",
    shell:
        """
        mcl {input} \
        -te 8 \
        -I 2.0 \
        --abc \
        -o {output}
        """


rule virus_diversity_vgenus_anlysis:
    input:
        clusters=results + "07_VIRUS_DIVERSITY/02_genera_clustering/genus_clusters.tsv",
    output:
        report(
            results + "07_VIRUS_DIVERSITY/virus_diversity_vgenus_figure.png",
            caption="../report/05_virus_clustering.rst",
            category="Step 05: Virus clustering",
        ),
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/07_VIRUS_DIVERSITY/virus_diversity_vgenus_anlysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/07_virus_diversity_vgenus_analysis.py"



# -----------------------------------------------------
# 03 Gene consensus
# -----------------------------------------------------
# download ictv taxonomy release
rule build_ictv:
    output:
        fix_taxdump_script=resources + "ictv/fix_taxdump.py",
        get_taxids_script=resources + "ictv/get_ictv_taxids.py"
    params:
        ictv_dir=resources + "ictv/",
        ictv_xlsx="ictv.xlsx",
        ictv_tsv=resources + "ictv/ictv.tsv",
    conda:
        "../envs/ictv.yml"
    shell:
        """
        rm -rf {params.ictv_dir}
        git clone https://github.com/apcamargo/ictv-mmseqs2-protein-database.git {params.ictv_dir}
        cd {params.ictv_dir}
        """


rule download_and_prepare_ictv_data:
    input:
        fix_taxdump_script=resources + "ictv/fix_taxdump.py",
        get_taxids_script=resources + "ictv/get_ictv_taxids.py",
    output:
        nr_faa=resources + "ictv/nr.virus.faa.gz",
        acc2taxid=resources + "ictv/nr.virus.accession2taxid.ictv"
    params:
        ictv_dir=resources + "ictv/",
    conda:
        "../envs/ictv.yml"
    script:
        "../scripts/ictv_taxa_download.sh"


rule prepare_ictv_mmseqs_db:
    input:
        nr_faa=resources + "ictv/nr.virus.faa.gz",
        acc2taxid=resources + "ictv/nr.virus.accession2taxid.ictv"
    output:
        resources + "ictv/virus_tax_db/virus_tax_db.mmsdb",
    params:
        ictv_dir=resources + "ictv/",
        virus_tax_db_dir=resources + "ictv/virus_tax_db/",
        virus_tax_db=resources + "ictv/virus_tax_db/virus_tax_db",
        ictv_taxdump=resources + "ictv/ictv-taxdump",
    conda:
        "../envs/ictv.yml"
    shell:
        """
        cd {params.ictv_dir}

        # Create the MMSeqs2 database
        mkdir {params.virus_tax_db_dir}
        mmseqs createdb --dbtype 1 {input.nr_virus_faa} {params.virus_tax_db}
        mmseqs createtaxdb {params.virus_tax_db} tmp --ncbi-tax-dump {params.ictv_taxdump}--tax-mapping-file {input.nr_virus_acc2taxid}
        rm -rf tmp
        """


rule mmseqs2:
    input:
        viruses=results + "06_VIRUS_QUALITY/02_quality_filter/quality_filtered_proteins.faa",
        mmsdb=resources + "ictv/virus_tax_db/virus_tax_db.mmsdb",
    output:
        results + "09_VIRUS_TAXONOMY/04_mmseqs_ictv/test",
    params:
        virus_tax_db=resources + "ictv/virus_tax_db/virus_tax_db",
        out_dir=results + "09_VIRUS_TAXONOMY/04_mmseqs2_ictv/",
    conda:
        "../envs/ictv.yml"
    shell:
        """
        mmseqs easy-taxonomy {input.viruses} \
        {params.virus_tax_db} \
        {params.out_dir} tmp \
        -e 1e-5 \
        -s 6 \
        --blacklist "" \
        --tax-lineage 1
        """


# -----------------------------------------------------
# 01 vOTU with RefSeq
# -----------------------------------------------------
# rule blast_viruses_w_refseq:
#     input:
#         viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
#         refseq=resources + "viral_refseq/refseq_genomes.fna",
#     output:
#         combined=results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/viruses_w_refseq.fna",
#         blast=results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/viruses_blast.tsv",
#     params:
#         blastdb=results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/viruses_w_refseq_blastdb",
#     # conda:
#     #     "../envs/blast:2.12.0--h3289130_3.yml"
#     container:
#         "docker://quay.io/biocontainers/blast:2.12.0--h3289130_3"
#     benchmark:
#         "benchmark/09_VIRUS_TAXONOMY/blast_viruses_w_refseq.tsv"
#     resources:
#         runtime="04:00:00",
#         mem_mb="10000",
#     threads: config["virus_diversity"]["blast_threads"]
#     shell:
#         """
#         # combine viruses with refseq
#         cat {input.viruses} {input.refseq} > {output.combined}

#         # make a blast db from phage contigs
#         makeblastdb -in {output.combined} -out {params.blastdb} -dbtype nucl

#         # all against all blast
#         blastn -query {output.combined} -db {params.blastdb} -out {output.blast} -num_threads {threads} -outfmt '6 std qlen slen' -max_target_seqs 25000 -perc_identity 90
#         """


# rule cluster_viruses_w_refseq:
#     input:
#         fasta=results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/viruses_w_refseq.fna",
#         blast=results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/viruses_blast.tsv",
#     output:
#         results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/refseq_clusters.tsv",
#     params:
#         blastani_script=resources + "mgv/ani_cluster/blastani.py",
#         cluster_script=resources + "mgv/ani_cluster/cluster.py",
#         ani_tsv=results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/viruses_ani.tsv",
#         min_ani=config["virus_diversity"]["min_ani"],
#         min_qcov=config["virus_diversity"]["min_qcov"],
#         min_tcov=config["virus_diversity"]["min_tcov"],
#     conda:
#         "../envs/jupyter.yml"
#     benchmark:
#         "benchmark/09_VIRUS_TAXONOMY/cluster_viruses_w_refseq.tsv"
#     resources:
#         runtime="00:10:00",
#         mem_mb="1000",
#     shell:
#         """
#         # calculate ani and af from blast results
#         python {params.blastani_script} -i {input.blast} -o {params.ani_tsv}

#         # cluster phage genomes based on 99% ani and 99% af
#         python {params.cluster_script} --fna {input.fasta} --ani {params.ani_tsv} --out {output} --min_ani {params.min_ani} --min_qcov {params.min_qcov} --min_tcov {params.min_tcov}
#         """


# rule refseq_votu_taxonomy:
#     input:
#         viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_clusters.tsv",
#         refseq=results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/refseq_clusters.tsv",
#     output:
#         results + "09_VIRUS_TAXONOMY/01_votu_taxonomy/virus_taxonomy.tsv",
#     benchmark:
#         "benchmark/09_VIRUS_TAXONOMY/refseq_votu_taxonomy.tsv"
#     resources:
#         runtime="00:10:00",
#         mem_mb="1000",
#     conda:
#         "../envs/jupyter.yml"
#     script:
#         "../scripts/09_refseq_votu_taxonomy.py"


rule download_refseq:
    output:
        resources + "viral_refseq/refseq_genomes.fna",
        resources + "viral_refseq/refseq_proteins.faa",
    params:
        refseq_dir=resources + "viral_refseq/",
        refseq_gen_gz=resources + "viral_refseq/refseq_genomes.fna.gz",
        refseq_prot_gz=resources + "viral_refseq/refseq_proteins.faa.gz",
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/combine_viruses.tsv"
    resources:
        runtime="00:06:00",
        mem_mb="0600",
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz -P {params.refseq_dir}
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.1.genomic.fna.gz -P {params.refseq_dir}
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.1.genomic.fna.gz -P {params.refseq_dir}
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.1.genomic.fna.gz -P {params.refseq_dir}
        cat {params.refseq_dir}*.fna.gz > {params.refseq_gen_gz}
        gunzip {params.refseq_gen_gz}

        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz -P {params.refseq_dir}
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.protein.faa.gz -P {params.refseq_dir}
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.3.protein.faa.gz -P {params.refseq_dir}
        wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.4.protein.faa.gz -P {params.refseq_dir}
        cat {params.refseq_dir}*.faa.gz > {params.refseq_prot_gz}
        gunzip {params.refseq_prot_gz}
        """

rule build_refseq_diamond_db:
    input:
        resources + "viral_refseq/refseq_proteins.faa",
    output:
        results + "09_VIRUS_TAXONOMY/02_diamond_refseq/viral_refseq.dmnd"
    # conda:
    #     "../envs/diamond:2.0.15--hb97b32f_1.yml"
    container:
        "docker://quay.io/biocontainers/diamond:2.0.15--hb97b32f_1"
    threads: config["virus_taxonomy"]["diamond_threads"]
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/votu_consensus_taxonomy.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        diamond makedb \
        --in {input} \
        --db {output} \
        --threads {threads}
        """


rule viral_refseq_diamond:
    input:
        viruses=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives_proteins.faa",
        refseq_db=results + "09_VIRUS_TAXONOMY/02_diamond_refseq/viral_refseq.dmnd"
    output:
        results + "09_VIRUS_TAXONOMY/02_diamond_refseq/viral_refseq_diamond.tsv"
    # conda:
    #     "../envs/diamond:2.0.15--hb97b32f_1.yml"
    container:
        "docker://quay.io/biocontainers/diamond:2.0.15--hb97b32f_1"
    benchmark:
        "benchmark/09_VIRUS_TAXONOMY/votu_consensus_taxonomy.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        diamond blastp \
        --query {input.proteins} \
        --db {input.db} \
        --out {output} \
        --outfmt 6 \
        --evalue 1e-5 \
        --max-target-seqs 10000 \
        --query-cover 50 \
        --subject-cover 50 \
        --unal 1
        """


# rule viral_refseq_diamond_taxonomy:
#     input:
#         results + "09_VIRUS_TAXONOMY/03_diamond_refseq/viral_refseq_diamond.tsv"
#     output:
#         results + "09_VIRUS_TAXONOMY/03_diamond_refseq/viral_refseq_diamond_taxonomy.tsv"
#     # conda:
#     #     "../envs/diamond:2.0.15--hb97b32f_1.yml"
#     container:
#         "docker://quay.io/biocontainers/diamond:2.0.15--hb97b32f_1"
#     benchmark:
#         "benchmark/09_VIRUS_TAXONOMY/votu_consensus_taxonomy.tsv"
#     resources:
#         runtime="00:10:00",
#         mem_mb="10000",
#     shell:
#         """
#         """


# -----------------------------------------------------
# 03 VIBRANT
# -----------------------------------------------------
# add sample column to genomad output
rule add_sample_to_vibrant_lifestyle:
    input:
        results + "04_VIRUS_IDENTIFICATION/04_vibrant/{sample}/VIBRANT_{sample}_contigs_dereplicated/VIBRANT_results_{sample}_contigs_dereplicated/VIBRANT_integrated_prophage_coordinates_{sample}_contigs_dereplicated.tsv",
    output:
        results + "10_VIRUS_LIFESTYLE/03_vibrant/{sample}_vibrant_lifestyle.tsv",
    shell:
        """
        # add sample column to each checkv output
        s={wildcards.sample}
        sed -i "s/$/\t$s/" {input}
        sample="sample"
        sed -i "1s/$s/$sample/" {input}

        cat {input} > {output}
        """


# parse checkV to determine prophage integration
rule combine_vibrant_lifestyles:
    input:
        expand(results + "10_VIRUS_LIFESTYLE/03_vibrant/{sample}_vibrant_lifestyle.tsv", sample=samples)
    output:
        results + "10_VIRUS_LIFESTYLE/03_vibrant/combined_vibrant_lifestyle.tsv",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# -----------------------------------------------------
# 04 CheckV
# -----------------------------------------------------
