# -----------------------------------------------------
# Virus Identification Module (If input_data = "reads" or "contigs")
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
# Virus identification rules
# -----------------------------------------------------
localrules:
    symlink_contigs,
    symlink_preprocessed_reads,
    combine_reports_across_samples,
    merge_reports_within_samples,
    merge_viral_contigs_within_samples,


# -----------------------------------------------------
# 00 Determine input for module
# -----------------------------------------------------
# symlink contigs if contigs are input
rule symlink_contigs:
    input:
        lambda wildcards: samples_df[(samples_df["sample"]) == wildcards.sample][
            "contigs"
        ].iloc[0],
    output:
        results + "00_INPUT/{sample}_contigs.fasta",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/symlink_contigs_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # symlink input contigs to renamed files
        ln -s {input} {output}
        """


# filter symlinked contigs based on contig length
rule filter_symlinked_contigs:
    input:
        results + "00_INPUT/{sample}_contigs.fasta",
    output:
        results + "00_INPUT/{sample}_contigs_dereplicated.fasta",
    params:
        min_length=config["read_assembly"]["min_contig_length"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/filter_symlinked_contigs_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/03_contig_length_filter.py"


# if reads are input type, use output from 03_read_assembly
if config["input_data"] == "reads":
    assembly = (
        results
        + "03_READ_ASSEMBLY/02_contig_filters/{sample}/{sample}_contigs_dereplicated.fasta",
    )
# if contigs are input, use symlinked contigs
elif config["input_data"] == "contigs":
    assembly = (results + "00_INPUT/{sample}_contigs_dereplicated.fasta",)


# symlink preprocessed reads for external hits/abundances
rule symlink_preprocessed_reads:
    input:
        R1=lambda wildcards: samples_df[(+samples_df["sample"]) == wildcards.sample][
            "R1"
        ].iloc[0],
        R2=lambda wildcards: samples_df[(+samples_df["sample"]) == wildcards.sample][
            "R2"
        ].iloc[0],
    output:
        R1=results + "00_INPUT/{sample}_proprcessed_1.fastq.gz",
        R2=results + "00_INPUT/{sample}_preprocessed_2.fastq.gz",
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # symlink input reads to renamed files
        ln -s {input.R1} {output.R1}
        ln -s {input.R2} {output.R2}
        """


# if reads are input, use preprocessed reads from 01_read_preprocessing for external hits/abundances
if config["input_data"] == "reads":
    R1 = results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq.gz"
    R2 = results + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq.gz"
# if contigs, vls, or viruses are input then symlink reads for external hits/abundances
else:
    R1 = results + "00_INPUT/{sample}_proprocessed_1.fastq.gz"
    R2 = results + "00_INPUT/{sample}_preprocessed_2.fastq.gz"


# -----------------------------------------------------
# 01 MGV (& VirFinder)
# -----------------------------------------------------
# download mgv repo and HMM files
rule download_mgv_databases:
    message:
        "Downloading MGV HMM databases"
    output:
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
    params:
        mgv_dir=resources + "mgv",
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm.gz",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm.gz",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/download_mgv_databases.tsv"
    resources:
        runtime="00:30:00",
        mem_mb="1000",
    shell:
        """
        # download MGV HMM databases
        wget -O {params.imgvr_hmm} https://img.jgi.doe.gov//docs/final_list.hmms.gz
        wget -O {params.pfam_hmm} ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz

        # unzip MGV HMM databases
        gunzip {params.imgvr_hmm}
        gunzip {params.pfam_hmm}
        """


# use prodigal to identify ORFs
rule mgv_prodigal:
    message:
        "Predicting {wildcards.sample} ORFs using Prodigal"
    input:
        assembly,
    output:
        mgv_fna=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.fna",
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.faa",
        mgv_ffn=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.ffn",
    # conda:
    #     "../envs/prodigal:2.6.3--h779adbc_3.yml"
    container:
        "docker://quay.io/biocontainers/prodigal:2.6.3--h779adbc_3"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_prodigal_{sample}.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    shell:
        """
        # predict ope reading frames
        prodigal -i {input} \
        -p meta \
        -a {output.mgv_faa} \
        -d {output.mgv_ffn}

        # create a symlink to original contigs fasta files
        ln -s {input} {output.mgv_fna}
        """


# use hmmsearch to find imgvr HMM hits
rule mgv_imgvr_hmmsearch:
    message:
        "Searching IMGVR HMMs for {wildcards.sample} hits"
    input:
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.faa",
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_imgvr.out",
    # conda:
    #     "../envs/hmmer:3.1b2--2.yml"
    container:
        "docker://quay.io/biocontainers/hmmer:3.1b2--2"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_imgvr_hmmsearch_{sample}.tsv"
    resources:
        runtime="15:00:00",
        mem_mb="10000",
    threads: config["virus_identification"]["mgv_threads"]
    shell:
        """
        # identify imgvr hits using hmm search
        hmmsearch \
        -Z 1 \
        --cpu {threads} \
        --noali \
        --tblout {output} \
        {input.imgvr_hmm} {input.mgv_faa}
        """


# use hmmsearch to find pfam HMM hits
rule mgv_pfam_hmmsearch:
    message:
        "Searching PFAM HMMs for {wildcards.sample} hits"
    input:
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.faa",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_pfam.out",
    # conda:
    #     "../envs/hmmer:3.1b2--2.yml"
    container:
        "docker://quay.io/biocontainers/hmmer:3.1b2--2"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_pfam_hmmsearch_{sample}.tsv"
    resources:
        runtime="10:00:00",
        mem_mb="10000",
    threads: config["virus_identification"]["mgv_threads"]
    shell:
        """
        # identify pfam hits using hmmsearch
        hmmsearch \
        -Z 1 \
        --cpu {threads} \
        --noali \
        --tblout {output} \
        {input.pfam_hmm} {input.mgv_faa}
        """


# use mgv count_hmm_hits.py script to determine bacterial/viral HMM hits
rule mgv_count_hmm_hits:
    message:
        "Counting HMM hits for {wildcards.sample}"
    input:
        contigs=assembly,
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.faa",
        mgv_imgvr=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_imgvr.out",
        mgv_pfam=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_pfam.out",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_hmm_hits.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    # conda:
    #     "../envs/biopython:1.78.yml"
    container:
        "docker://quay.io/biocontainers/biopython:1.78"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_count_hmm_hits_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # change to mgv directory so scripts work correctly
        cd {params.mgv_dir}

        # count hmm hits using mgv script
        python count_hmm_hits.py \
        {input.contigs} {input.mgv_faa} \
        {input.mgv_imgvr} {input.mgv_pfam} > {output}
        """


# run virfinder to identify viral kmers
rule mgv_virfinder:
    message:
        "Running VirFinder on {wildcards.sample} contigs"
    input:
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
        contigs=assembly,
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_virfinder.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    # conda:
    #     "../envs/r-virfinder:1.1--r36he1b5a44_0.yml"
    container:
        "docker://quay.io/biocontainers/r-virfinder:1.1--r36he1b5a44_0"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_virfinder_{sample}.tsv"
    resources:
        runtime="05:00:00",
        mem_mb="10000",
    shell:
        """
        # change to mgv directory so scripts work correctly
        cd {params.mgv_dir}

        # run mgv virfinder script
        Rscript virfinder.R \
        {input.contigs} {output}

        # convert NA to none
        sed -i 's/\tNA\tNA/\t0.0\t1.0/g' {output}
        """


# calculate the strand switch rate using mgv strand_switch.py script
rule mgv_strand_switch:
    message:
        "Calculating strand switch rate for {wildcards.sample} contigs"
    input:
        imgvr_hmm=resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
        pfam_hmm=resources + "mgv/viral_detection_pipeline/input/pfam.hmm",
        contigs=assembly,
        mgv_faa=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.faa",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_strand_switch.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    # conda:
    #     "../envs/biopython:1.78.yml"
    container:
        "docker://quay.io/biocontainers/biopython:1.78"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_strand_switch_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # change to mgv directory so scripts work correctly
        cd {params.mgv_dir}

        # run mgv strand switch script
        python strand_switch.py \
        {input.contigs} {input.mgv_faa} > {output}
        """


# create master table using mgv master_table.py script
rule mgv_master_table:
    message:
        "Creating MGV master table for {wildcards.sample} analysis"
    input:
        mgv_hmm_hits=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_hmm_hits.tsv",
        mgv_vf=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_virfinder.tsv",
        mgv_strand_switch=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_strand_switch.tsv",
    output:
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_master_table.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_master_table_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # change to mgv directory so scripts run correctly
        cd {params.mgv_dir}

        # run mgv master table script
        python master_table.py \
        {input.mgv_hmm_hits} {input.mgv_vf} {input.mgv_strand_switch} > {output}
        """


# predict viral contigs using HMM hits, strand switch rate, and virfinder results
rule mgv_viral_classify:
    message:
        "Classifying {wildcards.sample} to identify viral contigs"
    input:
        mgv_fna=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}.fna",
        mgv_master_table=results
        + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_master_table.tsv",
    output:
        fna=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_final.fna",
        tsv=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_final.tsv",
    params:
        mgv_dir=resources + "mgv/viral_detection_pipeline",
        input_base=results + "04_VIRUS_IDENTIFICATION/01_mgv/input/{sample}",
        output_base=results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_final",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mgv_viral_classify_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # change to mgv directory so scripts run correctly
        cd {params.mgv_dir}

        # run mgv viral classify script
        python viral_classify.py \
        --features {input.mgv_master_table} \
        --in_base {params.input_base} \
        --out_base {params.output_base}
        """


# -----------------------------------------------------
# 02 VirSorter
# -----------------------------------------------------
# download virsorter db
rule download_virsorter:
    message:
        "Downloading VirSorter database"
    output:
        resources + "virsorter/virsorter-data-v2.tar.gz",
    params:
        output_dir=resources + "virsorter/",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/download_virsorter.tsv"
    # conda:
    #     "../envs/virsorter:1.0.6--pl526h516909a_0.yml"
    container:
        "docker://quay.io/biocontainers/virsorter:1.0.6--pl526h516909a_0"
    resources:
        runtime="00:30:00",
        mem_mb="1000",
    shell:
        """
        # download virsorter database
        wget https://zenodo.org/record/1168727/files/virsorter-data-v2.tar.gz -P {params.output_dir}
        """


# extract virsorter db
rule extract_virsorter_db:
    message:
        "Extracting VirSorter database"
    input:
        resources + "virsorter/virsorter-data-v2.tar.gz",
    output:
        vsrm=resources + "virsorter/virsorter-data/VirSorter_Readme.txt",
    params:
        vs_dir=resources + "virsorter/",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/extract_virsorter_db.tsv"
    resources:
        runtime="00:30:00",
        mem_mb="1000",
    shell:
        """
        # download virsorter database
        tar -xvzf {input} -C {params.vs_dir}
        """


# identify viral contigs using virsorter
rule virsorter:
    message:
        "Running VirSorter on {wildcards.sample} to identify viral contigs"
    input:
        contigs=assembly,
        vsrm=resources + "virsorter/virsorter-data/VirSorter_Readme.txt",
    output:
        vs_translation=results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{sample}/fasta/input_sequences_id_translation.tsv",
        vs_results=results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{sample}/Metric_files/VIRSorter_phage_signal.tab",
    params:
        output_dir=results + "04_VIRUS_IDENTIFICATION/02_virsorter/{sample}",
        vs_dir=resources + "virsorter/virsorter-data/",
        extra_args=config["virus_identification"]["virsorter_arguments"],
    # conda:
    #     "../envs/virsorter:1.0.6--pl526h516909a_0.yml"
    container:
        "docker://quay.io/biocontainers/virsorter:1.0.6--pl526h516909a_0"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/virsorter_{sample}.tsv"
    resources:
        runtime="25:00:00",
        mem_mb="10000",
    threads: config["virus_identification"]["virsorter_threads"]
    shell:
        """
        # clear output directory if present
        rm -rf {params.output_dir}

        # run virsorter to identify viral contigs
        wrapper_phage_contigs_sorter_iPlant.pl \
        --fna {input.contigs} \
        --wdir {params.output_dir} \
        --ncpu {threads} \
        --data-dir {params.vs_dir} \
        {params.extra_args}
        """


# -----------------------------------------------------
# 03 VirSorter2
# -----------------------------------------------------
# download virsorter2 db
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
        runtime="00:30:00",
        mem_mb="1000",
    threads: config["virus_identification"]["virsorter2_threads"]
    shell:
        """
        # download virsorter2 database
        # remove the whole directory specified by -d
        rm -rf db

        # run setup
        virsorter setup -d {params.vs2_dir} -j {threads}
        """


# run virsorter2 to identifiy viral contigs
rule virsorter2:
    message:
        "Running VirSorter2 for {wildcards.sample} to identify viral contigs"
    input:
        contigs=assembly,
        vs2_db=resources + "virsorter2/Done_all_setup",
    output:
        vls=results
        + "04_VIRUS_IDENTIFICATION/03_virsorter2/{sample}/final-viral-score.tsv",
        dramv=results
        + "04_VIRUS_IDENTIFICATION/03_virsorter2/{sample}/for-dramv/viral-affi-contigs-for-dramv.tab",
    params:
        vs2_db=resources + "virsorter2",
        vs2_dir=results + "04_VIRUS_IDENTIFICATION/03_virsorter2/{sample}",
        extra_args=config["virus_identification"]["virsorter2_arguments"],
    conda:
        "../envs/virsorter:2.2.3--pyhdfd78af_1.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/virsorter2_{sample}.tsv"
    resources:
        runtime="20:00:00",
        mem_mb="100000",
    threads: config["virus_identification"]["virsorter2_threads"]
    shell:
        """
        # clear output directory
        rm -rf {params.vs2_dir}

        # run virsorter2
        virsorter run all --keep-original-seq \
        --db-dir {params.vs2_db} \
        -w {params.vs2_dir} \
        -i {input.contigs} \
        -j {threads} \
        --prep-for-dramv \
        --min-score 0.0 \
        --rm-tmpdir \
        {params.extra_args}
        """


# -----------------------------------------------------
# 04 VIBRANT
# -----------------------------------------------------
# download vibrant database
rule download_vibrant:
    message:
        "Downloading VIBRANT database"
    output:
        resources + "vibrant/db/databases/VIBRANT_setup.log",
    params:
        vb_dir=resources + "vibrant",
    conda:
        "../envs/vibrant:1.2.1--hdfd78af_2.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/download_vibrant.tsv"
    resources:
        runtime="02:00:00",
        mem_mb="1000",
    shell:
        """
        # clear output directory
        rm -rf {params.vb_dir}

        # download vibrant database
        cd $CONDA_PREFIX/share/vibrant-1.2.1/db/databases
        ./VIBRANT_setup.py

        # create vibrant directory in resources folder
        mkdir {params.vb_dir}

        # move vibrant databases to resources folder
        cp -r $CONDA_PREFIX/share/vibrant-1.2.1/db {params.vb_dir}/db
        """


# run vibrant to identify viral contigs
rule vibrant:
    message:
        "Running VIBRANT on {wildcards.sample} to identify viral contigs"
    input:
        contigs=assembly,
        vb_db=resources + "vibrant/db/databases/VIBRANT_setup.log",
    output:
        results=results
        + "04_VIRUS_IDENTIFICATION/04_vibrant/{sample}/VIBRANT_{sample}_contigs_dereplicated/VIBRANT_phages_{sample}_contigs_dereplicated/{sample}_contigs_dereplicated.phages_combined.txt",
        vb_prophage=results
        + "04_VIRUS_IDENTIFICATION/04_vibrant/{sample}/VIBRANT_{sample}_contigs_dereplicated/VIBRANT_phages_{sample}_contigs_dereplicated/{sample}_contigs_dereplicated.phages_lysogenic.fna",
    params:
        vb_db=resources + "vibrant/db/databases/",
        vb_files=resources + "vibrant/db/files/",
        vb_dir=results + "04_VIRUS_IDENTIFICATION/04_vibrant/{sample}",
        extra_args=config["virus_identification"]["vibrant_arguments"],
    conda:
        "../envs/vibrant:1.2.1--hdfd78af_2.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/vibrant_{sample}.tsv"
    resources:
        runtime="10:00:00",
        mem_mb="10000",
    threads: config["virus_identification"]["vibrant_threads"]
    shell:
        """
        # remove output directory
        # rm -rf {params.vb_dir}

        # run vibrant
        VIBRANT_run.py \
        -i {input.contigs} \
        -d {params.vb_db} \
        -m {params.vb_files} \
        -folder {params.vb_dir} \
        -t {threads} \
        {params.extra_args}
        """


# -----------------------------------------------------
# 05 DeepVirFinder
# -----------------------------------------------------
# clone deepvirfinder repo
rule build_deepvirfinder:
    message:
        "Building DeepVirFinder"
    output:
        resources + "deepvirfinder/dvf.py",
    params:
        dvf_dir=resources + "deepvirfinder/",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/build_deepvirfinder.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # run deepvirfinder
        git clone https://github.com/jessieren/DeepVirFinder {params.dvf_dir}

        # make script executable
        chmod 777 {params.dvf_dir}*
        """


# run deepvirfinder
rule deepvirfinder:
    message:
        "Running DeepVirFinder on {wildcards.sample} to identify viral contigs"
    input:
        dvf_script=resources + "deepvirfinder/dvf.py",
        contigs=assembly,
    output:
        results
        + "04_VIRUS_IDENTIFICATION/05_deepvirfinder/{sample}_contigs_dereplicated.fasta_gt1000bp_dvfpred.txt",
    params:
        output_dir=results + "04_VIRUS_IDENTIFICATION/05_deepvirfinder/",
        model_dir=resources + "deepvirfinder/models",
    conda:
        "../envs/deepvirfinder.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/deepvirfinder_{sample}.tsv"
    resources:
        runtime="10:00:00",
        mem_mb="10000",
    threads: config["virus_identification"]["virfinder_threads"]
    shell:
        """
        # run deepvirfinder
        {input.dvf_script} \
        -i {input.contigs} \
        -o {params.output_dir} \
        -l 1000 \
        -c {threads}
        """


# -----------------------------------------------------
# 06 Genomad
# -----------------------------------------------------
# download genomad database
rule download_genomad:
    message:
        "Downloading geNomad database"
    output:
        resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
    params:
        genomad_dir=resources + "genomad/",
    conda:
        "../envs/genomad:1.0.1.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/download_genomad.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="10000",
    shell:
        """
        # change to genomad directory
        cd {params.genomad_dir}

        # download genomad databases
        genomad download-database .
        """


# run genomad to identify viral contigs
rule genomad:
    message:
        "Running geNomad on {wildcards.sample} to identify viral contigs"
    input:
        genomad=resources + "genomad/genomad_db/virus_hallmark_annotation.txt",
        contigs=assembly,
    output:
        viruses=results
        + "04_VIRUS_IDENTIFICATION/06_genomad/{sample}/{sample}_contigs_dereplicated_summary/{sample}_contigs_dereplicated_virus_summary.tsv",
        provirus=results
        + "04_VIRUS_IDENTIFICATION/06_genomad/{sample}/{sample}_contigs_dereplicated_find_proviruses/{sample}_contigs_dereplicated_provirus.fna",
    params:
        out_dir=results + "04_VIRUS_IDENTIFICATION/06_genomad/{sample}/",
        genomad_dir=resources + "genomad/genomad_db",
    conda:
        "../envs/genomad:1.0.1.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/genomad_{sample}.tsv"
    resources:
        runtime="05:00:00",
        mem_mb="10000",
    threads: config["virus_identification"]["genomad_threads"]
    shell:
        """
        # run genomad
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
# 07 Identify external hits
# -----------------------------------------------------
# create a mash sketch of virusdb
rule mash_sketch_virusdb:
    message:
        "Creating mash sketch of {input}"
    input:
        config["virus_db"],
    output:
        config["virus_db"] + ".msh",
    # conda:
    #     "../envs/mash:2.3--ha61e061_0.yml"
    container:
        "docker://quay.io/biocontainers/mash:2.3--ha61e061_0"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/mash_sketch_virusdb.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="500000",
    threads: config["virus_identification"]["mash_threads"]
    shell:
        """
        # create a mash sketch of virusdb
        mash sketch \
        -p {threads} \
        -i {input}
        """


# screen reads to identify external viruses
rule screen_reads_against_virusdb:
    message:
        "Screening reads against {input.sketch} to identify external viruses present in {wildcards.sample}"
    input:
        R1=R1,
        R2=R2,
        sketch=config["virus_db"] + ".msh",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/virusdb_mash_screen.tab",
    params:
        combined=results
        + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/combined_reads.fastq",
    # conda:
    #     "../envs/mash:2.3--ha61e061_0.yml"
    container:
        "docker://quay.io/biocontainers/mash:2.3--ha61e061_0"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/screen_reads_against_virusdb_sketch_{sample}.tsv"
    resources:
        runtime="01:00:00",
        mem_mb="100000",
    threads: config["virus_identification"]["mash_threads"]
    shell:
        """
        # combine reads
        cat {input.R1} {input.R2} > {params.combined}

        # screen reads against virusdb
        mash screen \
        -p {threads} \
        {input.sketch} \
        {params.combined} > {output}

        # remove combined fastq to save space
        rm {params.combined}
        """


# filter to keep only external hits
rule extract_external_hits:
    message:
        "Extracting external viruses present in {wildcards.sample}"
    input:
        read_screen=results
        + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/virusdb_mash_screen.tab",
        virusdb=config["virus_db"],
    output:
        results + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/virusdb_hits.fna",
    params:
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/extract_virusdb_hits_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/04_extract_virusdb_hits.py"


# add external hits to contigs fasta
rule combine_external_hits_with_contigs:
    message:
        "Combining external viruses with {wildcards.sample} contigs"
    input:
        mash_hits=results
        + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/virusdb_hits.fna",
        contigs=assembly,
    output:
        results
        + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/contigs_external_hits_combined.fna",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/combine_mash_hits_with_contigs_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # combine hits fasta with contigs
        cat {input.contigs} {input.mash_hits} > {output}
        """


# -----------------------------------------------------
# 07 Combine outputs
# -----------------------------------------------------
# determine input files for detecting virus sequences
if config["virus_identification"]["run_mgv"]:
    mgv1 = (
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_master_table.tsv",
    )
    mgv2 = (results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_final.tsv",)
else:
    mgv1 = pd.DataFrame()
    mgv2 = pd.DataFrame()

if config["virus_identification"]["run_virfinder"]:
    virfinder = (
        results + "04_VIRUS_IDENTIFICATION/01_mgv/output/{sample}_master_table.tsv",
    )
else:
    virfinder = pd.DataFrame()

if config["virus_identification"]["run_virsorter"]:
    virsorter = (
        results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{sample}/Metric_files/VIRSorter_phage_signal.tab"
    )
    virsorter_translation = (
        results
        + "04_VIRUS_IDENTIFICATION/02_virsorter/{sample}/fasta/input_sequences_id_translation.tsv"
    )
else:
    virsorter = pd.DataFrame()
    virsorter_translation = pd.DataFrame()

if config["virus_identification"]["run_virsorter2"]:
    virsorter2 = (
        results + "04_VIRUS_IDENTIFICATION/03_virsorter2/{sample}/final-viral-score.tsv"
    )
else:
    virsorter2 = pd.DataFrame()

if config["virus_identification"]["run_vibrant"]:
    vibrant = (
        results
        + "04_VIRUS_IDENTIFICATION/04_vibrant/{sample}/VIBRANT_{sample}_contigs_dereplicated/VIBRANT_phages_{sample}_contigs_dereplicated/{sample}_contigs_dereplicated.phages_combined.txt"
    )
else:
    vibrant = pd.DataFrame()

if config["virus_identification"]["run_deepvirfinder"]:
    deepvirfinder = (
        results
        + "04_VIRUS_IDENTIFICATION/05_deepvirfinder/{sample}_contigs_dereplicated.fasta_gt1000bp_dvfpred.txt",
    )
else:
    deepvirfinder = pd.DataFrame()

if config["virus_identification"]["run_genomad"]:
    genomad = (
        results
        + "04_VIRUS_IDENTIFICATION/06_genomad/{sample}/{sample}_contigs_dereplicated_summary/{sample}_contigs_dereplicated_virus_summary.tsv"
    )
else:
    genomad = pd.DataFrame()

if config["virus_identification"]["run_external"] and config["input_data"] == "reads":
    external = (
        results
        + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/virusdb_mash_screen.tab"
    )
else:
    external = pd.DataFrame()


# combine outputs from all tools
rule merge_reports_within_samples:
    message:
        "Merging all virus identification reports within {wildcards.sample}"
    input:
        mgv_results=mgv1,
        mgv_viruses=mgv2,
        vf_results=virfinder,
        vs_results=virsorter,
        vs_translation=virsorter_translation,
        vs2_results=virsorter2,
        vb_results=vibrant,
        dvf_results=deepvirfinder,
        genomad_results=genomad,
        external_results=external,
    output:
        results
        + "04_VIRUS_IDENTIFICATION/08_combine_outputs/{sample}/combined_report.csv",
    params:
        run_mgv=config["virus_identification"]["run_mgv"],
        run_virfinder=config["virus_identification"]["run_virfinder"],
        run_virsorter=config["virus_identification"]["run_virsorter"],
        run_virsorter2=config["virus_identification"]["run_virsorter2"],
        run_vibrant=config["virus_identification"]["run_vibrant"],
        run_deepvirfinder=config["virus_identification"]["run_deepvirfinder"],
        run_genomad=config["virus_identification"]["run_genomad"],
        run_external=config["virus_identification"]["run_external"],
        assembly="{sample}",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/merge_reports_within_samples_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/04_merge_reports_within_samples.py"


# if running external identification, use combined contigs file
if config["virus_identification"]["run_external"] and config["input_data"] == "reads":
    combined = (
        results
        + "04_VIRUS_IDENTIFICATION/07_external_hits/{sample}/contigs_external_hits_combined.fna"
    )
else:
    combined = assembly


# combine viral contigs from all tool outputs using thresholds specified in config.yaml
rule merge_viral_contigs_within_samples:
    message:
        "Merging viral contigs meeting config.yaml criteria within {wildcards.sample}"
    input:
        contigs=combined,
        viral_report=results
        + "04_VIRUS_IDENTIFICATION/08_combine_outputs/{sample}/combined_report.csv",
    output:
        results
        + "04_VIRUS_IDENTIFICATION/08_combine_outputs/{sample}/combined_viral_contigs.fasta",
    params:
        run_mgv=config["virus_identification"]["run_mgv"],
        run_vf=config["virus_identification"]["run_virfinder"],
        vf_score=config["virus_identification"]["virfinder_min_score"],
        run_vs=config["virus_identification"]["run_virsorter"],
        vs_cat=config["virus_identification"]["virsorter_cat"],
        run_vs2=config["virus_identification"]["run_virsorter2"],
        vs2_score=config["virus_identification"]["virsorter2_min_score"],
        run_dvf=config["virus_identification"]["run_deepvirfinder"],
        dvf_score=config["virus_identification"]["deepvirfinder_min_score"],
        run_vb=config["virus_identification"]["run_vibrant"],
        run_genomad=config["virus_identification"]["run_genomad"],
        genomad_score=config["virus_identification"]["genomad_min_score"],
        genomad_fdr=config["virus_identification"]["genomad_max_fdr"],
        run_external=config["virus_identification"]["run_external"],
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_hashes=config["virus_identification"]["min_mash_hashes"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
        assembly="{sample}",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/merge_viral_contigs_within_samples_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/04_merge_viral_contigs_within_samples.py"


# -----------------------------------------------------
# Analyze combined virus data
# -----------------------------------------------------
# combine virus reports across samples
rule combine_reports_across_samples:
    message:
        "Combining viral reports across all samples"
    input:
        expand(
            results
            + "04_VIRUS_IDENTIFICATION/08_combine_outputs/{sample}/combined_report.csv",
            sample=samples,
        ),
    output:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.csv",
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/combine_reports_across_samples.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# plot virus counts by tool
rule virus_identification_analysis:
    message:
        "Visualizing virus identification"
    input:
        results + "04_VIRUS_IDENTIFICATION/virus_identification_report.csv",
    output:
        boxplot_svg=report(
            results + "04_VIRUS_IDENTIFICATION/virus_identification_boxplot.svg",
            category="Step 04: Virus identification",
        ),
        boxplot_html=results
        + "04_VIRUS_IDENTIFICATION/virus_identification_boxplot.html",
        upset=report(
            results + "04_VIRUS_IDENTIFICATION/virus_identification_upsetplot.png",
            category="Step 04: Virus identification",
        ),
    params:
        run_mgv=config["virus_identification"]["run_mgv"],
        run_vf=config["virus_identification"]["run_virfinder"],
        vf_score=config["virus_identification"]["virfinder_min_score"],
        run_vs=config["virus_identification"]["run_virsorter"],
        vs_cat=config["virus_identification"]["virsorter_cat"],
        run_vs2=config["virus_identification"]["run_virsorter2"],
        vs2_score=config["virus_identification"]["virsorter2_min_score"],
        run_dvf=config["virus_identification"]["run_deepvirfinder"],
        dvf_score=config["virus_identification"]["deepvirfinder_min_score"],
        run_vb=config["virus_identification"]["run_vibrant"],
        run_genomad=config["virus_identification"]["run_genomad"],
        genomad_score=config["virus_identification"]["genomad_min_score"],
        genomad_fdr=config["virus_identification"]["genomad_max_fdr"],
        run_external=config["virus_identification"]["run_external"],
        min_mash_score=config["virus_identification"]["min_mash_score"],
        min_mash_hashes=config["virus_identification"]["min_mash_hashes"],
        min_mash_multiplicity=config["virus_identification"]["min_mash_multiplicity"],
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/04_VIRUS_IDENTIFICATION/virus_identification_analysis.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="1000",
    script:
        "../scripts/04_virus_identification_analysis.py"
