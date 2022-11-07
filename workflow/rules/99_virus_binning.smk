# -----------------------------------------------------
# Virus Binning
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv("config/samples.tsv", sep="\t")
samples = samples_df['sample']
coassemblies = list(set(samples_df['assembly']))


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -----------------------------------------------------
# Virus binning rules
# -----------------------------------------------------
# -----------------------------------------------------
# 01 Combine checkV viruses
# -----------------------------------------------------
rule combine_checkv_output_single:
    input:
        checkv_proviruses=results + "07_VIRUS_QUALITY/01_checkv/single/proviruses.fna",
        checkv_viruses=results + "07_VIRUS_QUALITY/01_checkv/single/viruses.fna",
    output:
        results + "08_VIRUS_BINNING/01_align_viruses/single/checkv_viruses.fna",
    shell:
        """
        cat {input.checkv_proviruses} {input.checkv_viruses} > {output}
        """


rule combine_checkv_output_coassembly:
    input:
        checkv_proviruses=results + "07_VIRUS_QUALITY/01_checkv/coassembly/proviruses.fna",
        checkv_viruses=results + "07_VIRUS_QUALITY/01_checkv/coassembly/viruses.fna",
    output:
        results + "08_VIRUS_BINNING/01_align_viruses/coassembly/checkv_viruses.fna",
    shell:
        """
        cat {input.checkv_proviruses} {input.checkv_viruses} > {output}
        """

# Align reads to virus catalog using bowtie2
rule build_virus_bowtie2db_single:
    input:
        results + "08_VIRUS_BINNING/01_align_viruses/single/checkv_viruses.fna",
    output:
        results + "08_VIRUS_BINNING/01_align_viruses/single/virus_catalog.1.bt2",
    params:
        db=results + "08_VIRUS_BINNING/01_align_viruses/single/virus_catalog",
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input} {params.db} --threads {threads}
        """


# Align reads to virus catalog using bowtie2
rule virus_bowtie2_single:
    input:
        R1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
        db=results + "08_VIRUS_BINNING/01_align_viruses/single/virus_catalog.1.bt2",
    output:
        results
        + "08_VIRUS_BINNING/01_align_viruses/single/sam_files/{sample}.sam",
    params:
        db=results + "08_VIRUS_BINNING/01_align_viruses/single/virus_catalog",
        extra_args=config["virus_binning"]["bowtie2_arguments"],
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {output}
        """

rule identify_viruses_to_bin_single:
    input:
        checkv_viruses=results + "08_VIRUS_BINNING/01_align_viruses/single/checkv_viruses.fna",
        checkv_results=results + "07_VIRUS_QUALITY/01_checkv/single/quality_summary.tsv",
    output:
        results + "08_VIRUS_BINNING/01_align_viruses/single/viruses_to_bin.csv",
    params:
        max_completeness_to_bin=config["virus_binning"]["max_completeness_to_bin"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/08_identify_viruses_to_bin.py.ipynb"


# Align reads to virus catalog using bowtie2
rule build_virus_bowtie2db_coassembly:
    input:
        results + "08_VIRUS_BINNING/01_align_viruses/coassembly/checkv_viruses.fna",
    output:
        results + "08_VIRUS_BINNING/01_align_viruses/coassembly/virus_catalog.1.bt2",
    params:
        db=results + "08_VIRUS_BINNING/01_align_viruses/coassembly/virus_catalog",
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # make a bowtie2 db from virusdb
        bowtie2-build {input} {params.db} --threads {threads}
        """


# Align reads to virus catalog using bowtie2
rule virus_bowtie2_coassembly:
    input:
        R1=results
        + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_1.fastq",
        R2=results
        + "01_READ_PREPROCESSING/03_kneaddata/{sample}_paired_2.fastq",
        db=results + "08_VIRUS_BINNING/01_align_viruses/coassembly/virus_catalog.1.bt2",
    output:
        results
        + "08_VIRUS_BINNING/01_align_viruses/coassembly/sam_files/{sample}.sam",
    params:
        db=results + "08_VIRUS_BINNING/01_align_viruses/coassembly/virus_catalog",
        extra_args=config["virus_binning"]["bowtie2_arguments"],
    conda:
        "../envs/kneaddata.yml"
    threads: config["virus_abundance"]["metapop_threads"]
    shell:
        """
        # align reads to bowtie2 database
        bowtie2 \
        --threads {threads} \
        -x {params.db} \
        -1 {input.R1} \
        -2 {input.R2} \
        -S {output}
        """

rule identify_viruses_to_bin_coassembly:
    input:
        checkv_viruses=results + "08_VIRUS_BINNING/01_align_viruses/coassembly/checkv_viruses.fna",
        checkv_results=results + "07_VIRUS_QUALITY/01_checkv/coassembly/quality_summary.tsv",
    output:
        results + "08_VIRUS_BINNING/01_align_viruses/coassembly/viruses_to_bin.csv",
    params:
        max_completeness_to_bin=config["virus_binning"]["max_completeness_to_bin"],
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/08_identify_viruses_to_bin.py.ipynb"


# -----------------------------------------------------
# 02 vRhyme
# -----------------------------------------------------
# build vrhyme
rule build_vrhyme:
    output:
        resources + "vrhyme/vRhyme/vRhyme",
    params:
        vrhyme_dir=resources + "vrhyme",
        vrhyme_models=resources + "vrhyme/vRhyme/models/vRhyme_machine_model_ET.sav.gz"
    conda:
        "../envs/vrhyme.yml"
    shell:
        """
        # clone vhryme repo
        rm -rf {params.vrhyme_dir}
        git clone https://github.com/AnantharamanLab/vRhyme {params.vrhyme_dir}

        # unzip models for binning
        gunzip {params.vrhyme_models}

        # change properties for vrhyme
        cd {params.vrhyme_dir}/vRhyme
        chmod +x vRhyme scripts/*.py aux/*.py

        # test vhryme
        aux/test_vRhyme.py
        """


# build vrhyme
rule vrhyme_single:
    input:
        vrhyme_script=resources + "vrhyme/vRhyme/vRhyme",
        catalog=results + "08_VIRUS_BINNING/01_align_viruses/single/checkv_viruses.fna",
        bam_files=expand(results
        + "08_VIRUS_BINNING/01_align_viruses/single/sam_files/{sample}.sam", sample=samples),
        viruses_to_bin=results + "08_VIRUS_BINNING/01_align_viruses/single/viruses_to_bin.csv",
    output:
        vmags=results + "08_VIRUS_BINNING/02_vrhyme/single/vrhyme_vmags.fasta",
        unbinned=results + "08_VIRUS_BINNING/02_vrhyme/single/vrhyme_unbinned.fasta",
    params:
        vrhyme_dir=resources + "vrhyme",
        sam_dir=results
        + "08_VIRUS_BINNING/01_align_viruses/single/sam_files/*.sam",
        out_dir=results + "08_VIRUS_BINNING/02_vrhyme/single/",
        extra_args=config["virus_binning"]["vrhyme_arguments"],
        fasta_dir=results + "08_VIRUS_BINNING/02_vrhyme/single/vRhyme_best_bins_fasta/",
        linked_fasta_dir=results + "08_VIRUS_BINNING/02_vrhyme/single/vRhyme_linked_fasta_bins/",
    conda:
        "../envs/vrhyme.yml"
    threads:
        config["virus_binning"]["vrhyme_threads"]
    shell:
        """
        # run vrhyme
        rm -rf {params.out_dir}
        cd {params.vrhyme_dir}
        python vRhyme/vRhyme -i {input.catalog} \
        -s {params.sam_dir} \
        -o {params.out_dir} \
        -t {threads} \
        --verbose \
        --interest {input.viruses_to_bin} \
        {params.extra_args}

        python vRhyme/aux/link_bin_sequences.py \
        -i {params.fasta_dir} \
        -o {params.linked_fasta_dir} \
        -e fasta \
        -n 1000

        cat {params.linked_fasta_dir}*.fasta > {output.vmags}

        python vRhyme/aux/extract_unbinned_sequences.py \
        -i {params.out_dir}vRhyme_best_bins.*.membership.tsv \
        -f {input.catalog} \
        -o {output.unbinned}
        """


# build vrhyme
rule vrhyme_coassembly:
    input:
        vrhyme_script=resources + "vrhyme/vRhyme/vRhyme",
        catalog=results + "08_VIRUS_BINNING/01_align_viruses/coassembly/checkv_viruses.fna",
        bam_files=expand(results
        + "08_VIRUS_BINNING/01_align_viruses/coassembly/sam_files/{sample}.sam", sample=samples),
        viruses_to_bin=results + "08_VIRUS_BINNING/01_align_viruses/coassembly/viruses_to_bin.csv",
    output:
        vmags=results + "08_VIRUS_BINNING/02_vrhyme/coassembly/vrhyme_vmags.fasta",
        unbinned=results + "08_VIRUS_BINNING/02_vrhyme/coassembly/vrhyme_unbinned.fasta",
    params:
        vrhyme_dir=resources + "vrhyme",
        sam_dir=results
        + "08_VIRUS_BINNING/01_align_viruses/coassembly/sam_files/*.sam",
        out_dir=results + "08_VIRUS_BINNING/02_vrhyme/coassembly/",
        extra_args=config["virus_binning"]["vrhyme_arguments"],
        fasta_dir=results + "08_VIRUS_BINNING/02_vrhyme/coassembly/vRhyme_best_bins_fasta/",
        linked_fasta_dir=results + "08_VIRUS_BINNING/02_vrhyme/coassembly/vRhyme_linked_fasta_bins/",
    conda:
        "../envs/vrhyme.yml"
    threads:
        config["virus_binning"]["vrhyme_threads"]
    shell:
        """
        # run vrhyme
        rm -rf {params.out_dir}
        cd {params.vrhyme_dir}
        python vRhyme/vRhyme -i {input.catalog} \
        -s {params.sam_dir} \
        -o {params.out_dir} \
        -t {threads} \
        --verbose \
        --interest {input.viruses_to_bin} \
        {params.extra_args}

        python vRhyme/aux/link_bin_sequences.py \
        -i {params.fasta_dir} \
        -o {params.linked_fasta_dir} \
        -e fasta \
        -n 1000

        cat {params.linked_fasta_dir}*.fasta > {output.vmags}

        python vRhyme/aux/extract_unbinned_sequences.py \
        -i {params.out_dir}vRhyme_best_bins.*.membership.tsv \
        -f {input.catalog} \
        -o {output.unbinned}
        """


# -----------------------------------------------------
# 03 CheckV
# -----------------------------------------------------
# run checkv on viral contigs to determine genome quality
rule vmag_checkv_single:
    input:
        checkv_db=resources + "checkv/checkv-db-v1.2/README.txt",
        vmags=results + "08_VIRUS_BINNING/02_vrhyme/single/vrhyme_vmags.fasta",
        unbinned=results + "08_VIRUS_BINNING/02_vrhyme/single/vrhyme_unbinned.fasta",
    output:
        checkv_results=results + "08_VIRUS_BINNING/03_checkv/single/quality_summary.tsv",
        checkv_proviruses=results + "08_VIRUS_BINNING/03_checkv/single/proviruses.fna",
        checkv_viruses=results + "08_VIRUS_BINNING/03_checkv/single/viruses.fna",
    params:
        combined_fasta=results + "08_VIRUS_BINNING/03_checkv/single/vmags_and_unbinned.fasta",
        checkv_dir=results + "08_VIRUS_BINNING/03_checkv/single",
        checkv_db=resources + "checkv/checkv-db-v1.2",
    log:
        results + "00_LOGS/07_virus_binning.checkv_single.log",
    conda:
        "../envs/checkv.yml"
    threads: config["virus_quality"]["checkv_threads"]
    shell:
        """
        cat {input.vmags} {input.unbinned} > {params.combined_fasta}

        # run checkv to determine virus quality
        checkv end_to_end {params.combined_fasta} {params.checkv_dir} \
        -d {params.checkv_db} \
        -t {threads} > {log} 2>&1
        """

# run checkv on viral contigs to determine genome quality
rule vmag_checkv_coassembly:
    input:
        checkv_db=resources + "checkv/checkv-db-v1.2/README.txt",
        vmags=results + "08_VIRUS_BINNING/02_vrhyme/coassembly/vrhyme_vmags.fasta",
        unbinned=results + "08_VIRUS_BINNING/02_vrhyme/coassembly/vrhyme_unbinned.fasta",
    output:
        checkv_results=results + "08_VIRUS_BINNING/03_checkv/coassembly/quality_summary.tsv",
        checkv_proviruses=results + "08_VIRUS_BINNING/03_checkv/coassembly/proviruses.fna",
        checkv_viruses=results + "08_VIRUS_BINNING/03_checkv/coassembly/viruses.fna",
    params:
        combined_fasta=results + "08_VIRUS_BINNING/03_checkv/coassembly/vmags_and_unbinned.fasta",
        checkv_dir=results + "08_VIRUS_BINNING/03_checkv/coassembly",
        checkv_db=resources + "checkv/checkv-db-v1.2",
    log:
        results + "00_LOGS/08_virus_binning.checkv_coassembly.log",
    conda:
        "../envs/checkv.yml"
    threads: config["virus_quality"]["checkv_threads"]
    shell:
        """
        cat {input.vmags} {input.unbinned} > {params.combined_fasta}

        # run checkv to determine virus quality
        checkv end_to_end {params.combined_fasta} {params.checkv_dir} \
        -d {params.checkv_db} \
        -t {threads} > {log} 2>&1
        """

# -----------------------------------------------------
# Binning analysis
# -----------------------------------------------------
# rule virus_binning_analysis:
#     input:
#         single=results + "08_VIRUS_BINNING/01_checkv/single/quality_summary.tsv",
#         coassembly=results + "08_VIRUS_BINNING/01_checkv/coassembly/quality_summary.tsv",
#         single_bin=results + "08_VIRUS_BINNING/03_checkv/single/quality_summary.tsv",
#         coassembly_bin=results + "08_VIRUS_BINNING/03_checkv/coassembly/quality_summary.tsv",
#     output:
#         report(
#             results + "07_VIRUS_QUALITY/virus_quality_figure.png",
#             caption="../report/06_virus_quality_analysis.rst",
#             category="Step 06: Virus quality",
#         ),
#     conda:
#         "../envs/jupyter.yml"
#     notebook:
#         "../notebooks/06_virus_quality_analysis.py.ipynb"