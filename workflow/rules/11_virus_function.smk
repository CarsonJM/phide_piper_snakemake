# -----------------------------------------------------
# Virus function
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples_df"], sep="\t")
samples = samples_df['sample']


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -----------------------------------------------------
# Virus function rules
# -----------------------------------------------------


# -----------------------------------------------------
# 01 DRAMv
# -----------------------------------------------------
# build dram
rule build_dramv:
    output:
        test="test"
    params:
        dram_dir=resources + "dram/DRAM_data",
    conda:
        "../envs/dram:1.3.4--pyhdfd78af_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/dram:1.3.4--pyhdfd78af_0"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/build_dramv.tsv"
    resources:
        runtime="48:00:00",
        mem_mb="100000",
        partition="compute-hugemem",
        account="pedslabs"
    threads: config["virus_function"]["dramv_threads"]
    shell:
        """
        DRAM-setup.py prepare_databases --skip_uniref --output_dir {params.dram_dir}

        touch {output.test}
        """


rule update_dram:
    input:
        test="test"
    output:
        test="update_test"
    params:
        dram_dir=resources + "dram/DRAM_data",
    conda:
        "../envs/dram:1.3.4--pyhdfd78af_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/dram:1.3.4--pyhdfd78af_0"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/update_dram.tsv"
    resources:
        runtime="48:00:00",
        mem_mb="100000",
        partition="compute-hugemem",
        account="pedslabs",
    threads: config["virus_function"]["dramv_threads"]
    shell:
        """
        DRAM-setup.py update_description_db

        touch {output.test}
        """


# run virsorter2 to identifiy viral contigs
rule virsorter2_dram:
    input:
        contigs=results + "07_VIRUS_DIVERSITY/01_votu_clusters/votu_representatives.fna",
        vs2_db=resources + "virsorter2/Done_all_setup",
    output:
        tab=results
        + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/viral-affi-contigs-for-dramv.tab",
        fasta=results + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/final-viral-combined-for-dramv.fa",
    params:
        vs2_db=resources + "virsorter2",
        out_dir=results + "11_VIRUS_FUNCTION/01_virsorter2/",
        extra_args=config["virus_identification"]["virsorter2_arguments"],
    conda:
        "../envs/virsorter:2.2.3--pyhdfd78af_1.yml"
    threads: config["virus_identification"]["virsorter2_threads"]
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/virsorter2.tsv"
    resources:
        runtime="04:00:00",
        mem_mb="10000",
    shell:
        """
        rm -rf {params.out_dir}

        # run virsorter2
        virsorter run all --keep-original-seq \
        --db-dir {params.vs2_db} \
        -w {params.out_dir} \
        -i {input.contigs} \
        -j {threads} \
        --prep-for-dramv \
        --min-score 0.0 \
        --rm-tmpdir \
        --viral-gene-enrich-off \
        {params.extra_args}
        """


# run dramv
rule dramv:
    input:
        test="update_test",
        vs2=results + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/viral-affi-contigs-for-dramv.tab",
        viruses=results + "11_VIRUS_FUNCTION/01_virsorter2/for-dramv/final-viral-combined-for-dramv.fa",
    output:
        results + "dramv_test"
    params:
        out_dir = results + "11_VIRUS_FUNCTION/02_dramv/",
    conda:
        "../envs/dram:1.3.4--pyhdfd78af_0.yml"
    # container:
    #     "docker://quay.io/biocontainers/dram:1.3.4--pyhdfd78af_0"
    threads: config["virus_function"]["dramv_threads"]
    shell:
        """
        DRAM-v.py annotate \
        -i {input.viruses} \
        -v {input.vs2} \
        --verbose \
        --keep_tmp_dir \
        -o {params.out_dir} \
        --threads {threads}
        """


# -----------------------------------------------------
# 02 geNomad function
# -----------------------------------------------------
# add sample column to genomad output
rule add_sample_to_genomad_functions:
    input:
        results
        + "04_VIRUS_IDENTIFICATION/06_genomad/{sample}/{sample}_scaffolds_dereplicated_annotate/{sample}_scaffolds_dereplicated_genes.tsv",
    output:
        results + "11_VIRUS_FUNCTION/02_genomad/{sample}/gene_functions.tsv"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/add_sample_to_genomad_functions_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # add sample column to each checkv output
        s={wildcards.sample}
        sed -i "s/$/\t$s/" {input}
        sample="sample"
        sed -i "1s/$s/$sample/" {input}

        cat {input} > {output}
        """


# combine taxonomies across samples
rule combine_genomad_functions:
    input:
        expand(results + "11_VIRUS_FUNCTION/02_genomad/{sample}/gene_functions.tsv", sample=samples)
    output:
        results + "11_VIRUS_FUNCTION/02_genomad/gene_functions_combined.tsv",
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/combine_genomad_functions.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# filter genomad function to contain only viruses
rule filter_genomad_functions:
    input:
        results + "11_VIRUS_FUNCTION/02_genomad/gene_functions_combined.tsv",
    output:
        results + "11_VIRUS_FUNCTION/02_genomad/gene_functions_filtered.tsv",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/filter_genomad_functions.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "11_filter_genomad_functions.py"


# -----------------------------------------------------
# 03 VIBRANT
# -----------------------------------------------------
# add sample column to genomad output
rule add_sample_to_vibrant_functions:
    input:
        results + "04_VIRUS_IDENTIFICATION/04_vibrant/{sample}/VIBRANT_{sample}_contigs_dereplicated/VIBRANT_results_{sample}_contigs_dereplicated/VIBRANT_annotations_oral_enriched_contigs_dereplicated.tsv"
    output:
        results + "11_VIRUS_FUNCTION/03_vibrant/{sample}/gene_functions.tsv"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/add_sample_to_vibrant_functions_{sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # add sample column to each checkv output
        s={wildcards.sample}
        sed -i "s/$/\t$s/" {input}
        sample="sample"
        sed -i "1s/$s/$sample/" {input}

        cat {input} > {output}
        """


# combine taxonomies across samples
rule combine_vibrant_functions:
    input:
        expand(results + "11_VIRUS_FUNCTION/03_vibrant/{sample}/gene_functions.tsv", sample=samples)
    output:
        results + "11_VIRUS_FUNCTION/03_vibrant/gene_functions_combined.tsv",
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/combine_vibrant_functions.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """


# filter genomad function to contain only viruses
rule filter_vibrant_functions:
    input:
        results + "11_VIRUS_FUNCTION/03_vibrant/gene_functions_combined.tsv",
    output:
        results + "11_VIRUS_FUNCTION/03_vibrant/gene_functions_filtered.tsv",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/11_VIRUS_FUNCTION/filter_vibrant_functions.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "11_filter_vibrant_functions.py"