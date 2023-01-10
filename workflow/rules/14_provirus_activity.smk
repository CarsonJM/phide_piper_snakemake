# -----------------------------------------------------
# Provirus activity Module
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"


samples_df = pd.read_csv(config["samples_df"], sep="\t")
group_sample_dictionary = samples_df.groupby(["group"])["sample"].apply(list).to_dict()

# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "../report/workflow.rst"


# -----------------------------------------------------
# Provirus activity rules
# -----------------------------------------------------
localrules:
    identify_integrated_prophages,
    build_propagate,


# -----------------------------------------------------
# 00 Determine inputs for module
# -----------------------------------------------------
if (
    config["input_data"] == "reads"
    or config["input_data"] == "contigs"
    or config["input_data"] == "vls"
):
    viruses = (
        results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps_viruses.fasta",
    )
elif config["input_data"] == "viruses":
    viruses = results + "00_INPUT/{group_sample}_viruses.fasta"


# read symlink in 04_virus_identification.smk
if config["input_data"] == "reads":
    R1 = results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_1.fastq.gz"
    R2 = results + "01_READ_PREPROCESSING/03_kneaddata/{group_sample}_paired_2.fastq.gz"
elif (
    config["input_data"] == "contigs"
    or config["input_data"] == "viruses"
    or config["input_data"] == "processed_viruses"
):
    R1 = results + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R1.fastq.gz"
    R2 = results + "00_INPUT/01_merge_repicates/{group_sample}.preprocessed_R2.fastq.gz"


# -----------------------------------------------------
# 01 PropagAtE
# -----------------------------------------------------
rule build_propagate:
    output:
        resources + "propagate/example_output/dormant/test_run/test_run.tsv",
    params:
        propagate_dir=resources + "propagate",
    conda:
        "../envs/propagate:1.1.0.yml"
    benchmark:
        "benchmark/14_PROVIRUS_ACTIVITY/build_propagate.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    shell:
        """
        # clone propagate
        rm -rf {params.propagate_dir}
        git clone https://github.com/CarsonJM/PropagAtE_NaN.git {params.propagate_dir}

        # pip install propagate
        cd {params.propagate_dir}
        pip install .

        # test propagate
        cd example_output/dormant
        Propagate -f example_sequence.fasta -r sample_forward_reads.fastq.gz sample_reverse_reads.fastq.gz -v VIBRANT_integrated_prophage_coordinates_example.tsv -o test_run --clean -t {threads}
        """


rule identify_integrated_prophages:
    input:
        profiles=lambda wildcards: expand(
            results
            + "12_VIRUS_ABUNDANCE/02_instrain_profile/{{group}}_{sample}/output/{{group}}_{sample}_genome_info_w_abundance.tsv",
            sample=group_sample_dictionary[wildcards.group],
        ),
        votu_members=results + "07_VIRUS_DIVERSITY/01_votu_clustering/votu_clusters.tsv",
        derep_untrimmed=results
        + "06_VIRUS_DEREPLICATION/02_dereplicate_viruses/dereplicate_reps_untrimmed.fasta",
    output:
        report=results
        + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/integrated_prophage_report.tsv",
        coords=results
        + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/integrated_prophage_coords.tsv",
        sequences=results
        + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/untrimmed_prophage_sequences.fasta",
    params:
        group="{group}",
    conda:
        "../envs/jupyter.yml"
    benchmark:
        "benchmark/14_PROVIRUS_ACTIVITY/identify_integrated_prophages_{group}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    script:
        "../scripts/14_identify_integrated_prophages.py"


rule propagate:
    input:
        resources + "propagate/example_output/dormant/test_run/test_run.tsv",
        R1=R1,
        R2=R2,
        sequences=results
        + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/untrimmed_prophage_sequences.fasta",
        coords=results
        + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/integrated_prophage_coords.tsv",
    output:
        results
        + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/{group_sample}/{group_sample}.tsv",
    params:
        out_dir=results + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/{group_sample}",
        min_id=config["provirus_activity"]["min_id"],
        min_breadth=config["provirus_activity"]["min_breadth"],
        min_depth=config["provirus_activity"]["min_depth"],
        extra_args=config["provirus_activity"]["propagate_arguments"],
    conda:
        "../envs/propagate:1.1.0.yml"
    benchmark:
        "benchmark/14_PROVIRUS_ACTIVITY/propagate_{group}_{group_sample}.tsv"
    resources:
        runtime="00:10:00",
        mem_mb="10000",
    threads: config["provirus_activity"]["propagate_threads"]
    shell:
        """
        rm -rf {params.out_dir}

        Propagate -f {input.sequences} \
        -v {input.coords} \
        -r {input.R1} {input.R2} \
        -o {params.out_dir} \
        -t {threads} \
        -p {params.min_id} \
        --min {params.min_depth} \
        --breadth {params.min_breadth} \
        {params.extra_args}

        # add group column to each propagate output
        s={wildcards.group}
        sed -i "s/$/\t$s/" {output}
        group="group"
        sed -i "1s/$s/$group/" {output}

        # add sample column to each propagate output
        s={wildcards.group_sample}
        sed -i "s/$/\t$s/" {output}
        sample="sample"
        sed -i "1s/$s/$sample/" {output}
        """


rule combine_propagate_within_groups:
    input:
        lambda wildcards: expand(
            results
            + "14_PROVIRUS_ACTIVITY/01_propagate/{group}/{group_sample}/{group_sample}.tsv",
            group_sample=group_to_group_sample_dictionary[wildcards.group],
            group=wildcards.group,
        ),
    output:
        results + "14_PROVIRUS_ACTIVITY/01_propagate/{group}_combined_output.tsv",
    benchmark:
        "benchmark/14_PROVIRUS_ACTIVITY/propagate_{group}.tsv"
    resources:
        runtime=config["provirus_activity"]["propagate_runtime"],
        mem_mb=config["provirus_activity"]["propagate_memory"],
    shell:
        """
        # combine all outputs, only keeping header from one file
        awk 'FNR>1 || NR==1' {input} > {output}
        """
