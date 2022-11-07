# -----------------------------------------------------
# Spacer prediction
# -----------------------------------------------------
import pandas as pd


# Load sample information and validate
configfile: "config/config.yaml"
samples_df = pd.read_csv(config["samples_df"], sep="\t")


# load results path
results = config["results"]


# load resources path
resources = config["resources"]


# load report
report: "report/workflow.rst"


# -----------------------------------------------------
# Spacer prediction rules
# -----------------------------------------------------
# -----------------------------------------------------
# 01 CRT & PILER-CR
# -----------------------------------------------------
rule mgv_identify_crispr:
    input:
        resources + "mgv/viral_detection_pipeline/input/imgvr.hmm",
    output:
        spacers=results + "13_SPACER_PREDICTION/01_crt_piler/merged.spacers",
        arrays=results + "13_SPACER_PREDICTION/01_crt_piler/merged.arrays",
    params:
        out_dir=results + "13_SPACER_PREDICTION/01_crt_piler/",
        identify_crispr=resources + "mgv/crispr_spacers/identify_crispr.py",
        merge_crispr=resources + "mgv/crispr_spacers/merge_crispr.py",
        bacteria_db_dir=config['bacteria_db_dir'],
        spacer_dir = resources + "mgv/crispr_spacers/",
    conda:
        "../envs/mgv_spacer.yml"
    shell:
        """
        # change to mgv_dir
        cd {params.spacer_dir}

        touch {output.arrays}
        touch {output.spacers}

        for file in {params.bacteria_db_dir}/*.fasta
        do
        file_suffix=$(sed 's%{params.bacteria_db_dir}/%%g' <<< $file)
        python {params.identify_crispr} -i $file -o {params.out_dir}
        python {params.merge_crispr} {params.out_dir}crt {params.out_dir}pilercr {params.out_dir}/$file_suffix
        tail -n +2 {params.out_dir}"$file_suffix".arrays >> {output.arrays}
        tail -n +2 {params.out_dir}"$file_suffix".spacers >> {output.spacers}
        rm {params.out_dir}crt* {params.out_dir}pilercr* {params.out_dir}$file_suffix*
        done
        """

# -----------------------------------------------------
# Combine CRISPR predictions
# -----------------------------------------------------
rule combine_crispr_spacers:
    input:
        piler_crt_spacers=results + "13_SPACER_PREDICTION/01_crt_piler/merged.spacers",
        piler_crt_arrays=results + "13_SPACER_PREDICTION/01_crt_piler/merged.arrays",
    output:
        results + "13_SPACER_PREDICTION/02_combined_spacers/combined_spacers.fasta",
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/13_combine_predicted_spacers.py.ipynb"


rule add_crispr_spacers_to_spacer_db:
    input:
        results + "13_SPACER_PREDICTION/02_combined_spacers/combined_spacers.fasta",
    output:
        config["spacer_db"]
    conda:
        "../envs/jupyter.yml"
    notebook:
        "../notebooks/13_combine_predicted_spacers.py.ipynb"