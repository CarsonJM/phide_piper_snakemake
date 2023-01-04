#
# rule combine_propagate_results:
#     output:
#         resources + "propagate/propagate_build_complete",
#     params:
#         propagate_dir=resources + "propagate",
#     conda:
#         "../envs/propagate:1.1.0.yml"
#     benchmark:
#         "benchmark/10_VIRUS_LIFESTYLE/build_propagate.tsv"
#     resources:
#         runtime=config["virus_lifestyle"]["bacphlip_runtime"],
#         mem_mb=config["virus_lifestyle"]["bacphlip_memory"],
