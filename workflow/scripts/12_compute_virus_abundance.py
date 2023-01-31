import pandas as pd

# load instrain output
genome_present = pd.read_csv(str(snakemake.input), sep='\t')

# filter to retain only present viruses
# if snakemake.params.recover_low_abundance:
#     genome_present = genome[(genome['breadth_expected'] >= snakemake.params.min_breadth) & (genome['breadth_minCov']/genome['breadth_expected'] >= snakemake.params.min_breadth)]
# else:
#     genome_present = genome[genome['breadth_minCov'] >= snakemake.params.min_breadth]

# calculate relative abundance using coverage
total_coverage = sum(genome_present['coverage'])
genome_present['relative_abundance'] = genome_present['coverage']/total_coverage
genome_present.to_csv(str(snakemake.output), sep='\t', index=False)
