import pandas as pd

# load iphop data table
iphop = pd.read_csv(str(snakemake.input.iphop))
iphop.rename(columns = {'Virus':'viral_genome'}, inplace=True)

# load phist data table
phist = pd.read_csv(str(snakemake.input.phist))

# merge phist and iphop
merged = iphop.merge(phist, on='viral_genome', how='left')
merged.to_csv(str(snakemake.output), sep=',', index=False)