import pandas as pd

combined_df = pd.DataFrame()

for file in snakemake.input:
    df = pd.read_csv(file, sep='\t')
    combined_df = combined_df.append(df)

combined_df.to_csv(str(snakemake.output), index=False, sep='\t')
