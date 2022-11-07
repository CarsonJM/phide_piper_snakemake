# import genomad taxonomy
gn_taxonomy = pd.read_csv(str(snakemake.input), sep='\t')

# filter taxonomy based on genes and agreement
gn_taxonomy_filtered = gn_taxonomy[(gn_taxonomy['n_genes_with_taxonomy'] > snakemake.params.genomad_min_genes) & (gn_taxonomy['agreement'] > snakemake.params.genomad_min_agreement)]

# 