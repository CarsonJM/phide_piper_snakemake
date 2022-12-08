import pandas as pd

# load ani results
ani = pd.read_csv(str(snakemake.input.ani), sep='\t')
ani['accession'] = ani['tname'].str.split('.', expand=True)[0]

# load virus metadata
metadata = pd.read_csv(str(snakemake.input.meta))

# merge mash hits and metadata
ani_meta = ani.merge(metadata, on='accession', how='left')

# identify top hits
ani_id95 = ani_meta[ani_meta['pid'] >= snakemake.params.min_ani]
ani_id95_sort = ani_id95.sort_values('qcov', ascending=False)
ani_id95_top_qcov = ani_id95_sort.groupby('qname', as_index=False).first()

# identify species-level hits
ani_id95_top_qcov['novel_species'] = ani_id95_top_qcov.apply(lambda x: False if x.pid >= snakemake.params.min_ani and x.qcov >= snakemake.params.min_qcov and x.tcov >= snakemake.params.min_tcov else True, axis=1)

ani_id95_top_qcov.to_csv(str(snakemake.output), sep='\t', index=False)