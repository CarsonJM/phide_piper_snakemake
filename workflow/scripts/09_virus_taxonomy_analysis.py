import pandas as pd
import plotly.express as px

# load genomad taxonomy
genomad = pd.read_csv(str(snakemake.input.genomad), sep='\t')
genomad = genomad.add_prefix('genomad_')
genomad_filt = genomad[genomad['genomad_agreement'] > snakemake.params.min_genomad_agreement]
genomad_filt['seq_name'] = genomad_filt['genomad_seq_name']

# load mmseqs taxnomy
mmseqs = pd.read_csv(str(snakemake.input.mmseqs), sep='\t', names=['seq_name', 'mmseqs_taxid', 'mmseqs_rank', 'mmseqs_taxa', 'mmseqs_total_proteins', 'mmseqs_proteins_labeled', 'mmseqs_proteins_agreeing', 'mmseqs_agreement', 'mmseqs_lineage'])
mmseqs_filt = mmseqs[mmseqs['mmseqs_agreement'] > snakemake.params.min_mmseqs_agreement]

# merge and save
merged = mmseqs.merge(genomad_filt, on='seq_name', how='outer')
merged.to_csv(str(snakemake.output.report), sep='\t', index=False)

# generate figure of taxonomy outputs
merged['GeNomad'] = merged.apply(lambda x: str(x.genomad_lineage).split(';')[4] if str(x.genomad_lineage).count(';') > 2 else 'Unknown', axis=1)
merged['MMSeqs2'] = merged.apply(lambda x: str(x.mmseqs_lineage).split(';')[3] if str(x.mmseqs_lineage).count(';') > 2 else 'Unknown', axis=1)
merged['MMSeqs2'] = merged.apply(lambda x: str(x.MMSeqs2).split('_')[1] if str(x.mmseqs_lineage).count('_') > 0 else 'Unknown', axis=1)

merged_melt = merged.melt(id_vars=['seq_name'], value_vars=['MMSeqs2', 'GeNomad'])
genomad_group = merged.groupby('GeNomad', as_index=False).count()
mmseqs_group = merged.groupby('MMSeqs2', as_index=False).count()
plot = merged_melt.groupby(['value', 'variable'], as_index=False).count()
plot['proportion'] = plot['seq_name']/sum(plot['seq_name'])
ordered = plot.sort_values('value').drop_duplicates('value')
plot['value'] = pd.Categorical(plot['value'], ordered['value'])
plot.rename(columns={"value": "Order"}, inplace=True)

# plot kneaddata read counts and save
fig = px.bar(plot, x='variable', y='seq_name',  color='Order',
             labels={
                     "variable": "Taxonomy method",
                     "seq_name" : "Proportion of viruses"
                 })

fig.update_layout(title_text='Virus taxonomic assignment')
fig.write_html(str(snakemake.output.html))
fig.write_image(str(snakemake.output.svg))
