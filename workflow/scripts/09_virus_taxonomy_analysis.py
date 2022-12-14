import pandas as pd
import plotly.express as px

# load genomad taxonomy
genomad = pd.read_csv(str(snakemake.input.genomad), sep='\t')
genomad_filt = genomad[(genomad['agreement'] > snakemake.params.min_genomad_agreement) & (genomad['n_genes_with_taxonomy'] > snakemake.params.min_genomad_genes)]

# merge and save
genomad_filt.to_csv(str(snakemake.output.report), sep='\t', index=False)

# generate figure of taxonomy outputs
genomad_filt['GeNomad'] = genomad_filt.apply(lambda x: str(x.lineage).split(';')[4] if str(x.lineage).count(';') > 2 else 'Unknown', axis=1)

genomad_melt = genomad_filt.melt(id_vars=['seq_name'], value_vars=['GeNomad'])
plot = genomad_melt.groupby(['value', 'variable'], as_index=False).count()
plot['proportion'] = plot['seq_name']*2/sum(plot['seq_name'])
ordered = plot.sort_values('value').drop_duplicates('value')
plot['value'] = pd.Categorical(plot['value'], ordered['value'])
plot.rename(columns={"value": "Order"}, inplace=True)

# plot kneaddata read counts and save
fig = px.bar(plot, x='variable', y='proportion',  color='Order',
             labels={
                     "variable": "Taxonomy method",
                     "proportion" : "Proportion of viruses"
                 })

fig.update_layout(title_text='Virus taxonomic assignment')
fig.write_image(str(snakemake.output.svg))
