import pandas as pd
import plotly.express as px

# load iphop data table
iphop = pd.read_csv(str(snakemake.input))
# iphop = pd.read_csv('/home/carsonjm/results/phide_piper_test/08_VIRUS_HOST/01_iphop/Host_prediction_to_genus_m90.csv')
iphop.rename(columns = {'Virus':'viral_genome', 'Host genus':'iphop_genus'}, inplace=True)

# merge phist and iphop
iphop.to_csv(str(snakemake.output.report), sep='\t', index=False)

# generate figure of taxonomy outputs
iphop['iPHoP'] = iphop.apply(lambda x: str(x.iphop_genus).split(';')[1] if str(x.iphop_genus).count(';') > 0 else 'Unknown', axis=1)
iphop['iPHoP'] = iphop.apply(lambda x: str(x.iPHoP).split('p__')[1] if str(x.iPHoP).count('p__') > 0 else 'Unknown', axis=1)
iphop['iPHoP'] = iphop.apply(lambda x: str(x.iPHoP).split('_')[0] if str(x.iPHoP).count('_') > 0 else x.iPHoP, axis=1)


iphop_melt = iphop.melt(id_vars=['viral_genome'], value_vars=['iPHoP', 'PHIST'])
plot = iphop_melt.groupby(['value', 'variable'], as_index=False).count()
plot['proportion'] = plot['viral_genome']*2/sum(plot['viral_genome'])
ordered = plot.sort_values('value').drop_duplicates('value')
plot['value'] = pd.Categorical(plot['value'], ordered['value'])
plot.rename(columns={"value": "Phylum"}, inplace=True)

# plot host taxonomy and save
fig = px.bar(plot, x='variable', y='proportion',  color='Phylum',
             labels={
                     "variable": "Host taxonomy method",
                     "proportion" : "Proportion of viruses"
                 })

fig.update_layout(title_text='Virus host taxonomic assignment')
fig.write_image(str(snakemake.output))