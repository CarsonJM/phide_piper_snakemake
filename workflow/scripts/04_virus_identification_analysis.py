import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt

# read viral identification report
report = pd.read_csv(str(snakemake.input), sep='\t')

tools = []
tool_counts = pd.DataFrame()

genomad_report = report[(report['virus_score'] >= snakemake.params.genomad_score)]
tools.append('geNomad')
tool_counts[['assembly', 'geNomad']] = genomad_report.groupby(['assembly'], as_index=False).count()[['assembly', 'vls_id']]

external_report = report[(report['identity'] >= snakemake.params.min_mash_score) & (report['shared-hashes'] >= snakemake.params.min_mash_hashes) & (report['median-multiplicity'] >= snakemake.params.min_mash_multiplicity)]
external_report.rename(columns={'shared-hashes':'shared_hashes', 'median-multiplicity':'median_multiplicity'}, inplace=True)
tools.append('external')
tool_counts[['assembly','external']] = external_report.groupby(['assembly'], as_index=False).count()[['assembly', 'vls_id']]

id = tool_counts.melt(id_vars=['assembly'], value_vars=tools)

# plot kneaddata read counts and save
fig = px.box(id, x='variable', y="value", points="all", color="variable", hover_name='assembly',
             labels={
                     "variable": "Virus identification tool",
                     "value": "Number VLS identified"
                 })

fig.update_layout(title_text='Virus-like sequence (VLS) identification')
fig.update_layout(showlegend=False)
fig.write_image(str(snakemake.output))
