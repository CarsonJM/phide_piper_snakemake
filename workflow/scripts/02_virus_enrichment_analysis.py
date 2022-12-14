# import packages
import plotly.express as px
import pandas as pd

# read kneaddata read counts
vqc = pd.read_csv(str(snakemake.input), sep='\t')
vqc.rename(columns={'total enrichmnet score':'total enrichment score'}, inplace=True)
vqc_melt = vqc.melt(id_vars=['Sample'], value_vars=['total enrichment score'])

# plot kneaddata read counts and save
fig = px.box(vqc_melt, y="value", points="all", hover_name='Sample',
             labels={
                     "value": "Fold enrichment relative to bulk"
                 })
fig.update_layout(title_text='Virome enrichment')
fig.update_layout(showlegend=False)
fig.write_image(str(snakemake.output))