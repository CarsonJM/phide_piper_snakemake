import pandas as pd
import plotly.express as px

# read metaphlan data
lifestyles = pd.read_csv(str(snakemake.input), sep='\t')

# prepare for plotting
lifestyles['Classification'] = lifestyles.apply(lambda x: 'Virulent' if x.Virulent > snakemake.params.bacphlip_prob else 'Unknown', axis=1)
lifestyles['Classification'] = lifestyles.apply(lambda x: 'Temperate' if x.Temperate > snakemake.params.bacphlip_prob else x.Classification, axis=1)
lifestyles_group = lifestyles.groupby("Classification", as_index=False).count()
lifestyles_group['proportions'] = lifestyles_group['Virulent']/sum(lifestyles_group["Virulent"])
lifestyles_group["axis"] = "BACPHLIP"

# plot kneaddata read counts and save
fig = px.bar(lifestyles_group, x='axis', y='proportions', color='Classification',
             labels={
                     "axis": "BACPHLIP Classification",
                     "proportions" : "Proportion of viruses"
                 })

fig.update_xaxes(showticklabels=False) # Hide x axis ticks 

fig.update_layout(title_text='Virus lifestyle classification')
fig.write_html(str(snakemake.output.html))
fig.write_image(str(snakemake.output.svg))