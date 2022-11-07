import pandas as pd
import plotly.express as px

# read metaphlan data
lifestyles = pd.read_csv(str(snakemake.input))

# prepare for plotting
bacphlip_count = lifestyles[['viral_genome', 'BACPHLIP_class']]
bacphlip_count['tool'] = 'BACPHLIP'
bacphlip_count['Lifestyle classification'] = bacphlip_count['BACPHLIP_class']
bacphlip_count2 = bacphlip_count[['Lifestyle classification', 'tool']]
bacphlip_group = bacphlip_count2.groupby('Lifestyle classification', as_index=False).count()
bacphlip_group['axis'] = "BACPHLIP"
bacphlip_group['proportions'] = bacphlip_group['tool']/sum(bacphlip_group['tool'])

# plot kneaddata read counts and save
fig = px.bar(bacphlip_group, x='axis', y='proportions', color='Lifestyle classification',
             labels={
                     "axis": "BACPHLIP Classification",
                     "proportions" : "Proportion of viruses"
                 })

fig.update_xaxes(showticklabels=False) # Hide x axis ticks 

fig.update_layout(title_text='Virus lifestyle classification')
fig.write_html(str(snakemake.output.html))
fig.write_image(str(snakemake.output.svg))