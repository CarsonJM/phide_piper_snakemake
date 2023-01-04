import pandas as pd
import plotly.express as px

# read dereplicate rep data
derep_reps = pd.read_csv(str(snakemake.input.derep_reps), sep='\t', names=['viral_genome', 'nodes'])

# read metaphlan data
bacphlip = pd.read_csv(str(snakemake.input.bacphlip), sep='\t', names=['viral_genome', 'Virulent', 'Temperate'], header=0)
lifestyles = derep_reps.merge(bacphlip, on="viral_genome", how='outer')
lifestyles = lifestyles[['viral_genome', 'Virulent', 'Temperate']]


# prepare for plotting
lifestyles['Classification'] = lifestyles.apply(lambda x: 'Virulent' if float(x.Virulent) > snakemake.params.bacphlip_prob else 'Unknown', axis=1)
lifestyles['Classification'] = lifestyles.apply(lambda x: 'Temperate' if float(x.Temperate) > snakemake.params.bacphlip_prob else x.Classification, axis=1)
lifestyles['Classification'] = lifestyles.apply(lambda x: 'Temperate' if 'provirus' in x.viral_genome else x.Classification, axis=1)
lifestyles['Integration_status'] = lifestyles.apply(lambda x: 'integrated_provirus' if 'provirus' in x.viral_genome else 'unintegrated_virus', axis=1)

lifestyles.to_csv(str(snakemake.output.report), sep='\t', index=False)

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
fig.write_image(str(snakemake.output.svg))