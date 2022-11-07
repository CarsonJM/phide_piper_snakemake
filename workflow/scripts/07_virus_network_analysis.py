# import packages
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import plotly.express as px

# load vcontact2 network
G = nx.read_weighted_edgelist(str(snakemake.input.network))

# draw network
nx.draw(G, node_size=25, node_shape='.', edge_color='dimgray', alpha=0.7)
plt.savefig(str(snakemake.output.network))

# open clustering results
cluster_df = pd.read_csv(str(snakemake.input.clusters))
cluster_group = cluster_df.groupby("Size", as_index=False).count()

# plot kneaddata read counts and save
fig = px.bar(cluster_group, x='Size', y='VC',
             labels={
                     "Size": "Number of viruses in genus cluster",
                     "VC" : "Number of genus clusters"
                 })

fig.update_layout(title_text='Virus dereplication')
fig.update_layout(showlegend=False)
fig.write_html(str(snakemake.output.html))
fig.write_image(str(snakemake.output.svg))