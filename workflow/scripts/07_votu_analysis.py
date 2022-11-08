import pandas as pd
import plotly.express as px

# open clustering results
clusters = open(
    str(snakemake.input)
    # '/home/carsonjm/CarsonJM/phide_piper/results/05_VIRUS_CLUSTERING/03_cluster_viruses/virus_virusdb_combined_clusters.tsv'
    , 'r')

node_lengths = []
for line in clusters:
    stripped = line.strip()
    centroid, nodes = stripped.split('\t')
    nodes_split = nodes.split(",")
    node_lengths.append(len(nodes_split))

cluster_df = pd.DataFrame()
cluster_df['node_lengths'] = node_lengths
cluster_df["count"] = True
cluster_group = cluster_df.groupby("node_lengths", as_index=False).count()

# plot kneaddata read counts and save
fig = px.bar(cluster_group, x='node_lengths', y='count',
             labels={
                     "node_lengths": "Number of viruses vOTU",
                     "count" : "Number of vOTUs"
                 })

fig.update_layout(title_text='Virus dereplication')
fig.update_layout(showlegend=False)
fig.write_html(str(snakemake.output.html))
fig.write_image(str(snakemake.output.svg))