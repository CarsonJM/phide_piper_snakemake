import pandas as pd
import numpy as np
import plotly.express as px

# load virus abundance information
# virus_abund = pd.read_csv(str(snakemake.input.virus_abund), sep='\t')
virus_abund = pd.read_csv('/gscratch/scrubbed/carsonjm/results/phide_piper_test/12_VIRUS_ABUNDANCE/02_instrain_profile/combined_genome_info_w_abundance.tsv', sep='\t')
virus_abund.rename(columns={'genome':'viral_genome'}, inplace=True)
virus_abund_summary = virus_abund[["viral_genome", "sample", "relative_abundance"]]
virus_abund_summary['relative_abundance'] = virus_abund_summary['relative_abundance'] * 100

# load virus host information
# virus_host = pd.read_csv(str(snakemake.input.virus_host), sep='\t')
virus_host = pd.read_csv('/gscratch/scrubbed/carsonjm/results/phide_piper_test/08_VIRUS_HOST/virus_host_taxonomy_report.tsv', sep='\t')

# merge virus abundance with host info
virus_abund_host = virus_abund_summary.merge(virus_host, on='viral_genome', how='outer')
virus_abund_host['genus'] = virus_abund_host.loc[:, 'iphop_genus'].str.split(';g__', expand=True)[1]
virus_abund_host['viral_genome_genus'] = virus_abund_host.loc[:, 'viral_genome'].astype(str) + ';genus_' + virus_abund_host.loc[:, 'genus'].astype(str)
virus_abund_host_summary = virus_abund_host[['viral_genome_genus', 'sample', 'relative_abundance']]
virus_abund_host_grp = virus_abund_host_summary.groupby(['viral_genome_genus'], as_index=False).count()
virus_abund_host_persistent = virus_abund_host_summary[virus_abund_host_summary['viral_genome_genus'].isin(virus_abund_host_grp[virus_abund_host_grp['sample'] >= 3]['viral_genome_genus'])]
virus_abund_pivot = virus_abund_host_persistent.pivot('sample', 'viral_genome_genus', 'relative_abundance')
virus_abund_pivot.fillna(0, inplace=True)

# load host abundance information
# host_abund = pd.read_csv(str(snakemake.input.host_abun), sep='\t', header=1)
host_abund = pd.read_csv('/gscratch/scrubbed/carsonjm/results/phide_piper_test/13_HOST_ABUNDANCE/metaphlan_merged_profiles_gtdb.tsv', sep='\t', header=1)

# split host column into genera
host_abund_genus = host_abund[(host_abund["clade_name"].str.contains('s__') == True)]
host_abund_genus['genus'] = host_abund_genus['clade_name'].str.split(';g__', expand=True)[1]
host_abund_genus_summary = host_abund_genus.drop(['clade_name'], axis=1)

# format host abundance data for correlation analysis
host_abund_genus_melt = host_abund_genus_summary.melt(id_vars='genus')
host_abund_genus_melt.rename(columns={"variable":"sample"}, inplace=True)
host_abund_genus_melt.loc[:,'sample'] = host_abund_genus_melt.loc[:,'sample'].str.split('_gtdb', expand=True)[0]
host_abund_genus_melt.replace(0, np.nan, inplace=True)
host_abund_genus_grp = host_abund_genus_melt.groupby(['genus'], as_index=False).count()
host_abund_genus_persistent = host_abund_genus_melt[host_abund_genus_melt['genus'].isin(host_abund_genus_grp[host_abund_genus_grp['value'] >= 3]['genus'])]
host_abund_pivot = host_abund_genus_persistent.pivot('sample', 'genus', 'value')
host_abund_pivot.fillna(0, inplace=True)

# merge virus info with host abund
virus_host_merge = host_abund_pivot.merge(virus_abund_pivot, on='sample', how='outer')

# calculate pearson correlation
corr = virus_host_merge.corr()
corr.to_csv('test_corr.tsv', sep='\t')
spearman = virus_host_merge.corr('spearman')
spearman.to_csv('test_spearman.tsv', sep='\t')

virus_host_merge['sample'] = virus_host_merge.index
virus_host_merge_melt = virus_host_merge.melt(id_vars='sample')

# plot rothia phage
rothia_phage = virus_host_merge_melt[(virus_host_merge_melt['variable'].str.contains('Rothia'))]
fig = px.line(rothia_phage, x='sample', y='value', color='variable', markers=True)
fig.update_layout(legend=dict(
    yanchor="top",
    y=1.5,
    xanchor="left",
    x=0
))
fig.write_image('rothia_phage.png')

mogi_div_phage = virus_host_merge_melt[(virus_host_merge_melt['variable'].str.contains('Mogibacterium;s__Mogibacterium diversum')) | (virus_host_merge_melt['variable'].isin(spearman[spearman['Mogibacterium;s__Mogibacterium diversum'] <= -0.9].index) & ((virus_host_merge_melt['variable'].str.contains(';genus_'))))]
fig = px.line(mogi_div_phage, x='sample', y='value', color='variable', markers=True)
fig.update_layout(legend=dict(
    yanchor="top",
    y=1.5,
    xanchor="left",
    x=0
))
fig.write_image('mogi_div_phage.png')