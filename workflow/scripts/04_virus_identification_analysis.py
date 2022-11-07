import pandas as pd
import plotly.express as px
import upsetplot
import matplotlib.pyplot as plt

# read viral identification report
report = pd.read_csv(str(snakemake.input))

tools = []
tool_counts = pd.DataFrame()

# filter based on config paramters
if snakemake.params.run_mgv: 
    mgv_report = report[report['MGV_viral'] == 'Viral']
    tools.append('MGV')
    tool_counts['assembly'] = mgv_report.groupby(['assembly'], as_index=False).count()['assembly']
    tool_counts['MGV'] = mgv_report.groupby(['assembly'], as_index=False).count()['contig_id']

if snakemake.params.run_vf:
    vf_report = report[(report['VirFinder_score'] >= snakemake.params.vf_score)]
    tools.append('VirFinder')
    tool_counts['VirFinder'] = vf_report.groupby(['assembly'], as_index=False).count()['contig_id']

if snakemake.params.run_vs:
    report['VirSorter_cat'] = report.apply(lambda x: x.Category_number if x.Category_text == 'complete_phage' else x.Category_number + 3, axis=1)
    vs_report = report[report['VirSorter_cat'].isin(snakemake.params.vs_cat)]
    tools.append('VirSorter')
    tool_counts['VirSorter'] = vs_report.groupby(['assembly'], as_index=False).count()['contig_id']


if snakemake.params.run_vs2:
    vs2_report = report[report['VirSorter2_max_score'] >= snakemake.params.vs2_score]
    tools.append('VirSorter2')
    tool_counts['VirSorter2'] = vs2_report.groupby(['assembly'], as_index=False).count()['contig_id']

if snakemake.params.run_dvf:
    dvf_report = report[(report['DeepVirFinder_score'] >= snakemake.params.dvf_score)]
    tools.append('DeepVirFinder')
    tool_counts['DeepVirFinder'] = dvf_report.groupby(['assembly'], as_index=False).count()['contig_id']

if snakemake.params.run_vb:
    vb_report = report[report['VIBRANT_viruses'].notnull()]
    tools.append('VIBRANT')
    tool_counts['VIBRANT'] = vb_report.groupby(['assembly'], as_index=False).count()['contig_id']

if snakemake.params.run_genomad:
    genomad_report = report[(report['virus_score'] >= snakemake.params.genomad_score) & (report['fdr'] <= snakemake.params.genomad_fdr)]
    tools.append('geNomad')
    tool_counts['geNomad'] = genomad_report.groupby(['assembly'], as_index=False).count()['contig_id']

if snakemake.params.run_external:
    external_report = report[(report['identity'] >= snakemake.params.min_mash_score) & (report['shared-hashes'] >= snakemake.params.min_mash_hashes) & (report['median-multiplicity'] >= snakemake.params.min_mash_multiplicity)]
    report.rename(columns={'shared-hashes':'shared_hashes', 'median-multiplicity':'median_multiplicity'}, inplace=True)
    tools.append('external')
    tool_counts['external'] = external_report.groupby(['assembly'], as_index=False).count()['query-id']

id = tool_counts.melt(id_vars=['assembly'], value_vars=tools)

# plot kneaddata read counts and save
fig = px.box(id, x='variable', y="value", points="all", color="variable", hover_name='assembly',
             labels={
                     "variable": "Virus identification tool",
                     "value": "Number VLS identified"     
                 })

fig.update_layout(title_text='Virus-like sequence (VLS) identification')
fig.update_layout(showlegend=False)
fig.write_html(str(snakemake.output.boxplot_html))
fig.write_image(str(snakemake.output.boxplot_svg))


# make upset plot to show overlap between tools
report['Tools'] = 'Total'
if snakemake.params.run_mgv: 
    report['Tools'] = report.apply(lambda x: x.Tools + ',MGV' if x.MGV_viral == 'Viral' else x.Tools, axis=1)
if snakemake.params.run_vf:
    report['Tools'] = report.apply(lambda x: x.Tools + ',VirFinder' if x.VirFinder_score >= snakemake.params.vf_score else x.Tools, axis=1)
if snakemake.params.run_vs:
    report['Tools'] = report.apply(lambda x: x.Tools + ',VirSorter' if x.VirSorter_cat in snakemake.params.vs_cat else x.Tools, axis=1)
if snakemake.params.run_vs2:
    report['Tools'] = report.apply(lambda x: x.Tools + ',VirSorter2' if x.VirSorter2_max_score >= snakemake.params.vs2_score else x.Tools, axis=1)
if snakemake.params.run_vb:
    report['Tools'] = report.apply(lambda x: x.Tools + ',VIBRANT' if len(str(x.VIBRANT_viruses)) > 3 else x.Tools, axis=1)
if snakemake.params.run_dvf:
    report['Tools'] = report.apply(lambda x: x.Tools + ',DeepVirFinder' if x.DeepVirFinder_score >= snakemake.params.dvf_score else x.Tools, axis=1)
if snakemake.params.run_genomad:
    report['Tools'] = report.apply(lambda x: x.Tools + ',geNomad' if (x.virus_score >= snakemake.params.min_mash_score and x.fdr <= snakemake.params.genomad_fdr) else x.Tools, axis=1)
if snakemake.params.run_external:
    report['Tools'] = report.apply(lambda x: x.Tools + ',external' if (x.identity >= snakemake.params.min_mash_score and x.shared_hashes >= snakemake.params.min_mash_hashes and x.median_multiplicity >= snakemake.params.min_mash_multiplicity) else x.Tools, axis=1)

output_by_tool = upsetplot.from_memberships(report.Tools.str.split(','), data=report)

fig = plt.gcf()
upsetplot.plot(output_by_tool, sort_by='cardinality', fig=fig)
fig.savefig(str(snakemake.output.upset))
