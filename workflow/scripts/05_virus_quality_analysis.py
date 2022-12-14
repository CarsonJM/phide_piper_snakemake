import pandas as pd
import plotly.express as px

checkv = pd.read_csv(
    str(snakemake.input)
    # '/home/carsonjm/CarsonJM/phide_piper/results/06_VIRUS_QUALITY/01_checkv/quality_summary.tsv'
    , sep='\t')

checkv_group = checkv.groupby('checkv_quality', as_index=False).count()
checkv_group['checkv_quality'] = pd.Categorical(checkv_group['checkv_quality'], ['Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete'])

# plot kneaddata read counts and save
fig = px.bar(checkv_group, x='checkv_quality', y="contig_id", color="checkv_quality",
             category_orders={'checkv_quality':['Not-determined', 'Low-quality', 'Medium-quality', 'High-quality', 'Complete']},
             labels={
                     "checkv_quality": "CheckV quality category",
                     "contig_id": "Number of VLS"
                 })

fig.update_layout(title_text='Virus-like sequence (VLS) quality assessment')
fig.update_layout(showlegend=False)
fig.write_image(str(snakemake.output))