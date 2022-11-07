import pandas as pd
import plotly.express as px

checkv = pd.read_csv(
    str(snakemake.input)
    # '/home/carsonjm/CarsonJM/phide_piper/results/06_VIRUS_QUALITY/01_checkv/quality_summary.tsv'
    , sep='\t')

checkv_group = checkv.groupby('checkv_quality', as_index=False).count()

# plot kneaddata read counts and save
fig = px.bar(checkv_group, x='checkv_quality', y="contig_id", color="checkv_quality",
             labels={
                     "checkv_quality": "CheckV quality category",
                     "contig_id": "Number of VLS"
                 })

fig.update_layout(title_text='Virus-like sequence (VLS) quality assessment')
fig.update_layout(showlegend=False)
fig.write_html(str(snakemake.output.html))
fig.write_image(str(snakemake.output.svg))