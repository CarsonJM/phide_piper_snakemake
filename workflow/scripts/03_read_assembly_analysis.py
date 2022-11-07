import pandas as pd
from plotnine import *

assembly = pd.read_csv(str(snakemake.input), sep='\t')
assembly_count = assembly.copy()
assembly_len = assembly.copy()

assembly_count.rename(columns = {'# contigs (>= 0 bp)':'0',
'# contigs (>= 1000 bp)':'1000', 
'# contigs (>= 5000 bp)':'5000', 
'# contigs (>= 10000 bp)':'10000'},
inplace = True)

assembly_len.rename(columns = {'Total length (>= 0 bp)':'0', 
'Total length (>= 1000 bp)':'1000', 
'Total length (>= 5000 bp)':'5000', 
'Total length (>= 10000 bp)':'10000'},
inplace = True)

assembly_count_melt = assembly_count.melt(id_vars=['Assembly'], value_vars=['0', '1000', '5000', '10000'])
assembly_len_melt = assembly_len.melt(id_vars=['Assembly'], value_vars=['0', '1000', '5000', '10000'])
assembly_count_melt['variable'] = pd.Categorical(assembly_count_melt['variable'], categories=['0', '1000', '5000', '10000'], ordered = True)
assembly_len_melt['variable'] = pd.Categorical(assembly_len_melt['variable'], categories=['0', '1000', '5000', '10000'], ordered = True)
assembly_count_melt['type'] = "Contig count"
assembly_len_melt['type'] = "Combined contig length"

assembly_concat = pd.concat([assembly_count_melt, assembly_len_melt], axis=0)

assembly_plot = (
    ggplot(assembly_concat)
    + geom_boxplot(aes(x='variable', y='value'))
    + xlab("Minimum contig length")
    + ylab("")
    + facet_wrap('type', dir='v', ncol=1, scales='free')
)

assembly_plot.save(str(snakemake.output), dpi=600)