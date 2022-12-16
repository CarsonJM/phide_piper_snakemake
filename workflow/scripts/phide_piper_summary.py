import pandas as pd

# read samples.tsv
samples = pd.read_csv(str(snakemake.input.samples), sep='\t')

# read votu diversity report
diversity = pd.read_csv(str(snakemake.input.diversity), sep='\t')
diversity.rename(columns={'representative':'viral_genome', 'cluster_size':'votu_members'}, inplace=True)
diversity['sample'] = diversity['viral_genome'].str.split('|', expand=True)[0]
summary = samples.merge(diversity, on='sample', how='left')[['sample', 'viral_genome', 'votu_members']]

dereplication = pd.read_csv(str(snakemake.input.dereplication), sep='\t')
dereplication.rename(columns={'representative':'viral_genome', 'cluster_size':'replicates'}, inplace=True)
summary = summary.merge(dereplication, on='viral_genome', how='left')
summary = summary[['viral_genome', 'sample', 'votu_members', 'replicates']]

# enrichment info
if snakemake.params.include_enrichment:
    enrichment = pd.read_csv(str(snakemake.input.enrichment), sep="\t")
    enrichment_summary = enrichment[['sample', 'total enrichmnet score']]
    enrichment_summary.rename(columns={'total enrichmnet score': 'total_enrichment_score'}, inplace=True)
    summary = summary.merge(enrichment_summary, on="sample", how='left')

if snakemake.params.input_data == 'reads':
    assembly = pd.read_csv(str(snakemake.input.assembly), sep='\t')
    assembly_summary =assembly[['Assembly', '# contigs (>= 1000 bp)', 'Total length (>= 1000 bp)', 'GC (%)', 'N50', 'L50', "# N's per 100 kbp"]]
    assembly_summary.rename(columns={"Assembly":"sample", "# contigs (>= 1000 bp)":"num_contigs", "Total length (>= 1000 bp)":"total_length", "GC (%)":"GC_perc", "# N's per 100 kbp":"Ns_per_100kbp"}, inplace=True)
    summary = summary.merge(assembly_summary, on='sample', how='left')

    identification = pd.read_csv(str(snakemake.input.identification), sep='\t')
    identification_summary = identification[['vls_id', 'length', 'topology', 'coordinates', 'n_genes', 'genetic_code', 'virus_score', 'n_hallmarks', 'identity', 'shared-hashes', 'median-multiplicity']]
    identification_summary.rename(columns={"vls_id":"viral_genome", "length":"genomad_length", "topology":"genomad_topology", "coordinates":"genomad_coordinates", "n_genes":"genomad_gene_count", "n_hallmakrs":"genomad_hallmark_count", "shared-hashes":"shared_hashes", "median-multiplicity":"median_multiplicity"}, inplace=True)
    summary = summary.merge(identification_summary, on='viral_genome', how='left')

    quality = pd.read_csv(str(snakemake.input.quality), sep='\t')
    quality_summary = quality[['contig_id', 'provirus', 'proviral_length', 'gene_count', 'viral_genes', 'host_genes', 'checkv_quality', 'completeness', 'completeness_method', 'kmer_freq', 'warnings']]
    quality_summary.rename(columns={'contig_id':'viral_genome', 'provirus':'checkv_provirus', 'proviral_length':'checkv_proviral_length', 'gene_count':'checkv_gene_count', "viral_genes":"checkv_viral_genes", "host_genes":"checkv_host_genes"}, inplace=True)
    summary = summary.merge(quality_summary, on='viral_genome', how='left')

if snakemake.params.input_data == 'contigs':
    identification = pd.read_csv(str(snakemake.input.identification), sep='\t')
    identification_summary = identification[['vls_id', 'length', 'topology', 'coordinates', 'n_genes', 'genetic_code', 'virus_score', 'n_hallmarks', 'identity', 'shared-hashes', 'median-multiplicity']]
    identification_summary.rename(columns={"vls_id":"viral_genome", "length":"genomad_length", "topology":"genomad_topology", "coordinates":"genomad_coordinates", "n_genes":"genomad_gene_count", "n_hallmakrs":"genomad_hallmark_count", "shared-hashes":"shared_hashes", "median-multiplicity":"median_multiplicity"}, inplace=True)
    summary = summary.merge(identification_summary, on='viral_genome', how='left')

    quality = pd.read_csv(str(snakemake.input.quality), sep='\t')
    quality_summary = quality[['contig_id', 'provirus', 'proviral_length', 'gene_count', 'viral_genes', 'host_genes', 'checkv_quality', 'completeness', 'completeness_method', 'kmer_freq', 'warnings']]
    quality_summary.rename(columns={'contig_id':'viral_genome', 'provirus':'checkv_provirus', 'proviral_length':'checkv_proviral_length', 'gene_count':'checkv_gene_count', "viral_genes":"checkv_viral_genes", "host_genes":"checkv_host_genes"}, inplace=True)
    summary = summary.merge(quality_summary, on='viral_genome', how='left')


if snakemake.params.input_data == 'vls':
    quality = pd.read_csv(str(snakemake.input.quality), sep='\t')
    quality_summary = quality[['contig_id', 'provirus', 'proviral_length', 'gene_count', 'viral_genes', 'host_genes', 'checkv_quality', 'completeness', 'completeness_method', 'kmer_freq', 'warnings']]
    quality_summary.rename(columns={'contig_id':'viral_genome', 'provirus':'checkv_provirus', 'proviral_length':'checkv_proviral_length', 'gene_count':'checkv_gene_count', "viral_genes":"checkv_viral_genes", "host_genes":"checkv_host_genes"}, inplace=True)
    summary = summary.merge(quality_summary, on='viral_genome', how='left')

summary.to_csv(str(snakemake.output), sep='\t', index=False)