import pandas as pd

### Section for all analyses ###
# read samples.tsv
samples = pd.read_csv(str(snakemake.input.samples), sep='\t')
samples.loc[:,'sample'] = samples.loc[:,"group"] + "_" + samples.loc[:,"sample"]

# read dereplication report
dereplication = pd.read_csv(str(snakemake.input.dereplication), sep='\t')
dereplication.rename(columns={'representative':'viral_genome', 'cluster_size':'replicates'}, inplace=True)
dereplication['sample'] = dereplication['viral_genome'].str.split('|', expand=True)[0]

# remove checkv provirus suffix to merge with genomad results
dereplication["viral_genome_base"] = dereplication.apply(lambda x: x.viral_genome.split("|checkv_provirus")[0].rpartition('_')[0] if "|checkv_provirus" in str(x.viral_genome) else x.viral_genome, axis=1)
summary = samples.merge(dereplication, on='sample', how='outer')[['sample', 'viral_genome', 'viral_genome_base', 'replicates']]

# read votu diversity report
clusters = open(str(snakemake.input.diversity), 'r')

genome_w_rep = {}
rep_w_votu_size = {}

for line in clusters:
    stripped = line.strip()
    centroid, nodes = stripped.split('\t')
    nodes_list = nodes.split(",")
    for node in nodes_list:
        genome_w_rep[node] = centroid
    rep_w_votu_size[centroid] = len(nodes_list)

votu_size = pd.DataFrame(list(rep_w_votu_size.items()), columns=['votu_representative', 'votu_size'])
votu_rep = pd.DataFrame(list(genome_w_rep.items()), columns=['viral_genome', 'votu_representative'])
diversity = votu_rep.merge(votu_size, on='votu_representative')
summary = summary.merge(diversity, on='viral_genome', how='left')


### Section for read analyses ###
# read count info
if snakemake.params.input_data == 'reads':
    read_count = pd.read_csv(str(snakemake.input.read_count), sep="\t")
    read_count_summary = read_count[['Sample', 'final pair1']]
    read_count_summary.rename(columns={'Sample':'sample', 'final pair1': 'decontaminated_reads'}, inplace=True)
    summary = summary.merge(read_count_summary, on="sample", how='left')

# enrichment info
if snakemake.params.include_enrichment:
    enrichment = pd.read_csv(str(snakemake.input.enrichment), sep="\t")
    enrichment_summary = enrichment[['sample', 'total enrichmnet score']]
    enrichment_summary.rename(columns={'total enrichmnet score': 'total_enrichment_score'}, inplace=True)
    summary = summary.merge(enrichment_summary, on="sample", how='left')

if snakemake.params.input_data == 'reads':
# assembly info
    assembly = pd.read_csv(str(snakemake.input.assembly), sep='\t')
    assembly_summary =assembly[['Assembly', '# contigs (>= 1000 bp)', 'Total length (>= 1000 bp)', 'GC (%)', 'N50', 'L50', "# N's per 100 kbp"]]
    assembly_summary.rename(columns={"Assembly":"sample", "# contigs (>= 1000 bp)":"num_contigs", "Total length (>= 1000 bp)":"total_length", "GC (%)":"GC_perc", "# N's per 100 kbp":"Ns_per_100kbp"}, inplace=True)
    summary = summary.merge(assembly_summary, on='sample', how='left')

# identification info
    identification = pd.read_csv(str(snakemake.input.identification), sep='\t')
    identification_summary = identification[['vls_id', 'length', 'topology', 'coordinates', 'n_genes', 'genetic_code', 'virus_score', 'n_hallmarks', 'identity', 'shared-hashes', 'median-multiplicity', 'assembly']]
    identification_summary.rename(columns={"vls_id":"viral_genome_base", "length":"genomad_length", "topology":"genomad_topology", "coordinates":"genomad_coordinates", "n_genes":"genomad_gene_count", "n_hallmakrs":"genomad_hallmark_count", "shared-hashes":"shared_hashes", "median-multiplicity":"median_multiplicity", "assembly":"sample_w_external_hit"}, inplace=True)
    identification_summary.sort_values('identity', ascending=False, inplace=True)
    identification_summary_top = identification_summary.groupby('viral_genome_base', as_index=False).first()
    identification_summary_top["sample_w_external_hit"] = identification_summary_top.apply(lambda x: x.sample_w_external_hit if "external|" in str(x.viral_genome_base) else '', axis=1)
    summary = summary.merge(identification_summary_top, on='viral_genome_base', how='left')

# checkv info
    quality = pd.read_csv(str(snakemake.input.quality), sep='\t')
    quality_summary = quality[['contig_id', 'provirus', 'proviral_length', 'gene_count', 'viral_genes', 'host_genes', 'checkv_quality', 'completeness', 'completeness_method', 'kmer_freq', 'warnings']]
    quality_summary.rename(columns={'contig_id':'viral_genome_base', 'provirus':'checkv_provirus', 'proviral_length':'checkv_proviral_length', 'gene_count':'checkv_gene_count', "viral_genes":"checkv_viral_genes", "host_genes":"checkv_host_genes"}, inplace=True)
    quality_summary_first = quality_summary.groupby('viral_genome_base', as_index=False).first()
    summary = summary.merge(quality_summary_first, on='viral_genome_base', how='left')

### Section for contig analyses ###
if snakemake.params.input_data == 'contigs':
# genomad info
    identification = pd.read_csv(str(snakemake.input.identification), sep='\t')
    identification_summary = identification[['vls_id', 'length', 'topology', 'coordinates', 'n_genes', 'genetic_code', 'virus_score', 'n_hallmarks', 'identity', 'shared-hashes', 'median-multiplicity', 'assembly']]
    identification_summary.rename(columns={"vls_id":"viral_genome_base", "length":"genomad_length", "topology":"genomad_topology", "coordinates":"genomad_coordinates", "n_genes":"genomad_gene_count", "n_hallmakrs":"genomad_hallmark_count", "shared-hashes":"shared_hashes", "median-multiplicity":"median_multiplicity", "assembly":"sample_w_external_hit"}, inplace=True)
    identification_summary.sort_values('identity', ascending=False, inplace=True)
    identification_summary_top = identification_summary.groupby('viral_genome_base', as_index=False).first()
    identification_summary_top["sample_w_external_hit"] = identification_summary_top.apply(lambda x: x.sample_w_external_hit if "external|" in str(x.viral_genome_base) else '', axis=1)
    summary = summary.merge(identification_summary_top, on='viral_genome_base', how='left')

# checkv info
    quality = pd.read_csv(str(snakemake.input.quality), sep='\t')
    quality_summary = quality[['contig_id', 'provirus', 'proviral_length', 'gene_count', 'viral_genes', 'host_genes', 'checkv_quality', 'completeness', 'completeness_method', 'kmer_freq', 'warnings']]
    quality_summary.rename(columns={'contig_id':'viral_genome_base', 'provirus':'checkv_provirus', 'proviral_length':'checkv_proviral_length', 'gene_count':'checkv_gene_count', "viral_genes":"checkv_viral_genes", "host_genes":"checkv_host_genes"}, inplace=True)
    quality_summary_first = quality_summary.groupby('viral_genome_base', as_index=False).first()
    summary = summary.merge(quality_summary_first, on='viral_genome_base', how='left')

### Section for vls analyses ###
if snakemake.params.input_data == 'vls':
# checkv info
    quality = pd.read_csv(str(snakemake.input.quality), sep='\t')
    quality_summary = quality[['contig_id', 'provirus', 'proviral_length', 'gene_count', 'viral_genes', 'host_genes', 'checkv_quality', 'completeness', 'completeness_method', 'kmer_freq', 'warnings']]
    quality_summary.rename(columns={'contig_id':'viral_genome_base', 'provirus':'checkv_provirus', 'proviral_length':'checkv_proviral_length', 'gene_count':'checkv_gene_count', "viral_genes":"checkv_viral_genes", "host_genes":"checkv_host_genes"}, inplace=True)
    summary = summary.merge(quality_summary, on='viral_genome_base', how='left')

### Sections that are individually specified ###
# phage taxonomy info
if snakemake.params.include_host:
    taxonomy = pd.read_csv(str(snakemake.input.taxonomy), sep='\t')
    taxonomy_summary = taxonomy[['seq_name', 'n_genes_with_taxonomy', 'agreement', 'lineage']]
    taxonomy_summary.rename(columns={'seq_name':'viral_genome', 'agreement':'genomad_agreement', 'lineage':'genomad_taxonomy'}, inplace=True)
    summary = summary.merge(taxonomy_summary, on='viral_genome', how='left')

# phage taxonomy info
if snakemake.params.include_taxonomy:
    taxonomy = pd.read_csv(str(snakemake.input.taxonomy), sep='\t')
    taxonomy_summary = taxonomy[['seq_name', 'n_genes_with_taxonomy', 'agreement', 'lineage']]
    taxonomy_summary.rename(columns={'seq_name':'viral_genome', 'agreement':'genomad_agreement', 'lineage':'genomad_taxonomy'}, inplace=True)
    summary = summary.merge(taxonomy_summary, on='viral_genome', how='left')

# phage lifestyle info
if snakemake.params.include_lifestyle:
    lifestyle = pd.read_csv(str(snakemake.input.lifestyle), sep='\t')
    lifestyle.rename(columns={'Virulent':'BACPHLIP_virulent_score', 'Temperate':'BACPHLIP_temperate_score', 'Classification':'lifestyle_prediction'}, inplace=True)
    summary = summary.merge(lifestyle, on='viral_genome', how='left')

# phage abundance/activity info
if snakemake.params.include_analysis:
    analysis = pd.read_csv(str(snakemake.input.analysis), sep='\t')
    analysis_summary = analysis[["genome", "sample", "relative_abundance"]]
    analysis_summary.rename(columns={"genome":"viral_genome"}, inplace=True)
    analysis_summary_pivot = analysis_summary.pivot('viral_genome', 'sample', 'relative_abundance')
    summary = summary.merge(analysis_summary_pivot, on='viral_genome', how='left')

summary.to_csv(str(snakemake.output), sep='\t', index=False)