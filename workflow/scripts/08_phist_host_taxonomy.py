import pandas as pd

# load uhgg metadata and blast results
bacteria_db_metadata = pd.read_csv(
# '/home/carsonjm/resources/uhgg/genome-all_metadata.tsv'
str(snakemake.input.bacteria_db_metadata))

# read phist report line-by-line
phist = open(str(snakemake.input.phist), 'r')

Lines = phist.readlines()

# create lists to store data in
phage_col = []
host_col = []
kmer_col = []

# if there are multiple common kmers, write out phage, host, and # kmers
for line in Lines:
    if "kmer-length:" in line:
        viruses = line.split(',')
        viruses.pop(0)
        viruses.pop(0)
        viruses.pop()
    elif "query-samples" in line:
        continue
    elif ':' in line:
        bacterial_genome = line.split(',')[0]
        hits = line.split(',')[1:]
        hits.pop(0)
        hits.pop()
        for hit in hits:
            index, kmers = str(hit).split(':')
            if int(kmers) > 1:
                phage_col.append(str(viruses[int(index) - 1]))
                host_col.append(str(bacterial_genome))
                kmer_col.append(int(kmers))
            else:
                continue

phist.close()

phist_results = pd.DataFrame()
phist_results['phage'] = phage_col
phist_results['host'] = host_col
phist_results["common_kmers"] = kmer_col

# format blast genome so it can be merged with metadata, then merge
phist_results['Genome'] = phist_results['host'].str.split('.', expand=True)[0]
phist_results_metadata = phist_results.merge(bacteria_db_metadata, on='Genome')

# filter to only retain genomes with significant number of common kmers
phist_results_metadata_hq = phist_results_metadata[phist_results_metadata['common_kmers'] >=
snakemake.params.min_common_kmers
# 0.5
]

merged=phist_results_metadata_hq[['phage', 'host', 'common_kmers', 'Genus']]

merged2 = merged[merged['host'].notnull()]

# count total common kmers
phist_count = merged2.groupby(by=['phage'], as_index=False).sum()
phist_counts = phist_count.loc[:,['phage', 'common_kmers']]
phist_counts.rename(columns = {'common_kmers':'total_kmers'}, inplace = True)
phist_counts_taxonomy = merged2.merge(phist_counts, on='phage', how='left')

# determine if any genomes have > min_agreement agreement at taxonomic level
def determine_consensus(taxonomic_rank, phist_table):
    taxonomy_hits = phist_table.groupby(['phage', taxonomic_rank], as_index=False
            ).agg(phist_taxonomy_hits=('common_kmers' ,'sum'))
    taxonomy_hits.rename(columns={'phist_taxonomy_hits':taxonomic_rank + '_kmers'}, inplace=True)
    taxonomy_hits_merged = phist_table.merge(taxonomy_hits, on=['phage', taxonomic_rank], how='outer')
    taxonomy_hits_merged2 = taxonomy_hits_merged[taxonomy_hits_merged[taxonomic_rank + "_kmers"].notnull()]
    taxonomy_hits_merged2[taxonomic_rank + '_percent_agreement'] = taxonomy_hits_merged2[taxonomic_rank + '_kmers'].astype(int)/taxonomy_hits_merged2['total_kmers']
    taxonomy_hits_consensus = taxonomy_hits_merged2[taxonomy_hits_merged2[taxonomic_rank + '_percent_agreement']*100 >
    # 70
    snakemake.params.min_agreement
    ]
    taxonomy_hits_consensus = taxonomy_hits_consensus.groupby('phage', as_index=False).first()

    return taxonomy_hits_consensus

# determine genus level consensus
genus_consensus = determine_consensus('Genus', phist_counts_taxonomy)
genus_unannotated = phist_counts_taxonomy[~phist_counts_taxonomy['phage'].isin(
    genus_consensus['phage'])]


# format results file and save
final_consensus = genus_consensus[['phage', 'common_kmers', 'total_kmers','Genus','Genus_kmers','Genus_percent_agreement']]
final_consensus = final_consensus.add_prefix('phist_',)
final_consensus.rename(columns = {'phist_phage': 'viral_genome'}, inplace=True)
final_consensus['viral_genome'] = final_consensus['viral_genome'].str.split('.fna', expand=True)[0]
final_consensus.to_csv(str(snakemake.output.report), index=False)
host_results = final_consensus[['viral_genome', 'phist_Genus']]
host_results.to_csv(str(snakemake.output.taxonomy), index=False)