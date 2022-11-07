import pandas as pd

# load uhgg metadata and blast results
bacteria_db_metadata = pd.read_csv(
# '/home/carsonjm/resources/uhgg/genome-all_metadata.tsv'
str(snakemake.input.bacteria_db_metadata)
, sep='\t')

# split lineage column
bacteria_db_metadata[['superkingdom', 'phylum', 'class', 'order','family', 'genus', 'species']] = bacteria_db_metadata['Lineage'].str.split(';', expand=True)
bacteria_db_metadata['superkingdom'] = bacteria_db_metadata['superkingdom'].str.partition('d__')[2]
bacteria_db_metadata['phylum'] = bacteria_db_metadata['phylum'].str.partition('p__')[2]
bacteria_db_metadata['class'] = bacteria_db_metadata['class'].str.partition('c__')[2]
bacteria_db_metadata['order'] = bacteria_db_metadata['order'].str.partition('o__')[2]
bacteria_db_metadata['family'] = bacteria_db_metadata['family'].str.partition('f__')[2]
bacteria_db_metadata['genus'] = bacteria_db_metadata['genus'].str.partition('g__')[2]
bacteria_db_metadata['superkingdom'] = bacteria_db_metadata['superkingdom'].str.partition('_')[0]
bacteria_db_metadata['phylum'] = bacteria_db_metadata['phylum'].str.partition('_')[0]
bacteria_db_metadata['class'] = bacteria_db_metadata['class'].str.partition('_')[0]
bacteria_db_metadata['order'] = bacteria_db_metadata['order'].str.partition('_')[0]
bacteria_db_metadata['family'] = bacteria_db_metadata['family'].str.partition('_')[0]
bacteria_db_metadata['genus'] = bacteria_db_metadata['genus'].str.partition('_')[0]


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

merged=phist_results_metadata_hq[['phage', 'host', 'common_kmers', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus']]

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
genus_consensus = determine_consensus('genus', phist_counts_taxonomy)
genus_unannotated = phist_counts_taxonomy[~phist_counts_taxonomy['phage'].isin(
    genus_consensus['phage'])]


# determine family level consensus
family_consensus = determine_consensus('family', genus_unannotated)
family_consensus['genus'] = 'NA'
family_unannotated = genus_unannotated[~genus_unannotated['phage'].isin(
    family_consensus['phage'])]

# determine order level consensus
order_consensus = determine_consensus('order', family_unannotated)
order_consensus[['genus', 'family']] = 'NA'
order_unannotated = family_unannotated[~family_unannotated['phage'].isin(
    order_consensus['phage'])]

# determine class level consensus
class_consensus = determine_consensus('class', order_unannotated)
class_consensus[['genus', 'family', 'order']] = 'NA'
class_unannotated = order_unannotated[~order_unannotated['phage'].isin(
    class_consensus['phage'])]

# determine phylum level consensus
phylum_consensus = determine_consensus('phylum', class_unannotated)
phylum_consensus[['genus', 'family', 'order', 'class']] = 'NA'
phylum_unannotated = class_unannotated[~class_unannotated['phage'].isin(
    phylum_consensus['phage'])]

# determine superkingdom level consensus
superkingdom_consensus = determine_consensus('superkingdom', phylum_unannotated)
superkingdom_consensus[['genus', 'family', 'order', 'class', 'phylum']] = 'NA'
superkingdom_unannotated = phylum_unannotated[~phylum_unannotated['phage'].isin(
    superkingdom_consensus['phage'])]

# merge all results
final_consensus = pd.concat([genus_consensus, family_consensus, order_consensus,
    class_consensus, phylum_consensus, superkingdom_consensus])
final_consensus = final_consensus.fillna('NA')
final_consensus['taxonomy'] = (final_consensus['superkingdom'] + ';'
    + final_consensus['phylum'] + ';' + final_consensus['class'] + ';'
    + final_consensus['order'] + ';' + final_consensus['family'] + ';'
    + final_consensus['genus'])

# format results file and save
final_consensus = final_consensus[['phage', 'common_kmers', 'total_kmers',
'superkingdom','superkingdom_kmers','superkingdom_percent_agreement','phylum','phylum_kmers','phylum_percent_agreement',
'class','class_kmers','class_percent_agreement','order','order_kmers','order_percent_agreement',
'family','family_kmers','family_percent_agreement','genus','genus_kmers','genus_percent_agreement', 'taxonomy']]
final_consensus = final_consensus.add_prefix('phist_',)
final_consensus.rename(columns = {'phist_phage': 'viral_genome'}, inplace=True)
final_consensus['viral_genome'] = final_consensus['viral_genome'].str.split('.fna', expand=True)[0]
final_consensus.to_csv(str(snakemake.output.report), index=False)
host_results = final_consensus[['viral_genome', 'phist_taxonomy']]
host_results.to_csv(str(snakemake.output.taxonomy), index=False)