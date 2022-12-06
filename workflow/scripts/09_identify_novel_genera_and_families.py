import pandas as pd
from Bio import SeqIO

# load mash distances
mash_dist = pd.read_csv(str(snakemake.input.mash_dist), sep='\t', names=['viral_genome', 'accession_full', 'distance', 'p_value', 'shared_hashes'])

# merge mash hits and metadata
mash_sort = mash_dist.sort_values('distance')
mash_group = mash_sort.groupby('viral_genome', as_index=False).first()

# identify highly similar sequences to extract for species-level analysis
mash_hi_sim = mash_group[mash_group['distance'] <= snakemake.params.max_dist_species]
hi_sim = set(mash_hi_sim['accession_full'])

viruses_for_species_analysis = []

for record in SeqIO.parse(str(snakemake.input.virus_db), 'fasta'):
    if record.id in hi_sim:
        viruses_for_species_analysis.append(record)
    else:
        continue

SeqIO.write(viruses_for_species_analysis, str(snakemake.output.hi_sim_viruses), 'fasta')


# load virus metadata and merge with mash
metadata = pd.read_csv(str(snakemake.input.meta))
mash_group['accession'] = mash_group['accession_full'].str.split('.', expand=True)[0]
mash_meta = mash_group.merge(metadata, on='accession', how='left')
# filter to identify genus-level hits
mash_meta['novel_genus'] = mash_meta.apply(lambda x: False if x.distance <= snakemake.params.max_dist_genus else True, axis=1)
mash_meta['novel_family'] = mash_meta.apply(lambda x: False if x.distance <= snakemake.params.max_dist_family else True, axis=1)
# save top hits for each genome
mash_meta.to_csv(str(snakemake.output.mash_results), index=False)