import pandas as pd
from Bio import SeqIO

# combine instrain results within groups
group_df = pd.DataFrame()

for file in snakemake.input.profiles:
    df = pd.read_csv(file, sep='\t')
    df_summary = df[['genome', 'sample', 'relative_abundance']]
    df_summary.rename(columns={"genome": "votu_representative"}, inplace=True)
    group_df = group_df.append(df_summary)

# sort and retain 1 row for each genome present in at least 1 sample in the group
group_df_sort = group_df.sort_values('relative_abundance', ascending=False)
group_df_top = group_df_sort.groupby('votu_representative', as_index=False).first()

# load votu membership file
clusters = open(str(snakemake.input.votu_members), 'r')

genome_w_rep = {}

for line in clusters:
    stripped = line.strip()
    centroid, nodes = stripped.split('\t')
    nodes_list = nodes.split(",")
    for node in nodes_list:
        genome_w_rep[node] = centroid

votu_rep = pd.DataFrame(list(genome_w_rep.items()), columns=['viral_genome', 'votu_representative'])

# merge abundance info with votu membership
group_df_votu = group_df_top.merge(votu_rep, on='votu_representative', how='left')

# identify untrimmed prophages from a sample in the group
group_df_provirus = group_df_votu[(group_df_votu['viral_genome'].str.contains('provirus_')) & (group_df_votu['viral_genome'].str.startswith(snakemake.params.group + "_"))]
print(group_df_provirus)
group_proviruses_set = set(group_df_provirus['viral_genome'])

# parse untrimmed sequences and extract lengths
contig_lengths = {}

for record in SeqIO.parse(str(snakemake.input.derep_untrimmed), "fasta"):
    if record.id in group_proviruses_set:
        contig_lengths[record.id] = len(record.seq)

contig_length_df = pd.DataFrame(list(contig_lengths.items()), columns=['viral_genome', 'length'])

# merge prophage info with contig length
group_df_length = group_df_provirus.merge(contig_length_df, on='viral_genome')

# retain only the longest untrimmed contig from each votu
group_df_length_sort = group_df_length.sort_values('length', ascending=False)
group_df_length_top = group_df_length_sort.groupby('votu_representative', as_index=False).first()


# determine prophage coordinates
group_df_length_top['prophage_start'] = group_df_length_top.apply(lambda x: x.viral_genome.split('|checkv_provirus_')[1].split('_')[0] if "|checkv_provirus_" in str(x.viral_genome) else x.viral_genome.split('|provirus_')[1].split('_')[0], axis=1)
group_df_length_top['prophage_stop'] = group_df_length_top.apply(lambda x: x.viral_genome.split('|checkv_provirus_')[1].split('_')[1] if "|checkv_provirus_" in str(x.viral_genome) else x.viral_genome.split('|provirus_')[1].split('_')[1], axis=1)

# save report of integrated prophages identified within the group
group_df_length_top.to_csv(str(snakemake.output.report), index=False, sep='\t')

# create coordinates file for propagate
coord_file = group_df_length_top[['viral_genome', 'prophage_start', 'prophage_stop']]
coord_file['scaffold'] = coord_file['viral_genome'].str.split('|', expand=True)[1]
coord_file['fragment'] = coord_file['viral_genome']
coord_file.rename(columns={'prophage_start':'start', 'prophage_stop':'stop'}, inplace=True)
coord_file_final = coord_file[['scaffold', 'fragment', 'start', 'stop']]

coord_file_final.to_csv(str(snakemake.output.coords), index=False, sep="\t")

# retain only representative prophage sequences
rep_prophages_set = set(group_df_length_top['viral_genome'])
rep_prophages = []

for record in SeqIO.parse(str(snakemake.input.derep_untrimmed), "fasta"):
    if record.id in rep_prophages_set:
        record.description = record.id.split('|')[1]
        record.id = record.id.split('|')[1]
        rep_prophages.append(record)

SeqIO.write(rep_prophages, str(snakemake.output.sequences), "fasta")