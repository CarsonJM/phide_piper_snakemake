import pandas as pd
from Bio import SeqIO

# open clustering results
read_screen =  pd.read_csv(str(snakemake.input.read_screen), sep='\t', header=None, index_col=False,
                    names=['identity', 'shared-hashes', 'median-multiplicity', 'p-value', 'query-id', 'query-comment'])
read_screen['shared-hashes'] = read_screen['shared-hashes'].str.partition('/')[0]
read_screen_filt = read_screen[((read_screen['identity']) > snakemake.params.min_mash_score ) & (read_screen["median-multiplicity"] > snakemake.params.min_mash_multiplicity)]

read_screen_filt_queries = set(read_screen_filt['query-id'])

# extract representative sequences from fasta file
virusdb_hits = []

for record in SeqIO.parse(str(snakemake.input.virusdb), "fasta"):
    if record.id in read_screen_filt_queries:
        virusdb_hits.append(record)
        
# save all sequences to specified file
SeqIO.write(virusdb_hits, str(snakemake.output), "fasta")