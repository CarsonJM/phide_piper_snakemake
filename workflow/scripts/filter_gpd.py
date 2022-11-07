import pandas as pd
from Bio import SeqIO

# load gpd metadata
gpd_meta = pd.read_csv('/gscratch/scrubbed/carsonjm/resources/gpd/GPD_metadata.tsv', sep='\t')
gpd_meta_hq = gpd_meta[gpd_meta['checkV_completion'] >= 90]

gpd_meta_hq_set = set(gpd_meta_hq['GPD_id'])
gpd_hq = []

# read in the fasta
for record in SeqIO.parse('/gscratch/scrubbed/carsonjm/resources/gpd/GPD_sequences.fa', "fasta"):
    # retain contigs in hq list
    if record.id in gpd_meta_hq_set:
        # add sequences passing filter to list 
        gpd_hq.append(record)

# save saved sequences to specified file
SeqIO.write(gpd_hq, '/gscratch/scrubbed/carsonjm/resources/gpd/GPD_hq_sequences.fa', "fasta")