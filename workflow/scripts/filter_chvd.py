import pandas as pd
from Bio import SeqIO

# load chvd metadata
chvd_meta = pd.read_excel('/gscratch/scrubbed/carsonjm/resources/chvd/HV3_table2_master_table.xlsx')
chvd_meta_hq = chvd_meta[chvd_meta['CheckV_estimated_completeness'] >= 90]

chvd_meta_hq_set = set(chvd_meta_hq['contig_name'])
chvd_hq = []

# read in the fasta
for record in SeqIO.parse('/gscratch/scrubbed/carsonjm/resources/chvd/CHVD_clustered_mash99_v1.fasta', "fasta"):
    # retain contigs in hq list
    record.id = record.id.partition('@')[0]
    if record.id in chvd_meta_hq_set:
        # add sequences passing filter to list 
        chvd_hq.append(record)

# save saved sequences to specified file
SeqIO.write(chvd_hq, '/gscratch/scrubbed/carsonjm/resources/chvd/CHVD_hq.fasta', "fasta")