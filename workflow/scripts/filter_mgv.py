import pandas as pd
from Bio import SeqIO

# load mgv metadata
mgv_meta = pd.read_csv('/gscratch/scrubbed/carsonjm/resources/mgv_db/MGV_v1.0_2021_07_08/mgv_contig_info.tsv', sep='\t')
mgv_meta_hq = mgv_meta[mgv_meta['completeness'] >= 90]

mgv_meta_hq_set = set(mgv_meta_hq['contig_id'])
mgv_hq = []

# read in the fasta
for record in SeqIO.parse('/gscratch/scrubbed/carsonjm/resources/mgv_db/MGV_v1.0_2021_07_08/mgv_contigs.fna', "fasta"):
    # retain contigs in hq list
    if record.id in mgv_meta_hq_set:
        # add sequences passing filter to list 
        mgv_hq.append(record)

# save saved sequences to specified file
SeqIO.write(mgv_hq, '/gscratch/scrubbed/carsonjm/resources/mgv_db/MGV_v1.0_2021_07_08/mgv_hq_contigs.fna', "fasta")