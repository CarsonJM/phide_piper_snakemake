import pandas as pd
from Bio import SeqIO

# load imgvr_6 metadata
imgvr_6_meta = pd.read_csv('/gscratch/scrubbed/carsonjm/resources/imgvr_6/IMGVR_all_Sequence_information.tsv', sep='\t')
imgvr_6_meta_hq = imgvr_6_meta[imgvr_6_meta['Estimated completeness'] >= 90]

imgvr_6_meta_hq_set = set(imgvr_6_meta_hq['UVIG'])
imgvr_6_hq = []

# read in the fasta
for record in SeqIO.parse('/gscratch/scrubbed/carsonjm/resources/imgvr_6/IMGVR_all_nucleotides.fna', "fasta"):
    # retain contigs in hq list
    record.id = record.id.partition('|')[0]
    if record.id in imgvr_6_meta_hq_set:
        # add sequences passing filter to list 
        imgvr_6_hq.append(record)

# save saved sequences to specified file
SeqIO.write(imgvr_6_hq, '/gscratch/scrubbed/carsonjm/resources/imgvr_6/IMGVR_hq_nucleotides.fna', "fasta")