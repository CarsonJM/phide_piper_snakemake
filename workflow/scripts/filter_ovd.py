import pandas as pd
from Bio import SeqIO

# load ovd metadata
ovd_meta = pd.read_excel('/gscratch/scrubbed/carsonjm/resources/oral_virome_db/1-s2.0-S2589004222006897-mmc4.xlsx', header=3)
ovd_meta_hq = ovd_meta[ovd_meta['Completeness'] >= 90]

ovd_meta_hq_set = set(ovd_meta_hq['Unnamed: 0'])
ovd_hq = []

# read in the fasta
for record in SeqIO.parse('/gscratch/scrubbed/carsonjm/resources/oral_virome_db/OVD-genomes.fa', "fasta"):
    # retain contigs in hq list
    if record.id in ovd_meta_hq_set:
        # add sequences passing filter to list 
        ovd_hq.append(record)

# save saved sequences to specified file
SeqIO.write(ovd_hq, '/gscratch/scrubbed/carsonjm/resources/oral_virome_db/OVD_hq_genomes.fa', "fasta")