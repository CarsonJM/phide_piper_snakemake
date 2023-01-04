import pandas as pd
from Bio import SeqIO

# open clustering results
clusters = pd.read_csv(str(snakemake.input.clusters), sep='\t', names=['cluster_rep', 'cluster_members'])

cluster_reps = set(clusters['cluster_rep'])

# extract representative sequences from fasta file
rep_sequences = []
for record in SeqIO.parse(str(snakemake.input.viruses), "fasta"):
    if '|' in record.id:
        if record.id.split('|')[1] in cluster_reps:
            rep_sequences.append(record)

# save all sequences to specified file
SeqIO.write(rep_sequences, str(snakemake.output), "fasta")