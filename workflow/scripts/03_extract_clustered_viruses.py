import pandas as pd
from Bio import SeqIO

# open clustering results
clusters = open(str(snakemake.input.clusters), 'r')

votu_reps= []
for line in clusters:
    stripped = line.strip()
    centroid, nodes = stripped.split('\t')
    votu_reps.append(centroid)

# extract representative sequences from fasta file
votu_rep_sequences = []
for record in SeqIO.parse(str(snakemake.input.viruses), "fasta"):
    if record.id in votu_reps:
        votu_rep_sequences.append(record)
    else:
        continue

# save all sequences to specified file
SeqIO.write(votu_rep_sequences, str(snakemake.output), "fasta")