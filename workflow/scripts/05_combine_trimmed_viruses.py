import pandas as pd
from Bio import SeqIO

combined = []

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.checkv_viruses), "fasta"):
    combined.append(record)

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.checkv_proteins), "fasta"):
    record.id = record.id.rpartition('_')[0]
    combined.append(record)

# save all sequences to specified file
SeqIO.write(combined, str(snakemake.output), "fasta")