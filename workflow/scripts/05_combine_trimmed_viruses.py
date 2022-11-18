import pandas as pd
from Bio import SeqIO

combined = []

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.checkv_viruses), "fasta"):
    combined.append(record)

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.checkv_proviruses), "fasta"):
    combined.append(record)

# save all sequences to specified file
SeqIO.write(combined, str(snakemake.output), "fasta")