import os
from Bio import SeqIO

if not str(snakemake.input).endswith(".fasta"):
    raise TypeError("Input contigs file must be in fasta format")
if not str(snakemake.output).endswith(".fasta"):
    raise TypeError("Output contigs file must be in fasta format")
if not os.path.isfile(str(snakemake.input)):
    raise TypeError("Input contigs file does not exist")

# create list to store filtered contigs
long_contigs = []

# read in the fasta
for record in SeqIO.parse(str(snakemake.input), "fasta"):
    # remove contigs shorter than min_len
    if len(record.seq) >= snakemake.params.min_length:
        # add sequences passing filter to list
        long_contigs.append(record)

# save saved sequences to specified file
SeqIO.write(long_contigs, str(snakemake.output), "fasta")