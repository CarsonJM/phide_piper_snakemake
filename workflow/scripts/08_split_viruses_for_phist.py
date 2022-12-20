import pandas as pd
from Bio import SeqIO

filtered_viruses = []

# parse through and combine provirus sequences for each sample
for record in SeqIO.parse(str(snakemake.input), "fasta"):

    if 'chunk' in record.id:
        # save all sequences to specified file
        SeqIO.write(record, str(snakemake.params.fasta_dir) + record.id + ".fna", "fasta")

output = open(str(snakemake.output), 'x')
output.close()