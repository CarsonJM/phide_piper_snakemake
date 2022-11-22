import os
from Bio import SeqIO

# create list to store contig names in
contig_names = []

# read in the fasta
for record in SeqIO.parse(str(snakemake.input), "fasta"):
    contig_names.append(record.id)

stb_df = pd.DataFrame()
stb_df['scaffold'] = contig_names
stb_df['bin'] = stb_df['scaffold']
stb_df.to_csv(str(snakemake.output), sep='\t', index=False, header=False)