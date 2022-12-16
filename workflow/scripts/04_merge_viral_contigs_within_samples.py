from Bio import SeqIO

combined_sequences = []
# extract viral contigs
for record in SeqIO.parse(str(snakemake.input), "fasta"):
    record.id = snakemake.params.assembly + '|' + record.id
    combined_sequences.append(record)

SeqIO.write(combined_sequences, str(snakemake.output), "fasta")