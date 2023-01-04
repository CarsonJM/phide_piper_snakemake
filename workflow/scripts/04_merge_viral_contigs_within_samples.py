from Bio import SeqIO

combined_sequences = []
# extract viral contigs
for record in SeqIO.parse(str(snakemake.input), "fasta"):
    if 'external|' in record.id:
        combined_sequences.append(record)
    else:
        record.id = snakemake.params.assembly + '|' + record.id
        combined_sequences.append(record)

SeqIO.write(combined_sequences, str(snakemake.output), "fasta")