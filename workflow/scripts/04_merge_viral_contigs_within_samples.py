from Bio import SeqIO

combined_sequences = []
# extract viral contigs
if snakemake.params.run_genomad or snakemake.params.run_external:
    for record in SeqIO.parse(str(snakemake.input[0]), "fasta"):
        record.id = snakemake.params.assembly + '_' + record.id
        combined_sequences.append(record)

if snakemake.params.run_genomad and snakemake.params.run_external:
    for record in SeqIO.parse(str(snakemake.input[1]), "fasta"):
        record.id = snakemake.params.assembly + '_' + record.id
        combined_sequences.append(record)

# prophage = []

# for record in SeqIO.parse(str(snakemake.input.gn_prophage), "fasta"):
#     record.id = snakemake.params.assembly + '_' + record.id
#     prophage.append(record.id.partition('|')[0])
#     combined_sequences.append(record)

# for record in SeqIO.parse(str(snakemake.input.vb_prophage), "fasta"):
#     record.id = snakemake.params.assembly + '_' + record.id
#     name = (record.id.partition('_fragment')[0])
#     if name in prophage:
#         continue
#     else:
#         combined_sequences.append(record)

SeqIO.write(combined_sequences, str(snakemake.output), "fasta")