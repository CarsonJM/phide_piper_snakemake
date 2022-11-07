from Bio import SeqIO

contig = ""
protein_num = 1
updated_proteins = []

for record in SeqIO.parse(str(snakemake.input), 'fasta'):
    if record.id == contig:
        record.id = record.id + "_" + str(protein_num)
        protein_num += 1
        updated_proteins.append(record)
    else:
        contig = record.id
        protein_num = 1
        record.id = record.id + "_" + str(protein_num)
        protein_num += 1
        updated_proteins.append(record)

SeqIO.write(updated_proteins, str(snakemake.output), "fasta")