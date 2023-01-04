import pandas as pd
from Bio import SeqIO

contigs_to_virus_dict = {}

virus_seq_ids = []
virus_seqs = []

untrimmed_hq_viruses = []

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.viruses), "fasta"):
    if record.id in set(virus_seq_ids):
        continue

    virus_seq_ids.append(record.id)
    virus_seqs.append(record)

    if '|checkv_' in record.id and '|provirus' not in record.id:
        provirus_base = record.id.split('|')[1]
        untrimmed_hq_viruses.append(provirus_base.rpartition('_')[0])
        contigs_to_virus_dict[provirus_base.rpartition('_')[0]] = record.id
    else:
        untrimmed_hq_viruses.append(record.id.split('|')[1])
        contigs_to_virus_dict[record.id.split('|')[1]] = record.id

untrimmed_virus_seq_ids = []
untrimmed_hq_virus_seqs = []
# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.untrimmed), "fasta"):
    if record.id in set(untrimmed_virus_seq_ids):
        continue
    untrimmed_virus_seq_ids.append(record.id)
    if '|' in record.id:
        record.id = record.id.split("|")[1]
    if record.id in untrimmed_hq_viruses:
        record.id = contigs_to_virus_dict[record.id]
        untrimmed_hq_virus_seqs.append(record)

# save all sequences to specified file
SeqIO.write(virus_seqs, str(snakemake.output.viruses), "fasta")
SeqIO.write(untrimmed_hq_virus_seqs, str(snakemake.output.untrimmed), "fasta")