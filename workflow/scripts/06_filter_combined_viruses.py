import pandas as pd
from Bio import SeqIO

contigs_to_virus_dict = {}

virus_seq_ids = set()
virus_seqs = []

untrimmed_hq_viruses = []

seq_num = 0
# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.viruses), "fasta"):
    seq_num += 1
    if seq_num%1000 == 0:
        print(seq_num)
    if record.id in virus_seq_ids:
        continue

    virus_seq_ids.add(record.id)
    virus_seqs.append(record)

    if '|checkv_' in record.id:
        provirus_base_checkv = record.id.split('|checkv')[0]
        untrimmed_hq_viruses.append(provirus_base_checkv)
        contigs_to_virus_dict[provirus_base_checkv] = record.id
    else:
        provirus_base_checkv = record.id

    if '|provirus' in record.id:
        provirus_base_genomad = record.id.split('|provirus')[0]
        untrimmed_hq_viruses.append(provirus_base_genomad)
        contigs_to_virus_dict[provirus_base_genomad] = record.id
    else:
        provirus_base_genomad = provirus_base_checkv

    if 'external|' in provirus_base_genomad:
        untrimmed_hq_viruses.append(provirus_base_genomad.split('external|')[1])
        contigs_to_virus_dict[provirus_base_genomad.split('external|')[1]] = record.id
    else:
        untrimmed_hq_viruses.append(provirus_base_genomad)
        contigs_to_virus_dict[provirus_base_genomad] = record.id

print('Writing untrimmed virus seqs')
SeqIO.write(virus_seqs, str(snakemake.output.viruses), "fasta")

untrimmed_virus_seq_ids = set()
untrimmed_hq_virus_seqs = []

seq_num = 0
# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.untrimmed), "fasta"):
    seq_num += 1
    if seq_num%1000 == 0:
        print(seq_num)
    if record.id in untrimmed_virus_seq_ids:
        continue
    untrimmed_virus_seq_ids.add(record.id)
    if 'external|' in record.id:
        record.id = record.id.split("external|")[1]
    if record.id in untrimmed_hq_viruses:
        record.id = contigs_to_virus_dict[record.id]
        untrimmed_hq_virus_seqs.append(record)

# save all sequences to specified file
SeqIO.write(untrimmed_hq_virus_seqs, str(snakemake.output.untrimmed), "fasta")