import pandas as pd
from Bio import SeqIO

# open clustering results
clusters = open(str(snakemake.input.clusters), 'r')

derep_reps= []
for line in clusters:
    stripped = line.strip()
    centroid, nodes = stripped.split('\t')
    derep_base = centroid
    if len(derep_base.split('|provirus')) > 1:
        derep_base = derep_base.split('|provirus')[0]
    if len(derep_base.split('|checkv_provirus')) > 1:
        derep_base = derep_base.split('|checkv_provirus')[0]
    derep_reps.append(derep_base)

derep_reps_set = set(derep_reps)

# extract representative sequences from fasta file
derep_rep_sequences = []

for record in SeqIO.parse(str(snakemake.input.viruses), "fasta"):
    record_id_base = record.id
    if len(record_id_base.split('|provirus')) > 1:
        record_id_base = record.id.split('|provirus')[0]
    if len(record_id_base.split('|checkv_provirus')) > 1:
        record_id_base = record.id.split('|checkv_provirus')[0]
    if record_id_base in derep_reps_set:
        derep_rep_sequences.append(record)

# save all sequences to specified file
SeqIO.write(derep_rep_sequences, str(snakemake.output), "fasta")
