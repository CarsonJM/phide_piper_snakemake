import pandas as pd
from Bio import SeqIO

# open clustering results
clusters = open(str(snakemake.input.clusters), 'r')

votu_reps= []
for line in clusters:
    stripped = line.strip()
    centroid, nodes = stripped.split('\t')
    votu_reps.append(centroid)

votu_reps_set = set(votu_reps)

# extract representative sequences from fasta file
votu_rep_sequences = []
votu_rep_ids = []
rep = 1
for record in SeqIO.parse(str(snakemake.input.viruses), "fasta"):
    if record.id in votu_reps_set:
        if record.id in set(votu_rep_ids):
            record.id = record.id + "_" + str(rep)
            rep += 1
        votu_rep_ids.append(record.id)
        votu_rep_sequences.append(record)

votu_reps_set = set(votu_reps)

# extract representative sequences from fasta file
votu_rep_sequences_untrimmed = []
votu_rep_ids_untrimmed = []
rep = 1
for record in SeqIO.parse(str(snakemake.input.untrimmed_viruses), "fasta"):
    if record.id in votu_reps_set:
        if record.id in set(votu_rep_ids_untrimmed):
            record.id = record.id + "_" + str(rep)
            rep += 1
        votu_rep_ids_untrimmed.append(record.id)
        votu_rep_sequences_untrimmed.append(record)


# extract representative sequences from fasta file
votu_rep_prot_ids= []
votu_rep_prots = []
rep = 1
for record in SeqIO.parse(str(snakemake.input.proteins), "fasta"):
    contig_id = record.id.rpartition('_')[0]
    if contig_id in votu_reps_set:
        if contig_id in set(votu_rep_prot_ids):
            record.id = record.id + "_" + str(rep)
            rep += 1
        votu_rep_prot_ids.append(record.id)
        votu_rep_prots.append(record)

# save all sequences to specified file
SeqIO.write(votu_rep_sequences, str(snakemake.output.viruses), "fasta")
SeqIO.write(votu_rep_sequences_untrimmed, str(snakemake.output.untrimmed_viruses), "fasta")
SeqIO.write(votu_rep_prots, str(snakemake.output.proteins), "fasta")