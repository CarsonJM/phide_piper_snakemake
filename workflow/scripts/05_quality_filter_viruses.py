import pandas as pd
from Bio import SeqIO

# load checkv results
checkv_df = pd.read_csv(str(snakemake.input.checkv_results), sep="\t")

# filter checkv results based on input:
if snakemake.params.min_completeness == "":
    checkv_filtered = checkv_df[(checkv_df["viral_genes"] >= snakemake.params.min_viral_genes)
                                & (checkv_df["host_genes"] <= snakemake.params.max_bacterial_genes)]

elif snakemake.params.min_completeness != "":
    checkv_filtered = checkv_df[(checkv_df["completeness"] >= snakemake.params.min_completeness)
                                & (checkv_df["viral_genes"] >= snakemake.params.min_viral_genes)
                                & (checkv_df["host_genes"] <= snakemake.params.max_bacterial_genes)]

if snakemake.params.remove_proviruses:
    checkv_filtered = checkv_filtered[(checkv_filtered["provirus"] != 'Yes')]

hq_viruses = set(checkv_filtered["contig_id"])
hq_virus_seqs = []
hq_virus_prots = []
hq_untrimmed_seqs = []

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.checkv_viruses), "fasta"):
    if record.id in hq_viruses:
        hq_virus_seqs.append(record)

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.checkv_proteins), "fasta"):
    if record.id.rpartition('_')[0] in hq_viruses:
        hq_virus_prots.append(record)


for record in SeqIO.parse(str(snakemake.input.untrimmed_viruses), "fasta"):
    if record.id in hq_viruses:
        hq_untrimmed_seqs.append(record)

# save all sequences to specified file
SeqIO.write(hq_virus_seqs, str(snakemake.output.viruses), "fasta")
SeqIO.write(hq_virus_prots, str(snakemake.output.proteins), "fasta")
SeqIO.write(hq_untrimmed_seqs, str(snakemake.output.untrimmed_viruses), "fasta")