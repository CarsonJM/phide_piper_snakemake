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

# parse through and combine virus sequences for each sample
for record in SeqIO.parse(str(snakemake.input.checkv_viruses), "fasta"):
    if record.id in hq_viruses:
        hq_virus_seqs.append(record)

for record in SeqIO.parse(str(snakemake.input.checkv_proviruses), "fasta"):
    if record.id.rpartition('_')[0] in hq_viruses:
        if '|' in record.id:
            genomad_provirus = record.id.split('|')[1]
            genomad_start = genomad_provirus.split('_')[1]
            checkv_provirus = record.description.split(' ')[1]
            checkv_provirus_coords = checkv_provirus.split('/')[0]
            checkv_start, checkv_stop = checkv_provirus_coords.split('-')
            checkv_start_total = int(checkv_start) + int(genomad_start)
            checkv_stop_total = int(checkv_stop) + int(genomad_start)
            record.id = record.id + "|checkv_provirus_" + str(checkv_start_total) + "_" + str(checkv_stop_total)
        else:
            checkv_provirus = record.description.split(' ')[1]
            checkv_provirus_coords = checkv_provirus.split('/')[0]
            checkv_start, checkv_stop = checkv_provirus_coords.split('-')
            checkv_start_total = int(checkv_start) + int(genomad_start)
            checkv_stop_total = int(checkv_stop) + int(genomad_start)
            record.id = record.id + "|checkv_provirus_" + str(checkv_start_total) + "_" + str(checkv_stop_total)

        hq_virus_seqs.append(record)


# save all sequences to specified file
SeqIO.write(hq_virus_seqs, str(snakemake.output), "fasta")